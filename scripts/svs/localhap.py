#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Usage:
    localhap call --sample=STR --tool=STR \
        --tumor_purity=FLOAT --pure_tumor_depth=FLOAT \
        --out_dir=DIR \
        --cn_fn=IN_FILE \
        [--tran_psl_fn=IN_FILE] [--gene_list=IN_FILE] \
        [--enhancer_fn=IN_FILE] [--super_enhancer_fn=IN_FILE] [--enhancer_source=STR]
    localhap -h | --help

Options:
    OPTIONS EXPAINATION

    -v | --verbose	Show log information on console.
    -h | --help		Show this screen.

Others:
    tenxtools.localhap
    ~~~~~~~~~~~~~~~~~~

    @Input: DESCRIPTION
    @Return:

    @Copyright: (c) 2018-10 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
import os
import copy

from tenxtools.mydocopt import docopt
from tenxtools.complexsv import base
from tenxtools.complexsv import complexsv
from tenxtools.cnv import cnv
from tenxtools.bplot import localhap as bplot_localhap

from LocalHaplotypeSolver import LocalHap
from LocalHaplotypeSolver.JuncGraph import Graph
from LocalHaplotypeSolver import SysSettings


def read_graph_from_region_group(region_group,
                                 sample=None, tumor_purity=None, pure_tumor_depth=None,
                                 cn_fn=None, **args):

    cnv_reader = cnv.CNVReader(in_fn=cn_fn)

    graph = Graph.Graph()
    graph.mName = sample
    graph.mAvgCov = 2.0
    graph.mAvgCovRaw = 2.0
    graph.mAvgPloidy = 2.0

    graph.mPurity = tumor_purity

    # graph_groups
    assign_graph_group_type(region_group)
    graph.mGroupPriority = region_group.region_priority

    for region in region_group.region_list:
        cnv_records = list(cnv_reader.fetch(region.chrom, region.start, region.end))
        # full, total, minor
        # print(cnv_records[0][1].fullCN, cnv_records[0][1].Cn, cnv_records[0][1].mCn)
        if not cnv_records:
            ploidy = [1, 1]  # ???
        else:
            source_cnv = cnv_records[0][1]
            sink_cnv = cnv_records[-1][1]
            if source_cnv.Cn > sink_cnv.Cn:
                cnv_record = sink_cnv
            else:
                cnv_record = source_cnv
            ploidy = [cnv_record.Cn-cnv_record.mCn, cnv_record.mCn]  # ???
        graph.mPloidy.update({region.graph_group_type: ploidy})
        graph.mPloidyRaw.update({region.graph_group_type: ploidy})
        # segment
        for i, seg in enumerate(region.cn_list):
            seg_type = region.graph_group_type
            seg_cov = seg.cn
            seg_pos = '{}:{}-{}'.format(seg.chrom, seg.start, seg.end)
            seg_cred = 1.0  # need to imporve
            seg_id = i + 1
            s = graph.addSegment(seg_type, seg_id, seg_cov, cred=seg_cred, position=seg_pos)
            s.mLowerBoundLimit = True
            s.mGroup = region.graph_group_type

            if i == 0:  # source
                graph.mSource.update({s.mGroup: s})
            if i == len(region.cn_list) - 1:
                graph.mSink.update({s.mGroup: s})

    sv_ids = set({})
    for region in region_group.region_list:
        # junction
        for i, sv in enumerate(region.sv_list):
            if sv.id in sv_ids:
                continue
            jun_cov = (sv.junc_reads + sv.span_reads) / float(pure_tumor_depth)
            j = graph.addJunction(*(list(get_vertex_info(sv, region_group)) + [jun_cov]))
            j.aJuncId = i + 1
            j.preferred = False
            j.cred = 1.0  # ???
            j.mLowerBoundLimit = True
            j.join_type = sv.orientation

            sv_ids.add(sv.id)

    SysSettings.HOST_SEG_LP_COE = float(1)
    SysSettings.VIRUS_SEG_LP_COE = float(1)
    SysSettings.JUNCTION_LP_COE = float(1)
    SysSettings.INFER_COV_COE = float(1)
    SysSettings.INFERRED_JUNCTION_WEIGHT = float(0)

    # graph._adjustCovByPurity() ???
    return graph


def get_vertex_info(sv, region_group):
    bs, bs_5p, bs_3p = region_group.has_sv(sv)
    is_5p = [i for i, b in enumerate(bs_5p) if b]
    is_3p = [i for i, b in enumerate(bs_3p) if b]
    if len(is_5p) > 1 or len(is_3p) > 1:
        raise Exception('is_5p > 1or is_3p > 1')
    region_5p = region_group.region_list[is_5p[0]]
    region_3p = region_group.region_list[is_3p[0]]

    seg_5p = []
    for i, seg in enumerate(region_5p.cn_list):
        start_offset = abs(seg.start - sv.bkpos_5p)
        end_offset = abs(sv.bkpos_5p - seg.end)
        offset = min(start_offset, end_offset)
        seg_5p.append((offset, i))
    seg_5p = sorted(seg_5p, key=lambda p: p[0])
    seg_5p = sorted(seg_5p[:2], key=lambda p: p[1])
    if sv.strand_5p == '+':
        v_id_5p = seg_5p[0][1] + 1
    else:
        v_id_5p = seg_5p[1][1] + 1

    seg_3p = []
    for i, seg in enumerate(region_3p.cn_list):
        start_offset = abs(seg.start - sv.bkpos_3p)
        end_offset = abs(sv.bkpos_3p - seg.end)
        offset = min(start_offset, end_offset)
        seg_3p.append((offset, i))
    seg_3p = sorted(seg_3p, key=lambda p: p[0])
    seg_3p = sorted(seg_3p[:2], key=lambda p: p[1])
    if sv.strand_3p == '+':
        v_id_3p = seg_3p[1][1] + 1
    else:
        v_id_3p = seg_3p[0][1] + 1

    v_type_5p = region_5p.graph_group_type
    v_type_3p = region_3p.graph_group_type

    v_dir_5p = sv.strand_5p
    v_dir_3p = sv.strand_3p
    return v_type_5p, v_id_5p, v_dir_5p, v_type_3p, v_id_3p, v_dir_3p


def assign_graph_group_type(region_group):
    # support AlphaBeta order
    region_group.region_priority = []
    for i, region in enumerate(region_group.region_list):
        region.graph_group_type = chr(ord('A') + i)
        region_group.region_priority.append(region.graph_group_type)


def run_call(**args):

    reader_cls = base.RegionGroupReader
    region_groups = complexsv.ComplexSVRegionGroupGenerator(
        reader_cls=reader_cls,
        write=False,  # True,
        **args
    ).evaluate()

    cnt = 0
    graphs = []
    for meta, region_group in region_groups:
        if meta == 'group':
            cnt += 1
            if cnt in [2, 3, 8, 11, 12]:
            # if cnt != 1:
                continue
            print(cnt)
            out_fn = os.path.join(args['out_dir'], meta,
                                  '{}{}'.format(meta, cnt),
                                  '{}.{}{}.ucyc.txt'.format(args['sample'], meta, cnt))

            graph = read_graph_from_region_group(region_group, **args)
            graph.printGraphInfo()
            localhap = LocalHap.LocalHaplotypeSolver(graph=graph)
            localhap.create_graph()
            localhap.balance_graph()
            g = copy.deepcopy(localhap.g)
            graphs.append((meta, g))

            try:
                localhap.get_cycles(output_path=out_fn, **args)
            except Exception:
                pass

    out_fn = os.path.join(args['out_dir'], '{}.{}.graph.svg'.format(args['sample'], meta))
    bplot_localhap.run(draw=True,
                       graphs=graphs,
                       out_fn=out_fn,
                       **args)


def run(call=None, **args):
    if call:
        run_call(**args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    run(**args)
