#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    SVAS.complexsv
    ~~~~~~~~~~~~~~~~~~~~~~~

    @Input: DESCRIPTION
    @Return:

    @Copyright: (c) 2020-05 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.

Usage:
    complexsv.py call --sv_fn=IN_FILE --out_dir=OUT_DIR [--sample=STR]
    complexsv.py -h | --help

Options:
    -h --help           Show this screen.
    --version           Show version.
    --sv_fn=IN_FILE     Path of SvABA sv vcf file.
    --out_dir=OUT_DIR   Path of output directory.
    --sample            Sample name.
"""

import os
import docopt

from biotool import constants

import sv

from compsv import base
from compsv import bplot_complexsv


class ComplexSVRegionGroupGenerator():

    def __init__(self, sv_records=None, groups=None, bam_fn=None, 
        id2genes=None, *args, **kwargs):
        self.sv_records = sv_records or []
        self.groups = groups or []
        self.bam_fn = bam_fn
        self.id2genes = id2genes or {}

    def call(self):
        self._group_sv()
        cnt_dict = {}
        for i, group in enumerate(self.groups):

            # filter
            group_type = group.group_type

            if group_type in cnt_dict:
                cnt_dict[group_type] += 1
            else:
                cnt_dict[group_type] = 1

            group.group_name = '{}{}'.format(group_type, cnt_dict[group_type])
            group.expand_region()

            # linkage information
            length = 0
            for region in group.region_list:
                length += region.end - region.start

            '''
            if length <= 10000:
                linkage_heatmap.run_call(bam_fn=self.bam_fn,
                                         regions=group.region_list,
                                         out_dir=out_dir)
            '''
            yield group

    def _assign_cn_list(self, cn_fn, region):
        if not os.path.isfile(cn_fn):
            return
        for line in open(cn_fn, 'r'):
            region_str, cn = line.strip().split('\t')
            cn = float(cn)
            chrom = region_str.split(':')[0]
            start = int(region_str.split(':')[1].split('-')[0])
            end = int(region_str.split(':')[1].split('-')[1])

            if region.chrom == chrom and start >= region.start and end <= region.end and end - start > 1:
                region.cn_list.append(base.Region(chrom=chrom, start=start, end=end, cn=cn))

    def read(self):
        for group_type in os.listdir(os.path.join(self.out_dir)):
            if not os.path.isdir(group_type):
                continue
            for group_id in os.listdir(os.path.join(self.out_dir, group_type)):
                in_dir = os.path.join(self.out_dir, group_type, group_id)
                prefix_fn = os.path.join(in_dir, '{}.{}'.format(self.sample, group_id))
                sv_fn = prefix_fn + '.sv'
                cn_fn = prefix_fn + '.cn'
                enhancer_fn = prefix_fn + '.enhancer'
                super_enhancer_fn = prefix_fn + '.super_enhancer'
                region_fn = os.path.join(in_dir, '{}.{}.region'.format(self.sample, group_id))
                region_list = []
                for line in open(region_fn, 'r'):
                    chrom = line.split(':')[0]
                    start = int(line.split(':')[1].split('-')[0])
                    end = int(line.split(':')[1].split('-')[1])
                    region = base.Region(chrom=chrom, start=start, end=end, sv_fn=sv_fn, sample=self.sample)
                    self._assign_cn_list(cn_fn, region)
                    self._assign_enhancer_list(enhancer_fn, region)
                    self._assign_super_enhancer_list(super_enhancer_fn, region)
                    region_list.append(region)

                group = base.RegionGroup(region_list=region_list, group_name=group_id)
                yield group_type, group

    def _group_sv(self, minimal_distance=1000):  # 1k
        for i, record in enumerate(self.sv_records):
            if record.chrom_5p not in constants.chrs or record.chrom_3p not in constants.chrs:
                continue
            if record.chrom_5p == record.chrom_3p and \
                    abs(record.bkpos_5p - record.bkpos_3p) < minimal_distance:
                continue
            genes = list(self.id2genes.get(record.id[:-2], []))
            record.genes = genes
            self._append_sv_to_groups(record)

    def _append_sv_to_groups(self, record):  # code review
        merge_groups = [group for group in self.groups if True in group.has_sv(record)[0]]
        groups = [group for group in self.groups if True not in group.has_sv(record)[0]]
        if merge_groups:
            group = merge_groups[0]
            group.append_sv(record)
            for g in merge_groups[1:]:
                g.append_sv(record)
                group.merge(g)
            groups.append(group)
        else:  # new group
            group = base.RegionGroup()
            group.append_sv(record)
            groups.append(group)
        self.groups = groups

    def _bkp_in_regions(self, bkp):
        bs = [bkp.inside_region(region) for region in self.regions]
        return True in bs, bs

    def _sv_in_regions(self, record):
        bkp_5p = base.Breakpoint(record.chrom_5p, record.bkpos_5p,
                                 record.bkpos_5p_orientation, record.orientation)
        bkp_3p = base.Breakpoint(record.chrom_3p, record.bkpos_3p,
                                 record.bkpos_3p_orientation, record.orientation)
        _, bs_5p = self._bkp_in_regions(bkp_5p)
        _, bs_3p = self._bkp_in_regions(bkp_3p)
        bs = [b1 or b2 for b1, b2 in zip(bs_5p, bs_3p)]
        return True in bs, bs

    def _append_sv_to_regions(self, record, bs):
        if True in bs:
            for i, b in enumerate(bs):
                if b:
                    self.regions[i].append_sv(record, self.tran_psl_reader)


def run_call(**args):
    regions = [base.Region(chrom=chrom, start=0, end=constants.hg19_fai_bp[chrom]) for chrom in constants.chrs]
    groups = [base.RegionGroup(region_list=regions)]
    groups = []
    sv_records = sv.read_vcf(args['sv_fn'], **args)

    groups = ComplexSVRegionGroupGenerator(
        groups=groups,
        sv_records=sv_records,
        **args
    ).call()
    groups = list(groups)

    write_groups(groups, args['out_dir'], args['sample'])

    data = {}
    for group in groups:
        meta = group.group_type
        if meta in data:
            data[meta].append((meta, group))
        else:
            data[meta] = [(meta, group)]

    for meta, group_list in data.items():
        out_prefix = os.path.join(args['out_dir'], '{}.{}.region'.format(args['sample'], meta))
        bplot_complexsv.run(draw=True,
                            groups=group_list,
                            out_prefix=out_prefix,
                            **args)
    return groups, None


def write_groups(groups, out_dir, sample):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    sv_ids = set()
    out_prefix = os.path.join(out_dir, '{}'.format(sample))
    print('=====', out_prefix)
    with open(out_prefix + '.sv.txt', 'w') as sv_f, \
        open(out_prefix + '.region.txt', 'w') as region_f:
        header = ['chrom_5p', 'bkpos_5p', 'strand_5p', 'chrom_3p', 'bkpos_3p', 'strand_3p', 'junc_reads', 'group', 'type', 'clone', 'genes', 'sample']
        sv_f.write('{}\n'.format('\t'.join(header)))

        for group in groups:
            for region in group.region_list:
                region_f.write('{}\t{}:{}-{}\n'.format(group.group_name, region.chrom, region.start, region.end))
                # write sv
                for sv_re in region.sv_list:
                    if sv_re.id not in sv_ids:
                        clone_prev = sv_re.meta_info.get('CLONE_PREV', {})
                        clone_prev = ','.join([ '{}={}'.format(c,f) for c, f in clone_prev.items()])
                        sv_str = '\t'.join(map(str, [
                            sv_re.chrom_5p, sv_re.bkpos_5p, sv_re.strand_5p,
                            sv_re.chrom_3p, sv_re.bkpos_3p, sv_re.strand_3p,
                            sv_re.junc_reads, group.group_name, sv_re.sv_type2, clone_prev,
                            ','.join(sv_re.genes), sample
                            ]))
                        sv_f.write('{}\n'.format(sv_str))
                        sv_ids.add(sv_re.id)


def run(call=None, **args):
    if call:
        run_call(**args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    new_args = {}
    for k, v in args.items():
        new_args[k.replace('--', '')] = v
    run(**new_args)
