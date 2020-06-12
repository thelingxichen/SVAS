#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Usage:
    complexsv spectrum --sample=STR --tool=STR \
        --sv_fn=IN_FILE [--cnv_fn=IN_FILE] --out_dir=DIR \
        [--tran_psl_fn=IN_FILE] [--gene_list=IN_FILE] \
        [--enhancer_fn=IN_FILE] [--super_enhancer_fn=IN_FILE] [--enhancer_source=STR]
    complexsv call --sample=STR --tool=STR \
        --bam_fn=IN_FILE \
        --sv_fn=IN_FILE [--cnv_fn=IN_FILE] --out_dir=DIR \
        [--chromthripsis_region_fn=IN_FILE] \
        [--tran_psl_fn=IN_FILE] [--gene_list=IN_FILE] \
        [--enhancer_fn=IN_FILE] [--super_enhancer_fn=IN_FILE] [--enhancer_source=STR]
    complexsv draw --sample=STR --tool=STR \
        --out_dir=DIR \
        [--tran_psl_fn=IN_FILE] [--gene_list=IN_FILE] \
        [--enhancer_fn=IN_FILE] [--super_enhancer_fn=IN_FILE] [--enhancer_source=STR]
    complexsv -h | --help

Options:
    OPTIONS EXPAINATION

    -v | --verbose	Show log information on console.
    -h | --help		Show this screen.

Others:
    tenxtools.complexsv
    ~~~~~~~~~~~~~~~~~~~~~~~

    @Input: DESCRIPTION
    @Return:

    @Copyright: (c) 2020-05 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
import os
import importlib

from tenxtools.mydocopt import docopt
from tenxtools.utils import module
from tenxtools.utils import myio
from tenxtools.utils import constants
from tenxtools.utils import psl_tran
from tenxtools.annotation import enhancer
from tenxtools.sv import base as sv_base
from csv import base
from csv import bplot_complexsv
import sv


class ComplexSVRegionGroupGenerator(module.Module):

    def __init__(self, reader_cls=None, writer_cls=myio.Writer, groups=None,
                 bam_fn=None,
                 tran_psl_fn=None, gene_list=None,
                 enhancer_fn=None, super_enhancer_fn=None, enhancer_source=None,
                 *args, **kwargs):
        self.groups = groups or []
        self.bam_fn = bam_fn
        if tran_psl_fn:
            self.tran_psl_reader = psl_tran.Reader(tran_psl_fn, gene_list)
        else:
            self.tran_psl_reader = None

        for group in self.groups:
            group.tran_psl_reader = self.tran_psl_reader

        if enhancer_fn:
            self.enhancer_reader = enhancer.EnhancerAltasReader(in_fn=enhancer_fn)
        else:
            self.enhancer_reader = None
        if super_enhancer_fn:
            self.super_enhancer_reader = enhancer.dbSUPERReader(in_fn=super_enhancer_fn)
        else:
            self.super_enhancer_reader = None
        self.enhancer_source = enhancer_source

        super(ComplexSVRegionGroupGenerator, self).__init__(reader_cls=reader_cls, writer_cls=writer_cls, *args, **kwargs)

    def __iter__(self):
        if issubclass(self.reader_cls, sv_base.SVReader):
            return self.call()

        if issubclass(self.reader_cls, base.RegionGroupReader):
            return self.read()

    def call(self):
        self._group_sv()
        cnt_dict = {}
        for i, group in enumerate(self.groups):
            self._collect_enhancer(group)

            # filter
            group_type = group.group_type
            '''
            if group.enhancer_num + group.super_enhancer_num + group.gene_num < 1:
                continue
            '''

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
            self._write_group(group)
            yield group.group_type, group

    def _write_group(self, group):
        out_dir = os.path.join(self.out_dir, group.group_type, group.group_name)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        sv_ids = set()
        out_prefix = os.path.join(out_dir, '{}.{}'.format(self.sample, group.group_name))
        with open(out_prefix + '.sv', 'w') as sv_f, \
            open(out_prefix + '.region', 'w') as region_f, \
            open(out_prefix + '.enhancer', 'w') as enhancer_f, \
            open(out_prefix + '.super_enhancer', 'w') as super_enhancer_f:
            for region in group.region_list:
                region_f.write('{}:{}-{}\n'.format(region.chrom, region.start, region.end))
                # write sv
                for sv in region.sv_list:
                    if sv.id not in sv_ids:
                        sv_f.write('{}\n'.format(sv))
                        sv_ids.add(sv.id)
                # write enhancer
                for meta, e in region.enhancer_list:
                    enhancer_f.write(str(e) + '\n')
                for meta, e in region.super_enhancer_list:
                    super_enhancer_f.write(str(e) + '\n')

    def _assign_enhancer_list(self, enhancer_fn, region):
        if not os.path.isfile(enhancer_fn):
            return
        region.enhancer_list = list(enhancer.EnhancerAltasReader(in_fn=enhancer_fn))

    def _assign_super_enhancer_list(self, super_enhancer_fn, region):
        if not os.path.isfile(super_enhancer_fn):
            return
        region.super_enhancer_list = list(enhancer.dbSUPERReader(in_fn=super_enhancer_fn))

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

    def _collect_enhancer(self, group):
        def condition_func(x):
            return self.enhancer_source in x.source
        for region in group.region_list:
            if self.enhancer_reader:
                region.enhancer_list = list(self.enhancer_reader.fetch(region.chrom, region.start, region.end, condition_func))

            if self.super_enhancer_reader:
                region.super_enhancer_list = list(self.super_enhancer_reader.fetch(region.chrom, region.start, region.end, condition_func))

    def _group_sv(self, linkage_distance=1000000.0):  # 1M
        for i, sv in enumerate(self.in_data):
            if sv.chrom_5p not in constants.chrs or sv.chrom_3p not in constants.chrs:
                continue
            if sv.chrom_5p == sv.chrom_3p and \
                    abs(sv.bkpos_5p - sv.bkpos_3p) < 1000:
                continue
            # print(sv.chrom_5p, sv.bkpos_5p/linkage_distance, sv.chrom_3p, sv.bkpos_3p/linkage_distance)
            self._append_sv_to_groups(sv)

    def _append_sv_to_groups(self, sv):
        merge_groups = [group for group in self.groups if True in group.has_sv(sv)[0]]
        groups = [group for group in self.groups if True not in group.has_sv(sv)[0]]
        if merge_groups:
            group = merge_groups[0]
            group.append_sv(sv)
            for g in merge_groups[1:]:
                g.append_sv(sv)
                group.merge(g)
            groups.append(group)
        else:  # new group
            group = base.RegionGroup(tran_psl_reader=self.tran_psl_reader)
            group.append_sv(sv)
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


def run_spectrum(**args):
    regions = [base.Region(chrom=chrom, start=0, end=constants.hg19_fai_bp[chrom]) for chrom in constants.chrs]
    run_call(regions=regions, **args)


def run_call(chromthripsis_region_fn=None, group=None, sv_module=None, regions=None, **args):
    region = base.Region(chrom='chr7', start=13500000, end=15000000)
    groups = [ base.RegionGroup(region_list=[region]),
              ]
    groups = []

    sv_records = sv.read_vcf(args['sv_fn'])

    reader_cls = sv_module.SVReader
    groups = ComplexSVRegionGroupGenerator(
        reader_cls=reader_cls,
        groups=groups,
        in_data=sv_records,
        write=False,  # True,
        suffix='bkp_region',
        **args
    ).evaluate()
    groups = list(groups)
    data = {}
    for meta, group in groups:
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


def run_draw(**args):
    reader_cls = base.RegionGroupReader
    groups = ComplexSVRegionGroupGenerator(
        reader_cls=reader_cls,
        write=False,  # True,
        **args
    ).evaluate()

    groups = list(groups)
    data = {}
    for meta, group in groups:
        if meta in data:
            data[meta].append((meta, group))
        else:
            data[meta] = [(meta, group)]

    for meta, group_list in data.items():
        out_fn = os.path.join(args['out_dir'], '{}.{}.region.svg'.format(args['sample'], meta))
        bplot_complexsv.run(draw=True,
                            sample=args['sample'],
                            groups=group_list,
                            out_fn=out_fn)
    return groups, None


def run(spectrum=None, call=None, draw=None, tool=None, **args):
    sv_module = importlib.import_module('tenxtools.sv.{}'.format(tool))
    if spectrum:
        run_spectrum(sv_module=sv_module, **args)
    if call:
        run_call(sv_module=sv_module, **args)
    if draw:
        run_draw(sv_module=sv_module, **args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    run(**args)
