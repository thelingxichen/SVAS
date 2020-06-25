#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    complexsv.base
    ~~~~~~~~~~~~~~~~~~~


    @Copyright: (c) 2018-08 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""


# from tenxtools.sv import base as sv_base

from biotool import constants
from biotool import myio

from compsv import chromothripsis


class RegionGroup(myio.Record):

    fields = ['region_list']

    def __init__(self, **kwargs):
        super(RegionGroup, self).__init__(**kwargs)

    def set(self, group_name=None, region_list=None, linkage_distance=None,
            tran_psl_reader=None, **kargs):
        self.group_name = group_name
        self.region_list = region_list or []
        self.linkage_distance = linkage_distance or 1000000.0
        self.tran_psl_reader = tran_psl_reader

    def append_region(self, region):
        region_list = self.region_list + [region]
        region_list = sorted(region_list)
        prev = None
        res = []
        for region in region_list:
            if not prev:
                prev = region
                continue

            if prev.overlap(region):
                prev.merge(region)
            else:
                res.append(prev)
                prev = region
        res.append(prev)

        self.region_list = res

    def merge(self, other):
        for region in other.region_list:
            self.append_region(region)

    def expand_region(self):
        for region in self.region_list:
            region.expand()

    @property
    def group_type(self):
        if self.sv_num <= 1:
            if self.region_num <= 1:
                return 'single'
            else:
                if str(self) == '[chr4:128.568158-128.57017.0.002012M, chr15:75.880429-75.882593.0.002164M]':
                    print(self)
                    print(self.sv_num)
                    print([region.sv_ids for region in self.region_list])
                    print(sum([region.sv_ids for region in self.region_list], []))
                    for region in self.region_list:
                        print(region)
                        for sv in region.sv_list:
                            print(sv.id)
                return 'single_translocation'
        else:
            chromo_metrics = chromothripsis.evaluate_region(regions=self.region_list, call='group', search=False)
            if chromo_metrics['cluster'] and chromo_metrics['complete_walk']:
                return 'chromothripsis'
            return 'group'

    @property
    def enhancer_num(self):
        return len([x for region in self.region_list for x in region.enhancer_list])

    @property
    def super_enhancer_num(self):
        return len([x for region in self.region_list for x in region.super_enhancer_list])

    @property
    def gene_num(self):
        return len([x.tran for region in self.region_list for x in region.bkp_list if x.tran])

    @property
    def sv_num(self):
        return len(set(sum([region.sv_ids for region in self.region_list], [])))

    @property
    def region_num(self):
        return len(self.region_list)

    def has_sv(self, sv):
        bkp_5p = Breakpoint(sv.chrom_5p, sv.bkpos_5p,
                            sv.bkpos_5p_orientation, sv.orientation,
                            sv.chrom_3p, sv.bkpos_3p, '5p')
        bkp_3p = Breakpoint(sv.chrom_3p, sv.bkpos_3p,
                            sv.bkpos_3p_orientation, sv.orientation,
                            sv.chrom_5p, sv.bkpos_5p, '3p')
        bs_5p = [bkp_5p.inside_region(region) for region in self.region_list]
        bs_3p = [bkp_3p.inside_region(region) for region in self.region_list]
        bs = [b1 or b2 for b1, b2 in zip(bs_5p, bs_3p)]
        return bs, bs_5p, bs_3p

    def append_sv(self, sv):
        if sv.chrom_5p == sv.chrom_3p and \
                abs(sv.bkpos_5p - sv.bkpos_3p) <= self.linkage_distance:
            start = min(sv.bkpos_5p, sv.bkpos_3p)
            end = max(sv.bkpos_5p, sv.bkpos_3p)
            region = Region(chrom=sv.chrom_5p, start=start, end=end)
            region.append_sv(sv, self.tran_psl_reader)
            self.append_region(region)
        else:
            region = Region(chrom=sv.chrom_5p, start=sv.bkpos_5p, end=sv.bkpos_5p)
            region.append_sv(sv, self.tran_psl_reader)
            self.append_region(region)

            region = Region(chrom=sv.chrom_3p, start=sv.bkpos_3p, end=sv.bkpos_3p)
            region.append_sv(sv, self.tran_psl_reader)
            self.append_region(region)

    def __repr__(self):
        return str(self.region_list)

    def __str__(self):
        return str(self.region_list)


class Region(myio.Record):

    def __init__(self, **kwargs):
        super(Region, self).__init__(**kwargs)

    def set(self, chrom=None, arm=None, start=None, end=None,
            bkp_list=None, sv_list=None, sv_fn=None,
            cn_list=None, cn=None,
            enhancer_list=None, super_enhancer_list=None,
            linkage_distance=None,
            margin=None,
            sample=None,
            **kwargs):
        self.chrom = chrom
        self.arm = arm or ''
        self.start = start
        self.end = end
        self.linkage_distance = linkage_distance or 1000000
        self.margin = margin or 1000
        self.sample = sample

        self.bkp_list = bkp_list or []
        self.sv_list = sv_list or []

        if not self.sv_list and sv_fn:
            self.read_sv_fn(sv_fn)

        self.cn_list = cn_list or []
        self.cn = cn

        self.enhancer_list = enhancer_list or []
        self.super_enhancer_list = super_enhancer_list or []

    '''
    def read_sv_fn(self, sv_fn):
        for meta, sv in sv_base.SVReader(in_fn=sv_fn, sample=self.sample, file_type='csv'):
            not_5p = sv.chrom_5p != self.chrom or \
                sv.bkpos_5p < self.start or sv.bkpos_5p > self.end
            not_3p = sv.chrom_3p != self.chrom or \
                sv.bkpos_3p < self.start or sv.bkpos_3p > self.end
            if not_5p and not_3p:
                continue
            else:
                self.append_sv(sv)
    '''

    def expand(self):
        start = self.start - self.margin
        if start < 1:
            self.start = 1
        else:
            self.start = start
        end = self.end + self.margin
        if end > constants.hg19_fai_bp[self.chrom]:
            self.end = constants.hg19_fai_bp[self.chrom]
        else:
            self.end = end

    @property
    def length(self):
        if self.start is None or self.end is None:
            return ''
        return '{}M'.format((self.end-self.start)/1000000.0)

    @property
    def sv_ids(self):
        return list(set([sv.id for sv in self.sv_list]))

    @property
    def min_bkp_distance(self):
        xs = sorted(list(set([bkp.bkpos for bkp in self.bkp_list if bkp.chrom == self.chrom] + [self.start, self.end])))
        ds = [x - y for x, y in zip(xs[1:], xs[:-1])]
        min_ds = 0
        if ds:
            min_ds = min(ds)
        return min_ds

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if self.start is not None and self.end is not None:
            range = ':{}-{}'.format(self.start/1000000.0, self.end/1000000.0)
        else:
            range = ''
        return '{}{}{}|{}'.format(self.chrom, self.arm, range, self.length)

    def __lt__(self, other):
        if constants.chrs.index(self.chrom) < constants.chrs.index(other.chrom):
            res = True
        elif self.chrom == other.chrom and self.start < other.start:
            res = True
        else:
            res = False
        return res

    def overlap(self, other):
        '''
            ss     se
            --------
                ---------
                os      oe
            or
                ss     se
                --------
            ---------
            os      oe
        '''
        if self.chrom != other.chrom:
            return False
        elif self.end >= other.start - self.linkage_distance and self.start - self.linkage_distance <= other.end:
            return True
        elif self.start - self.linkage_distance <= other.end and self.end >= other.start - self.linkage_distance:
            return True
        else:
            return False

    def merge(self, other):
        self.start = min(self.start, other.start)
        self.end = max(self.end, other.end)

        self.bkp_list = self.bkp_list + other.bkp_list
        self.sv_list = self.sv_list + other.sv_list

        self.enhancer_list = self.enhancer_list + other.enhancer_list
        self.super_enhancer_list = self.super_enhancer_list + other.enhancer_list

    def _append_bkp(self, bkp, tran_psl_reader=None):
        if tran_psl_reader:
            trans_bkp = tran_psl_reader.fetch_genes_tran_by_point(bkp.chrom, bkp.bkpos)
            bkp.tran = tran_psl_reader.choose_tran_from_group(trans_bkp)
            if bkp.tran:
                func_reg, func_loc = bkp.tran.mark_functional_region(bkp.bkpos)
                bkp.tran.func_loc = func_loc

        bkp_list = self.bkp_list
        bkp_list.append(bkp)
        self.bkp_list = sorted(bkp_list)
        return bkp

    def append_sv(self, sv, tran_psl_reader=None):
        if sv.id in self.sv_ids:
            return

        bkp_5p = Breakpoint(sv.chrom_5p, sv.bkpos_5p,
                            sv.bkpos_5p_orientation, sv.orientation,
                            sv.chrom_3p, sv.bkpos_3p,
                            '5p')
        bkp_3p = Breakpoint(sv.chrom_3p, sv.bkpos_3p,
                            sv.bkpos_3p_orientation, sv.orientation,
                            sv.chrom_5p, sv.bkpos_5p,
                            '3p')

        bkp_5p = self._append_bkp(bkp_5p, tran_psl_reader=tran_psl_reader)
        bkp_3p = self._append_bkp(bkp_3p, tran_psl_reader=tran_psl_reader)

        sv.tran_5p = bkp_5p.tran
        sv.tran_3p = bkp_3p.tran

        self.sv_list.append(sv)

    def has_chrom(self, chrom):
        return chrom == self.chrom

    def has_chrom_arm(self, chrom, arm):
        return chrom == self.chrom and arm == self.arm

    def cut_region(self, direction='5p', step=50):
        if direction == '5p':
            self.bkp_list = self.bkp_list[step:]

        if not self.bkp_list:
            return False

        self.start = self.bkp_list[0].bkpos
        self.end = self.bkp_list[-1].bkpos
        return True


class Breakpoint(myio.Record):

    def __init__(self, chrom, bkpos,
                 orientation, join_type,
                 mate_chrom, mate_bkpos,
                 prime, tran=None):
        self.chrom = chrom
        self.bkpos = bkpos
        self.orientation = orientation
        self.join_type = join_type
        self.mate_chrom = mate_chrom
        self.mate_bkpos = mate_bkpos
        self.prime = prime
        self.tran = tran

    def __lt__(self, other):  # less than
        return self.bkpos < other.bkpos

    def inside_region(self, region):
        return self.nearby_region(region, linkage_distance=0)

    def nearby_region(self, region, linkage_distance=1000000):  # 1M
        res = self.chrom == region.chrom \
            and self.bkpos >= region.start - linkage_distance \
            and self.bkpos <= region.end + linkage_distance
        return res

    def __str__(self):
        return '{}_{}_{}_{}_{}'.format(self.chrom, self.bkpos, self.orientation, self.join_type, self.tran)

    def __repr__(self):
        return self.__str__()
