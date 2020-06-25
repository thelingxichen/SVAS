#!/usr/bin/env python
# -*- coding: utf-8 -*-
import importlib
import numpy as np
import collections

from biotool.mymath import significant_test

from biotool import constants

from compsv import base

'''
class CNVState():

    def __init__(self, reader_cls=cnv.CNVReader, writer_cls=myio.Writer, *args, **kwargs):
        pass

    def __iter__(self):
        res = {}
        # look through
        for i, record in enumerate(self.in_iter):
            meta, record = record
            region = record.chr + record.arm
            state_set = res.get(region, set([]))
            state_set.add(record.Cn)
            res[region] = state_set

        # merge p-arm and q-arm if states are same
        for chrom in constants.chrs:
            p = '{}p'.format(chrom)
            q = '{}q'.format(chrom)
            if p in res and q in res:
                if res[p] == res[q]:
                    region = base.Region(chrom=chrom)
                    res[str(region)] = res[p]
                else:
                    res[str(base.Region(chrom=chrom, arm='p'))] = res[p]
                    res[str(base.Region(chrom=chrom, arm='q'))] = res[q]
                del res[p]
                del res[q]

            elif p in res:
                res[str(base.Region(chrom=chrom, arm='p'))] = res[p]
                del res[p]
            elif q in res:
                res[str(base.Region(chrom=chrom, arm='q'))] = res[q]
                del res[q]

        for region, state_set in res.items():
            if len(state_set) == 2:
                args = {
                    'region': region,
                    'cn_state': list(state_set)
                }
                yield None, myio.Record(args=args)
'''


class SVBkpLocation():

    def __init__(self, *args, **kwargs):
        pass

    def __iter__(self):
        width = len(constants.chrs)+1
        matrix = np.zeros([width, width])
        for i, record in enumerate(self.in_iter):
            meta, record = record
            if record.chrom_5p in constants.chrs:
                i = constants.chrs.index(record.chrom_5p)
            else:
                i = width-1
            if record.chrom_3p in constants.chrs:
                j = constants.chrs.index(record.chrom_3p)
            else:
                j = width-1

            matrix[i][j] += 1
            matrix[j][i] += 1

        yield None, matrix


class SVBkpRegion():

    def __init__(self, regions=None, sv_records=None, if_search=None, *args, **kwargs):
        self.if_search = if_search
        self.regions = regions or []
        self.sv_records = sv_records or []

    def call_single(self):
        self._collect_svs()

        for region in self.regions:

            if self.if_search:
                metrics = self.search(region)
            else:
                metrics = self.chromothripsis(region)

            yield metrics

    def call_group(self):
        self._collect_svs()

        return self.chromothripsis()

    def _collect_svs(self):
        for i, record in enumerate(self.sv_records):
            meta, record = record
            if record.chrom_5p == record.chrom_3p and \
                    abs(record.bkpos_5p - record.bkpos_3p) < 1000:
                continue
            if record.chrom_5p not in constants.chrs or record.chrom_3p not in constants.chrs:
                continue
            if record.chrom_5p != 'chr7':
                pass

            b, bs = self._sv_in_regions(record)
            if b:
                self._append_sv_to_regions(record, bs)

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

    def chromothripsis(self, region=None):
        # check if breakpoint cluster together
        cluster, cluster_pvalue = self.bkp_cluster(region)

        # check if breakpoint join randomly
        random_join, random_join_pvalue = self.bkp_random_join(region)

        # check the ability to complete walk
        complete_walk, complete_walk_pvalue, complete_walk_runs = self.bkp_complete_walk(region)

        chromothripsis = cluster and random_join and complete_walk
        # genes = set([bp.tran.gene for bp in region.bkp_list if bp.tran])
        # print region, len(region.bkp_list), len(genes), chromothripsis, cluster, random_join, complete_walk, complete_walk_runs, region.bkp_join, ','.join(genes)
        # print (region, len(region.bkp_list), len(genes), ','.join(genes))

        args = {
            'region': region,
            'chromothripsis': chromothripsis,
            'cluster': cluster,
            'cluster_pvalue': cluster_pvalue,
            'random_join': random_join,
            'random_join_pvalue': random_join_pvalue,
            'complete_walk': complete_walk,
            'complete_walk_pvalue': complete_walk_pvalue,
            'complete_walk_runs': complete_walk_runs,
            # 'bkp_list': bkp_list,
            # 'bkp_join': bkp_join,
        }
        return args  # chromothripsis, myio.Record(args=args)

    def search(self, region):
        print('in search')
        chromothripsis, record = self.chromothripsis(region)

        if chromothripsis:  # base case
            return chromothripsis, record
        else:
            cut = region.cut_region()
            if not cut:  # base case
                return chromothripsis, None

            # recursive step
            return self.search(region)

    def bkp_cluster(self, region=None):
        '''
        The null model of random breakpoint locations implies
        that the distances between adjacent breakpoints,
        should be distributed according to an exponential distribution,
        which can be readily evaluated using a goodness-of-fit test.

        In our experience, chromothripsis is typically associated with
        a strong departure from this null distribution,
        although some situations of progressive rearrangements

        Reject: might exist chromothripsis
            not from exponetial distribution,
            breakpoint cluster together
        '''
        if region:
            bkp_list = [bkp.bkpos for bkp in region.bkp_list]
        else:
            bkp_list = [bkp.bkpos for r in self.regions for bkp in r.bkp_list]

        adjacent_distances = [x - y for x, y in zip(bkp_list[1:], bkp_list[:-1])]
        reject, pvalue = significant_test.exponent_distribution(adjacent_distances)
        if reject is None:
            reject = False
        return reject, pvalue

    def bkp_random_join(self, region=None):
        """
        Collect the counts of observed rearrangements that have a
        deletion-type, tandem duplication-type, head-to-head-inverted,
        and tail-to-tail-inverted orientation respectively.

        If more than one chromosome is involved, then interchromosomal rearrangements
        can be interpreted in the same four categories using
        orientation of the strands at the breakpoint.

        Then, in a region of chromothripsis,
        we would expect these counts to be distributed as a multinomial distribution
        with parameters n=\sum{r_i} and probability p_i=1/4.

        A departure from this distribution can be employed as evidence against
        the rearrangements arising from a chromothripsis process.

        Null hypothesis:
            the count of following join types fit multinomial distribution.
            intra-chromosome:   DEL-ht, DUP-th, INV-tt, INV-hh
            inter-chromosome:   TRX-ht, TRX-th, TRX-tt, TRX-hh
            =>
            the orientation "ht, th, hh, tt" fit multinomial distribution
        Accept: might exist chromothripsis.
            breakpoints joint randomly
        """
        if region:
            bkp_join = [bkp.join_type for bkp in region.bkp_list]
        else:
            bkp_join = [bkp.join_type for r in self.regions for bkp in r.bkp_list]

        bkp_join = collections.Counter(bkp_join)
        for orientation in 'ht,th,tt,hh'.split(','):
            if orientation not in bkp_join:
                bkp_join[orientation] = 0

        reject, pvalue = significant_test.multinomial_distribution(list(bkp_join.values()))

        if region:
            region.bkp_join = bkp_join
        return not reject, pvalue

    def bkp_complete_walk(self, region=None):
        '''
        chromothripsis pattern on reference
            |++A++||++B++||++C++||++D++||++E++|
                  h       t     h
                   t     h             h
                                 t      t
        -------------------------------------------
                  ht     ht     ht     ht
        => h, t alternating
        tumor:
            |++A++||++C++||++B++||--D--||++E++|
                  35     35     33     55
                  ht     ht     hh     tt
                  DEL    DUP        head-to-head INV
                          tail-to-tail INV

        accumulating pattern on reference:
            |++A++||++B++||++C++||++D++||++E++|
                  h       t     h
                   t            h      h
                          t             t
        --------------------------------------------
                         h       t
        impossible to append, cannot go complete walk   # formulate to a math puzzle and prove it???
        --------------------------------------------
                  ht     2t    2h      ht               # 2t, 2h should not be common in real case

        tumor:
            |++A++||++C++||++B++C++||--D--C--||++E++|
                  35     35        33        55
                  ht     ht        hh        tt
                  DEL    DUP           head-to-head INV
                             tail-to-tail INV

        '''
        if region:
            bkp_orientation = [bkp.orientation for bkp in region.bkp_list]
        else:
            bkp_orientation = [bkp.orientation for r in self.regions for bkp in r.bkp_list]

        reject, pvalue, n_runs = significant_test.wald_wolfowitz_test(bkp_orientation)
        if reject is None:
            reject = True
        return not reject, pvalue, n_runs


'''
def run_cnv_state(**args):

    cnv_fn = args['cnv_fn']
    cnv_states = CNVState(
        in_fn=cnv_fn,
        write=True,
        suffix='cn_state',
        **args
    ).evaluate()

    return cnv_states, None
'''


def run_bkp_location(sv_module, **args):
    reader_cls = sv_module.SVReader
    sv_fn = args['sv_fn']
    sv_records = sv_module.SVGenerator(
        reader_cls=reader_cls,
        in_fn=sv_fn,
        **args
    )
    SVBkpLocation(
        reader_cls=reader_cls,
        in_data=sv_records,
        write=True,
        suffix='bkp_location',
        **args
    ).evaluate()


def evaluate_region(regions=None, search=None, call='single', **args):
    sv_records = []
    for region in regions:
        sv_records += region.sv_list

    sv_bkp_region = SVBkpRegion(
        search=search,
        regions=regions,
        in_data=sv_records,
        write=False,  # True,
        suffix='bkp_region',
        **args
    )

    if call == 'single':
        return sv_bkp_region.call_single()
    else:
        return sv_bkp_region.call_group()


def parse_region(region_module, **args):

    regions = region_module.RegionGenerator(
        in_fn=args['region_fn'],
        **args
    )
    return regions, regions.out_fn


def run_call(sv_module=None, **args):
    # run_bkp_location(sv_module=sv_module, **args)
    '''
    cnv_states, _ = run_cnv_state(**args)

    regions = [x.region for m, x in cnv_states]
    '''
    regions = [base.Region(chrom=chrom, start=0, end=constants.hg19_fai_bp[chrom]) for chrom in constants.chrs]
    evaluate_region(sv_module=sv_module, regions=regions, search=False, **args)


def run_evaluate(sv_module=None, format=None, **args):
    region_module = importlib.import_module('tenxtools.chromothripsis.{}'.format(format))
    regions, out_fn = parse_region(region_module, **args)

    regions = [x[1] for x in regions]
    regions, out_fn = evaluate_region(sv_module=sv_module, regions=regions, search=False, **args)
