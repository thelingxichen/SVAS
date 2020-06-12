#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate CNV and depth of given genomic region

Usage:
    xx.py call --bam_fn=IN_FILE [--meta_fn=IN_FILE] [--region_fn=IN_FILE] [--region=STR]
    xx.py -h | --help

Options:
    -h --help               Show this screen.
    --version               Show version.
    --bam_fn=IN_FILE        Path of 10x scCNV bam file.
    --meta_fn=IN_FILE       Path of meta file includes single cell subgroup information.
    --bed_fn=IN_FILE        Path of bed file includes targeted genomic regions.
    --region=STR            Bed format of genomic region.
"""
import docopt
import pysam
import pandas as pd


def run_call(bam_fn=None, meta_fn=None,
             bed_fn=None, region=None,
             **args):

    if meta_fn:
        meta_df = pd.read_csv(meta_fn)
        bc2group = dict(zip(meta_df['cell_id'], meta_df['hcluster']))
        groups = list(set(bc2group.values()))
    else:
        bc2group = None
        groups = []

    samfile = pysam.AlignmentFile(bam_fn, "rb")

    regions = []
    if bed_fn:
        regions += open(bed_fn, 'r').read().split('\n')
    if region:
        regions.append(region)
    for region in regions:
        chrom, pos = region.split(':')
        start, end = map(int, pos.split('-'))
        call_depth(samfile, chrom, start, end, bc2group, groups)


def call_depth(samfile, chrom, start, end, bc2group, groups):
    print(chrom, start, end)
    groups.append('bulk')
    d_list = []
    for pileupcolumn in samfile.pileup(chrom, start, end,
                                       min_base_quality=5,
                                       min_mapping_quality=10):
        pos = pileupcolumn.pos
        group2depth = {g: 0 for g in groups}
        for i, pileupread in enumerate(pileupcolumn.pileups):
            '''
            print(not pileupread.is_del, not pileupread.is_refskip)
            if not pileupread.is_del and not pileupread.is_refskip:
                continue
                # query position is None if is_del or is_refskip is set.
            '''
            read = pileupread.alignment
            bc = read.get_tag('CB') if read.has_tag('CB') else None
            group2depth['bulk'] += 1
            if bc in bc2group:
                group2depth[bc2group[bc]] += 1

        group2depth['chrom'] = chrom
        group2depth['pos'] = pos
        d_list.append(group2depth)
    depth_df = pd.DataFrame(d_list)
    for g in groups:
        d = get_average_depth(depth_df[g], start, end)
        print(g, d)


def get_average_depth(depths, start, end):
    xs = depths.to_list()
    xs.sort()
    return sum(xs[1:-1])/(end-start)


def run(call=None, **args):
    if call:
        run_call(**args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    new_args = {}
    for k, v in args.items():
        new_args[k.replace('--', '')] = v
    run(**new_args)
