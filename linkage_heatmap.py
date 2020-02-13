#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Linkage Heatmap.

Usage:
    linkage_heatmap.py call --bam_fn=IN_FILE --sv_fn=IN_FILE [--out_dir=IN_DIR] [--tenx=STR]
    linkage_heatmap.py -h | --help

Options:
    -h --help           Show this screen.
    --version           Show version.
    --bam_fn=IN_FILE    Path of tumor bam file.
    --sv_fn=IN_FILE     Path of tumor sv vcf file.
    --out_dir=OUT_DIR   Path of out directory, [default: ./] 
    --tenx=BOOL         Boolean value, True if the protocal is 10x, otherwise False.
"""
import pysam
import os
import numpy as np

import docopt
import sv


class Region():

    def __init__(self, chrom, start, end, sv_record):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.sv_list = [sv_record]


def run_call(bam_fn=None, sv_fn=None, resolution=50,
             expand=1000, tenx=None,
             out_dir='./', **args):
    if tenx is None or tenx == 'True':
        tenx = True
    else:
        tenx = False

    regions = get_regions(sv_fn, expand)

    samfile = pysam.AlignmentFile(bam_fn)

    for i, region_pair in enumerate(regions):
        linkages = {}
        for r in region_pair:
            linkages_new = collect_linkages(samfile, r, resolution, tenx)
            linkages.update(linkages_new)

        if tenx:
            linkage_types = ['10x', 'pair_end']
        else:
            linkage_types = ['pair_end']

        for linkage_type in linkage_types:
            matrixs = []
            for x_region in region_pair:
                matrix_list = []
                for y_region in region_pair:
                    matrix = get_linkage_matrix(linkages,
                                                x_region, y_region, resolution,
                                                linkage_type=linkage_type)
                    matrix_list.append(matrix)
                matrixs.append(matrix_list)

            out_fn = os.path.join(out_dir, '{}.txt'.format(linkage_type))
            sv_list = sum([region.sv_list for region in region_pair], [])
            if i == 0:
                mode = 'w'
            else:
                mode = 'a'
            write_heatmap(sv_list, matrixs, region_pair, resolution, out_fn, mode,
                          linkage_type)


def get_regions(sv_fn, expand):
    for sv_record in sv.read_vcf(sv_fn):
        region_5p = Region(sv_record.chrom_5p,
                           sv_record.bkpos_5p-expand,
                           sv_record.bkpos_5p+expand,
                           sv_record)
        region_3p = Region(sv_record.chrom_3p,
                           sv_record.bkpos_3p-expand,
                           sv_record.bkpos_3p+expand,
                           sv_record)

        if sv_record.chrom_5p == sv_record.chrom_3p:
            if sv_record.bkpos_5p < sv_record.bkpos_3p:
                yield (region_5p, region_3p)
            else:
                yield (region_3p, region_5p)
        else:
            chroms = list(map(str, range(1, 23))) + ['X', 'Y']
            if chroms.index(sv_record.chrom_5p.replace('chr', '')) < chroms.index(sv_record.chrom_3p.replace('chr', '')):
                yield (region_5p, region_3p)
            else:
                yield (region_3p, region_5p)
            

def get_linkage_matrix(linkages, x_region, y_region, resolution, linkage_type=None):
    x_start = int(x_region.start/resolution)
    x_end = int(x_region.end/resolution)
    x_length = x_end - x_start

    y_start = int(y_region.start/resolution)
    y_end = int(y_region.end/resolution)
    y_length = y_end - y_start

    matrix = np.zeros((x_length, y_length))

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            x_linkages = linkages.get((linkage_type, x_region.chrom, x_start + i), set())
            y_linkages = linkages.get((linkage_type, y_region.chrom, y_start + j), set())
            shared = x_linkages & y_linkages
            if None in shared:
                shared.remove(None)
            matrix[i][j] = len(shared)

    return matrix


def get_read_linkage(read, linkage_type):
    if linkage_type == 'pair_end':
        # read id
        return read.query_name
    if linkage_type == '10x':
        # some read might not have 'BX:Z' tag
        if read.has_tag('BX:Z'):
            return read.get_tag('BX:Z')
        else:
            return None


def collect_linkages(samfile, region, resolution, tenx):
    linkages = {}
    for read in samfile.fetch(region.chrom, region.start, region.end):
        ref_start = read.pos
        ref_end = read.reference_end
        # igonre I and D first ???

        read_linkage = get_read_linkage(read, linkage_type='pair_end')
        if tenx:
            barcode_linkage = get_read_linkage(read, linkage_type='10x')
            linkage_pair = [('10x', barcode_linkage),
                            ('pair_end', read_linkage)]
        else:
            linkage_pair = [('pair_end', read_linkage)]
        for pos in range(ref_start, ref_end + 1):
            for linkage_type, linkage in linkage_pair:
                key = (linkage_type, region.chrom, pos/resolution)
                if key in linkages:
                    # update
                    linkages[key].add(linkage)
                else:
                    linkages[key] = set([linkage])

    return linkages


def write_heatmap(sv_list, matrixs, region_pair, resolution, out_fn, mode, linkage_type):
    with open(out_fn, mode) as f:
        f.write('#sv\n')
        f.write('{}\n'.format(sv_list[0]))
        f.write('#heatmap\tlinkage_type={}\n'.format(linkage_type))
        for j, m_list in enumerate(matrixs):
            y_region = region_pair[j]
            y_region_repr = '{}:{}-{}'.format(y_region.chrom, y_region.start, y_region.end)
            for i, m in enumerate(m_list):
                x_region = region_pair[i]
                x_region_repr = '{}:{}-{}'.format(x_region.chrom, x_region.start, x_region.end)

                f.write('h={}\tv={}\tresolution={}\taxis={},{}\n'.format(
                    y_region_repr, x_region_repr, resolution, i, j))
                for row in m:
                    f.write(','.join(str(x) for x in row) + '\n')


def run(call=None, **args):
    if call:
        run_call(**args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    new_args = {}
    for k, v in args.items():
        new_args[k.replace('--', '')] = v
    run(**new_args)
