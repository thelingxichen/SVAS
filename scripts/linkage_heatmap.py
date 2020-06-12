#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Linkage Heatmap.

Usage:
    linkage_heatmap.py callDNA --bam_fn=IN_FILE --sv_fn=IN_FILE [--out_prefix=STR] --linkage_type=STR
    linkage_heatmap.py callRNA --bam_fn=IN_FILE --sv_fn=IN_FILE [--out_prefix=STR] --linkage_type=STR
    linkage_heatmap.py -h | --help

Options:
    -h --help               Show this screen.
    --version               Show version.
    --bam_fn=IN_FILE        Path of tumor bam file.
    --sv_fn=IN_FILE         Path of tumor sv vcf file.
    --out_prefix=STR        Path of out directory, [default: ./] 
    --linkage_type=STR      Heatmap linkage type, choose from PE, 10x, LR, [default: PE]
"""
import pysam
import os
import numpy as np
import csv

import docopt
import sv


class Region():

    def __init__(self, chrom, start, end, sv_record):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.sv_list = [sv_record]


def run_call(bam_fn=None, sv_fn=None, regions=None, resolution=80,
             linkage_type=None, expand=1000, isDNA=True,
             out_prefix='./', **args):

    if isDNA:
        regions = get_regions_from_vcf(sv_fn, expand)
        seq = 'DNA-Seq'
    else:
        regions = get_regions_from_fusion_tsv(sv_fn, expand)
        seq = 'RNA-Seq'

    samfile = pysam.AlignmentFile(bam_fn)

    for i, region_pair in enumerate(regions):
        linkages = {}
        for r in region_pair:
            linkages_new = collect_linkages(samfile, r, resolution, linkage_type)
            linkages.update(linkages_new)

        for linkage_type in [linkage_type]:
            matrixs = []
            for x_region in region_pair:
                matrix_list = []
                for y_region in region_pair:
                    matrix = get_linkage_matrix(linkages,
                                                x_region, y_region, resolution,
                                                linkage_type=linkage_type)
                    matrix_list.append(matrix)
                matrixs.append(matrix_list)

            out_fn = os.path.join(out_prefix, '{}_res{}_{}.txt'.format(linkage_type, resolution, seq))
            sv_list = sum([region.sv_list for region in region_pair], [])
            if i == 0:
                mode = 'w'
            else:
                mode = 'a'
            write_heatmap(sv_list, matrixs, region_pair, resolution, out_fn, mode,
                          linkage_type, isDNA=isDNA)


def get_regions_from_vcf(sv_fn, expand):
    for sv_record in sv.read_vcf(sv_fn):
        yield get_region(sv_record, expand)
   

def get_regions_from_fusion_tsv(sv_fn, expand):
    for sv_record in read_fusion_tsv(sv_fn):
        yield get_region(sv_record, expand)


def read_fusion_tsv(sv_fn):
    for row in csv.DictReader(open(sv_fn), delimiter='\t'):
        sv_record = sv.SVRecord()
        sv_record.set(chrom_5p=row['up_chr'], bkpos_5p=row['up_pos'], strand_5p=row['up_strand'],
                      chrom_3p=row['down_chr'], bkpos_3p=row['down_pos'], strand_3p=row['down_strand'],
                      meta_info={'up_gene': row['#up_gene'],
                                 'down_gene': row['down_gene'],
                                 'VARTYPE': row['comments']})
        yield sv_record 

def get_region(sv_record, expand):
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
            return (region_5p, region_3p)
        else:
            return (region_3p, region_5p)
    else:
        chroms = list(map(str, range(1, 23))) + ['X', 'Y']
        if chroms.index(sv_record.chrom_5p.replace('chr', '')) < chroms.index(sv_record.chrom_3p.replace('chr', '')):
            return (region_5p, region_3p)
        else:
            return (region_3p, region_5p)
            

def get_linkage_matrix(linkages, x_region, y_region, resolution, linkage_type=None):
    x_length = int((x_region.end - x_region.start)/resolution)
    y_length = int((y_region.end - y_region.start)/resolution)

    matrix = np.zeros((x_length, y_length))

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            win_x_start = x_region.start + i*resolution
            win_x_end = win_x_start + resolution

            x_linkages = get_linkages_from_region(linkages, x_region.chrom, win_x_start, win_x_end, linkage_type, resolution)


            win_y_start = y_region.start + j*resolution
            win_y_end = win_y_start + resolution 
            y_linkages = get_linkages_from_region(linkages, y_region.chrom, win_y_start, win_y_end, linkage_type, resolution)


            '''
            x_linkages = x_linkages + linkages.get((linkage_type, x_region.chrom, x_start + i), set())
            y_linkages = linkages.get((linkage_type, y_region.chrom, y_start + j), set())
            '''
            shared = x_linkages & y_linkages
            if None in shared:
                shared.remove(None)
            matrix[i][j] = len(shared)

    return matrix

def get_linkages_from_region(linkages, chrom, start, end, linkage_type, resolution):
    res = set()
    for i in range(start, end):
        key = (linkage_type, chrom, i/resolution)
        res = res | linkages.get(key, set())
    return res
        

def get_read_linkage(read, linkage_type):
    if linkage_type != '10x':
        # read id
        return read.query_name
    if linkage_type == '10x':
        # some read might not have 'BX:Z' tag
        if read.has_tag('BX:Z'):
            return read.get_tag('BX:Z')
        else:
            return None


def collect_linkages(samfile, region, resolution, linkage_type):
    linkages = {}
    for read in samfile.fetch(region.chrom, region.start, region.end):
        if not read.reference_end:
            continue
        ref_start = read.pos
        ref_end = read.reference_end
        # igonre N first 
        pos_list = []
        current_pos = ref_start
        for i, p in enumerate(read.cigartuples):
            t, bp = p
            if t == 3:  # N
                current_pos = current_pos + bp
            elif i == 0 and t == 5:  # start with H
                pos_list += range(current_pos-bp, current_pos)
            else:  # M, D, I, S, endswith H
                pos_list += range(current_pos, current_pos + bp)
                current_pos = current_pos + bp

        linkage = get_read_linkage(read, linkage_type=linkage_type)
        for pos in pos_list:
            key = (linkage_type, region.chrom, pos/resolution)
            if key in linkages:
                # update
                linkages[key].add(linkage)
            else:
                linkages[key] = set([linkage])

    return linkages


def write_heatmap(sv_list, matrixs, region_pair, resolution, out_fn, mode, linkage_type, isDNA=True):
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

                if i == 0 and j == 0:
                    axis = 'left-top'
                if i == 1 and j == 0:
                    axis = 'right-top'
                if i == 0 and j == 1:
                    axis = 'left-bottom'
                if i == 1 and j == 1:
                    axis = 'right-bottom'
                f.write('v={}\th={}\tresolution={}\taxis={}\n'.format(
                    y_region_repr, x_region_repr, resolution, axis))
                for row in m:
                    f.write(','.join(str(x) for x in row) + '\n')


def run(callDNA=None, callRNA=None, **args):
    if callDNA:
        run_call(isDNA=True, **args)
    if callRNA:
        run_call(isDNA=False, expand=1000, resolution=50, **args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    new_args = {}
    for k, v in args.items():
        new_args[k.replace('--', '')] = v
    run(**new_args)
