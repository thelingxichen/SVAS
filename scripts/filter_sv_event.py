#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read Support.

Usage:
    filter_sv_event.py call --sv_fn=IN_FILE [--out_dir=OUT_DIR]
    filter_sv_event.py -h | --help

Options:
    -h --help           Show this screen.
    --version           Show version.
    --bam_fn=IN_FILE    Path of tumor bam file.
    --sv_fn=IN_FILE     Path of tumor sv vcf file.
    --out_dir=OUT_DIR   Path of out directory, [default: ./]
"""
import numpy as np
import docopt
import pysam
import collections
import os
import sv


def run_call(sv_fn=None, expand=500,
             out_dir='./', **args):
    out_fn = os.path.join(out_dir, 'sv.vcf')
    with open(out_fn, 'w') as f:
        sv_records = sv.read_vcf(sv_fn)
        pos_count = 0
        for i, record in enumerate(sv_records):
            if '-' not in [record.strand_5p, record.strand_3p]:
                if pos_count > 4:
                    continue
                if not record.inner_ins:
                    continue
                else:
                    pos_count += 1
            else:
                if not record.inner_ins:
                        continue
            '''
            print('{}:{}-{}'.format(record.chrom_5p, record.bkpos_5p-1000, record.bkpos_5p+1000))
            print('{}:{}-{}'.format(record.chrom_3p, record.bkpos_3p-1000, record.bkpos_3p+1000))
            '''
            print(record)



def run(call=None, **args):
    if call:
        run_call(**args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    new_args = {}
    for k, v in args.items():
        new_args[k.replace('--', '')] = v
    run(**new_args)

