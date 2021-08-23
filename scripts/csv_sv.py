#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Complex structure variation

Author: Chen Lingxi (chanlingxi@gmail.com)

Usage:
    xx.py call --sv_fn=IN_FILE --seg_fn=IN_FILE [--sample=STR] [--ref=STR]
    xx.py -h | --help

Options:
    -h --help               Show this screen.
    --version               Show version.
    --sv_fn=IN_FILE         Path of sv txt file.
    --seg_fn=IN_FILE        Path of segment txt file.
    --sample=STR            Sample name, [default: '-']
    --ref=STR               Reference version, [default: hg38]
"""
import os

import docopt
import sv
import localhap
from bioutensil import constants


def read_fai_length(fai_fn):
    contig2length = {}
    for line in open(fai_fn, 'r'):
        splits = line.split()
        contig, length = splits[0], int(splits[1])
        contig2length[contig] = length
    return contig2length


def run_call(sv_fn=None, seg_fn=None, sample=None, **args):
    if sample == '-':
        sample = ''
    out_dir, _ = os.path.split(sv_fn)
    out_dir = os.path.join(out_dir, 'group')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    seg2cn = localhap.read_seg_fn(seg_fn)

    contig2length = constants.hg19_fai_bp
    for group_i, group in sv.read_sv_group(sv_fn, 'group'):
        if group_i is None:
            group_i = 1
        contig_seg_list, contig2seg, bp2seg = localhap.get_segments_from_sv_group(group, contig2length)
        res = localhap.encode_graph(group_i, group, contig_seg_list, bp2seg, contig2length, contig2seg, seg2cn)
        if sample:
            lh_fn = os.path.join(out_dir, '{}.{}.lh'.format(sample, group_i))
        else:
            lh_fn = os.path.join(out_dir, '{}.lh'.format(group_i))
        with open(lh_fn, 'w') as f:
            f.write(res)

        res = localhap.encode_junc_db(group, contig2length)
        if sample:
            db_fn = os.path.join(out_dir, '{}.{}.junc.db'.format(sample, group_i))
        else:
            db_fn = os.path.join(out_dir, '{}.junc.db'.format(group_i))
        with open(db_fn, 'w') as f:
            f.write(res)


def run(call=None, **args):
    if call:
        run_call(**args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    new_args = {}
    for k, v in args.items():
        new_args[k.replace('--', '')] = v
    run(**new_args)
