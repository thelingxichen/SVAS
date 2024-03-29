#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Mobile element mediated structure variation

Author: Chen Lingxi (chanlingxi@gmail.com)

Usage:
    xx.py call --sv_fn=IN_FILE --ME_fai_fn=IN_FILE
    xx.py -h | --help

Options:
    -h --help               Show this screen.
    --version               Show version.
    --sv_fn=IN_FILE         Path of sv txt file.
    --ME_fai_fn=IN_FILE     Path of fasta index of mobile element.
"""
import os

import docopt
import sv
import localhap
from biotool import constants


def read_fai_length(fai_fn):
    contig2length = {}
    for line in open(fai_fn, 'r'):
        splits = line.split()
        contig, length = splits[0], int(splits[1])
        contig2length[contig] = length
    return contig2length


def run_call(sv_fn=None, ME_fai_fn=None, **args):
    out_dir, _ = os.path.split(sv_fn)
    out_dir = os.path.join(out_dir, 'localhap')
    print(out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    contig2length = read_fai_length(ME_fai_fn)
    contig2length.update(constants.hg19_fai_bp)
    for group_i, group in sv.read_sv_group(sv_fn, 'group'):
        contig_seg_list, contig2seg, bp2seg = localhap.get_segments(group, contig2length)
        res = localhap.encode_graph(group_i, group, contig_seg_list, bp2seg, contig2length, contig2seg)
        lh_fn = os.path.join(out_dir, 'group{}.lh'.format(group_i))
        with open(lh_fn, 'w') as f:
            f.write(res)

        res = localhap.encode_junc_db(group)
        db_fn = os.path.join(out_dir, 'group{}.junc.db'.format(group_i))
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
