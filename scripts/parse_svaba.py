#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parse SvABA sv vcf file

Usage:
    parse_svaba.py call --sv_fn=IN_FILE
    parse_svaba.py -h | --help

Options:
    -h --help           Show this screen.
    --version           Show version.
    --sv_fn=IN_FILE     Path of SvABA sv vcf file.
"""
import docopt
import sv


def run_call(sv_fn=None, **args):
    out_fn = sv_fn + '.txt'
    with open(out_fn, 'w') as f:
        sv_records = sv.read_vcf(sv_fn, precise=False)
        for i, record in enumerate(sv_records):
            f.write(str(record) + '\n')


def run(call=None, **args):
    if call:
        run_call(**args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    new_args = {}
    for k, v in args.items():
        new_args[k.replace('--', '')] = v
    run(**new_args)
