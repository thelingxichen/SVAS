#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Append tag BX:Z to 10x scCNV bam file.

Usage:
    parse_svaba.py call --bam_fn=IN_FILE
    parse_svaba.py -h | --help

Options:
    -h --help               Show this screen.
    --version               Show version.
    --bam_fn=IN_FILE        Path of 10x scCNV bam file.
"""
import docopt
import pysam
import subprocess


def run_call(bam_fn=None, **args):

    fn = bam_fn.replace('.bam', '.BX.bam')
    samfile = pysam.AlignmentFile(bam_fn, "rb")
    with pysam.AlignmentFile(fn, "wb", template=samfile) as bamfile:

        for read in samfile.fetch():
            if read.has_tag('CB'):
                read.tags += [('BX', read.get_tag('CB'))]
            bamfile.write(read)

    cmd = 'samtools index {}'.format(fn)
    subprocess.run(cmd, shell=True)


def run(call=None, **args):
    if call:
        run_call(**args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    new_args = {}
    for k, v in args.items():
        new_args[k.replace('--', '')] = v
    run(**new_args)
