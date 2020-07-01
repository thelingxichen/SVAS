#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Append tag BX:Z to 10x scCNV bam file.

Usage:
    parse_svaba.py call --bam_fn=IN_FILE [--replace_RG=STR]
    parse_svaba.py -h | --help

Options:
    -h --help               Show this screen.
    --version               Show version.
    --bam_fn=IN_FILE        Path of 10x scCNV bam file.
    --replace_RG=STR     Replace string in RG tag
"""
import docopt
import pysam
import subprocess


def run_call(bam_fn=None, replace_RG=None,**args):

    fn = bam_fn.replace('.bam', '.BX.bam')
    samfile = pysam.AlignmentFile(bam_fn, "rb")

    header = str(samfile.header)
    if replace_RG: 
        str1, str2 = replace_RG.split(',')
        header = header.replace(str1, str2)

    with pysam.AlignmentFile(fn, "wb", text=header) as bamfile:

        for i, read in enumerate(samfile.fetch()):

            if read.has_tag('CB'):
                read.tags += [('BX', read.get_tag('CB'))]
            if replace_RG and read.has_tag('RG'):
                tag = read.get_tag('RG').replace(str1, str2)
                read.set_tag('RG', tag)

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
