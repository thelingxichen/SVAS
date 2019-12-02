#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read Support.

Usage:
    read_support.py call --bam_fn=IN_FILE --sv_fn=IN_FILE [--out_dir=OUT_DIR]
    read_support.py -h | --help

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
import os
import sv


def run_call(bam_fn=None, sv_fn=None, expand=500,
             out_dir='./', **args):
    out_fn = os.path.join(out_dir, 'junc.reads.txt')
    with open(out_fn, 'w') as f:
        sv_records = sv.read_vcf(sv_fn)
        for i, record in enumerate(sv_records):
            if i != 4:
                pass

            split_reads = get_split_reads(record, bam_fn)
            if not split_reads:
                continue
            f.write(get_sv_info(record, int(len(split_reads)/2)))
            for read in split_reads:
                f.write(read)


def get_split_reads(sv_record, bam_fn):

    samfile = pysam.AlignmentFile(bam_fn, "rb")
    indexed_reads = pysam.IndexedReads(samfile)
    indexed_reads.build()

    split_reads = []
    tags = set()
    for prime in ['5p', '3p']:
        chrom = getattr(sv_record, 'chrom_{}'.format(prime))
        bkpos = getattr(sv_record, 'bkpos_{}'.format(prime))
        for read in samfile.fetch(chrom, bkpos-1, bkpos):
            if not is_split_read(read, bkpos):
                continue
            reads = list(indexed_reads.find(read.query_name))
            if len(reads) != 2:
                continue

            tag1, read_str1 = get_read_info(reads[0], sv_record, prime,
                                            is_split_read(reads[0], bkpos),
                                            reads[1])

            tag2, read_str2 = get_read_info(reads[1], sv_record, prime,
                                            is_split_read(reads[1], bkpos),
                                            reads[0])
            tag = '{}_{}'.format(tag1, tag2)
            if tag in tags:
                continue
            else:
                tags.add(tag)
            split_reads.append(read_str1)
            split_reads.append(read_str2)

    return split_reads


def is_split_read(read, bkpos, map_q=20, read_q=20, distance=20):
    if read.is_supplementary:
        return False

    if read.mapping_quality < map_q:
        return False

    if read.cigartuples[0][0] == 0 and read.cigartuples[-1][0] == 0:
        # start with M and end with M
        return False
    if read.cigartuples[0][0] in [4, 5]:
        '''
        114H36M or 61S89M
        from left to right
        '''
        if read.cigartuples[0][1] >= distance \
                and read.cigartuples[0][1] <= 150 - distance:
            return True
        else:
            return False
    if read.cigartuples[-1][0] in [4, 5]:
        '''
        reference_end points to one past the last aligned residue.
        114M36S or 61M89S
        from right to left
        '''
        if read.cigartuples[-1][1] >= distance \
                and read.cigartuples[-1][1] <= 150 - distance:
            return True
        else:
            return False

    if np.mean(read.query_qualities) < read_q:
        '''
        read sequence base qualities, including soft clipped bases
        no offset of 33 needs to be subtracted.
        '''
        return False
    return True


def get_read_info(read, record, prime, is_split, paired_read):
    query_qualities = ''.join(chr(x + 33) for x in read.query_qualities)
    flag = ''
    if is_split:
        flag += 'J'
    else:
        flag += 'P'
    if read.is_reverse:
        flag += 'r'
    if paired_read.is_reverse:
        flag += 'R'
    pos1 = read.reference_start+1 - record.bkpos_5p
    tmp = len(record.inner_ins) if record.inner_ins else 0
    pos2 = read.reference_start+1 - record.bkpos_3p - tmp
    if abs(pos1) < abs(pos2):
        pos = pos1
    else:
        pos = pos2

    res = '{}/{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        read.query_name,
        1 if read.is_read1 else 2,
        flag,
        pos,
        read.query_sequence,
        query_qualities,
        read.get_tag('BX')
    )
    tag = '{}_{}'.format(read.cigarstring, pos)
    return tag, res


def get_sv_info(record, junc_reads):
    res = '##{}\t{}\t{}\t{}\tInner-Ins:{}\t{}\n'.format(
        record.chrom_5p,
        record.bkpos_5p,
        record.chrom_3p,
        record.bkpos_3p,
        str(record.inner_ins).upper(),
        junc_reads
    )
    return(res)


def run(call=None, **args):
    if call:
        run_call(**args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    new_args = {}
    for k, v in args.items():
        new_args[k.replace('--', '')] = v
    run(**new_args)
