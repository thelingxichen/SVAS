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
import collections
import os
import sv


def run_call(bam_fn=None, sv_fn=None, expand=500,
             out_dir='./', **args):
    out_fn = os.path.join(out_dir, 'junc.reads.txt')
    with open(out_fn, 'w') as f:
        sv_records = sv.read_vcf(sv_fn)
        for i, record in enumerate(sv_records):
            if record.inner_ins and len(record.inner_ins) >= 10:
                continue
            print(i, '--- sv_record\n', record)

            split_reads = get_split_reads(record, bam_fn)
            if not split_reads or len(split_reads) < 5:
                continue
            f.write(get_sv_info(record, int(len(split_reads)/2)))
            for j, read in enumerate(split_reads):
                read, read_str = read
                f.write(read_str)


def get_reads(sv_record, samfile):

    for prime in ['5p', '3p']:
        chrom = getattr(sv_record, 'chrom_{}'.format(prime))
        bkpos = getattr(sv_record, 'bkpos_{}'.format(prime))
        for read in samfile.fetch(chrom, bkpos-1, bkpos):
            yield prime, bkpos, read


def get_split_reads(sv_record, bam_fn):

    samfile = pysam.AlignmentFile(bam_fn, "rb")
    indexed_reads = pysam.IndexedReads(samfile)
    indexed_reads.build()

    split_reads = []
    tags = set()
    homseq = sv_record.meta_info.get('HOMSEQ')
    pos = []
    for prime, bkpos, read in get_reads(sv_record, samfile):
        if not is_split_read_for_bkpos(prime, bkpos, read):
            continue
        print(read)
        reads = list(indexed_reads.find(read.query_name))
        if len(reads) == 1:
            continue

        if len(reads) == 3:
            res = deal_with_three_reads(read, reads, sv_record, prime)
        if len(reads) == 2:
            res = deal_with_two_reads(read, reads, sv_record, prime)

        if not res:
            continue
        tag, split_reads_tmp = res

        if tag in tags:
            continue
        else:
            tags.add(tag)

        split_reads = split_reads + split_reads_tmp

        for r, _ in split_reads_tmp:
            pos = pos + get_bkpos(r)

    most_common_bkpos = collections.Counter(pos).most_common(1)
    if not most_common_bkpos:
        return []
    most_common_bkpos = collections.Counter(pos).most_common(1)[0][0][1]
    if homseq and sv_record.bkpos_5p == most_common_bkpos:
        sv_record.bkpos_3p = sv_record.bkpos_3p + len(homseq)
    if homseq and sv_record.bkpos_3p == most_common_bkpos:
        sv_record.bkpos_5p = sv_record.bkpos_5p - len(homseq)
    return split_reads


def deal_with_three_reads(read, reads, sv_record, prime):
    bkpos = getattr(sv_record, 'bkpos_{}'.format(prime))
    if 'H' in read.cigarstring:
        return None

    pair_read = [r for r in reads if r.is_read1 != read.is_read1][0]
    split_read = [r for r in reads if r.is_read1 == read.is_read1
                  and (r.cigarstring != read.cigarstring)]
    if not split_read:
        reads = [read, pair_read]
        return deal_with_two_reads(read, reads, sv_record, prime)

    split_read = split_read[0]
    print(pair_read)
    print(split_read)

    split_prime, split_pos = get_bkpos(split_read)[0]
    split_bkpos = getattr(sv_record, 'bkpos_{}'.format(split_prime))
    print(split_bkpos - split_pos)
    if abs(split_bkpos - split_pos) > 100:
        return None

    tag1, read_str1 = get_read_info(read, sv_record, prime,
                                    is_split_read_for_bkpos(prime, bkpos, read),
                                    pair_read)

    tag2, read_str2 = get_read_info(pair_read, sv_record, prime,
                                    is_split_read_for_bkpos(prime, bkpos, pair_read),
                                    read)
    if tag1 is None or tag2 is None:
        return None
    tag = '{}_{}'.format(tag1, tag2)
    return tag, [(read, read_str1), (pair_read, read_str2)]


def deal_with_two_reads(reads, chrom, bkpos, prime):
    pass


def get_bkpos(read):
    pos = []
    if read.cigartuples[0][0] in [4, 5]:
        '''
        114H36M or 61S89M
        from left to right
        '''
        pos.append(('3p', read.pos + read.cigartuples[0][1]))
    if read.cigartuples[-1][0] in [4, 5]:
        '''
        reference_end points to one past the last aligned residue.
        114M36 or 61M89S
        from right to left
        '''
        pos.append(('5p', read.pos + len(read.query_sequence) - read.cigartuples[-1][1]))
    return pos


def is_split_read_for_bkpos(prime, bkpos, read, map_q=0, read_q=0, distance=10):
    if read.is_supplementary:
        return False

    if read.mapping_quality < map_q:
        return False

    if read.cigartuples[0][0] == 0 and read.cigartuples[-1][0] == 0:
        # start with M and end with M
        return False

    if np.mean(read.query_qualities) < read_q:
        '''
        read sequence base qualities, including soft clipped bases
        no offset of 33 needs to be subtracted.
        '''
        return False

    pos = get_bkpos(read)
    flag = False
    for p in pos:
        prime_tmp, p = p
        if prime_tmp == prime and abs(p - bkpos) < distance:
            flag = True

    return flag


def get_read_info(read, record, prime, is_split, paired_read):
    query_qualities = ''.join(chr(x + 33) for x in read.query_qualities)
    flag = ''
    if is_split:
        flag += 'j'
    else:
        flag += 'J'  #
    if read.is_reverse:
        flag += 'r'
    if paired_read.is_reverse:
        flag += 'R'
    pos1 = read.reference_start + 1 - record.bkpos_5p
    tmp = len(record.inner_ins) if record.inner_ins else 0
    pos2 = read.reference_start + 1 - record.bkpos_3p - tmp - 2
    if abs(pos1) < abs(pos2):
        pos = pos1
    else:
        pos = pos2
    if abs(pos) > 500:
        return None, None

    if 'BX' in read.get_tags():
        barcode = read.get_tag('BX')
    else:
        barcode = ''
    res = '{}/{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        read.query_name,
        1 if read.is_read1 else 2,
        flag,
        pos,
        read.query_sequence,
        query_qualities,
        barcode
    )
    tag = '{}_{}'.format(read.cigarstring, pos)
    return tag, res


def get_sv_info2(record, junc_reads):

    meta_info = 'Inner-ins:{}'.format(
        str(record.inner_ins).upper()
    )
    res = '##{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        record.chrom_5p,
        record.bkpos_5p,
        record.chrom_3p,
        record.bkpos_3p,
        meta_info,
        junc_reads
    )
    return(res)


def get_sv_info(record, junc_reads):
    barcodes = record.meta_info.get('BX')
    if barcodes:
        barcode = ','.join(barcodes)
    else:
        barcode = 'NONE'
    meta_info = 'Inner-ins:{},HM:{},BX:{}'.format(
        str(record.inner_ins).upper(),
        'NONE',
        'NONE'
    )
    res = '##{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        record.chrom_5p,
        record.bkpos_5p,
        record.strand_5p,
        record.chrom_3p,
        record.bkpos_3p,
        record.strand_3p,
        meta_info,
        junc_reads,
        barcode
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
