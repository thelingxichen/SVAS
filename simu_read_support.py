import numpy as np
import pybedtools
import random


def simu_pe(local_hap, sv_type, span):
    count = random.randint(10000, 99999)
    read_id = 'E00514:190:H3Y2MCCXY:4:1211:13321:{}'.format(
        count)

    read1_junc = np.random.randint(2)
    if read1_junc and sv_type in ['ht', 'th']:
        read1_flag, read2_flag = 'JR', 'jr'
        read1_pos = np.random.randint(-70, -30)
        read2_pos = np.random.randint(160, 310)
    elif not read1_junc and sv_type in ['ht', 'th']:
        read1_flag, read2_flag = 'jR', 'Jr'
        read1_pos = np.random.randint(-310, -160)
        read2_pos = np.random.randint(-70, -30)
    elif read1_junc and sv_type in ['hh', 'tt']:
        read1_flag, read2_flag = 'j', 'J'
        read1_pos = np.random.randint(-70, -30)
        read2_pos = np.random.randint(160, 310)
    elif not read1_junc and sv_type in ['hh', 'tt']:
        read1_flag, read2_flag = 'j', 'J'
        read1_pos = np.random.randint(-310, -160)
        read2_pos = np.random.randint(-70, -30)

    if np.random.randint(2):
        read1_flag += 'd'
        read2_flag += 'd'

    read1_seq = local_hap[span+read1_pos-1: span+read1_pos+150-1]
    read2_seq = local_hap[span+read1_pos-1: span+read2_pos+150-1]

    barcode = ''

    res = '{}/1\t{}\t{}\t{}\t{}\n'.format(
        read_id,
        read1_flag,
        read1_pos,
        read1_seq,
        'A'*150,
        barcode
    )
    res += '{}/2\t{}\t{}\t{}\t{}\n'.format(
        read_id,
        read2_flag,
        read2_pos,
        read2_seq,
        'A'*150,
        barcode
    )
    return res


def simu_local_hap(span, ref_fasta, chrom_5p, bkpos_5p, chrom_3p, bkpos_3p, sv_type, inner_ins):
    if sv_type in ['ht', 'th']:
        start_5p = bkpos_5p - span
        end_5p = bkpos_5p
        start_3p = bkpos_3p - 1  # bed start from 0
        end_3p = bkpos_3p + span
    elif sv_type == 'hh':
        start_5p = bkpos_5p - span
        end_5p = bkpos_5p
        start_3p = bkpos_3p - span
        end_3p = bkpos_3p
    elif sv_type == 'tt':
        start_5p = bkpos_5p - 1
        end_5p = bkpos_5p + span
        start_3p = bkpos_3p - 1
        end_3p = bkpos_3p + span

    bed_str = '{}\t{}\t{}\n{}\t{}\t{}\n'.format(
        chrom_5p, start_5p, end_5p,
        chrom_3p, start_3p, end_3p
    )
    bedtool = pybedtools.BedTool(bed_str, from_string=True)
    ref_fasta = pybedtools.example_filename(ref_fasta)
    bedtool.sequence(fi=ref_fasta)
    header_5p, seq_5p, header_3p, seq_3p = open(bedtool.seqfn, 'r').read().split()
    seq_5p, seq_3p = seq_5p.upper(), seq_3p.upper()
    inner_ins = '' if inner_ins == 'NONE' else inner_ins

    if sv_type in ['ht', 'th']:
        local_hap = seq_5p + inner_ins + seq_3p
    elif sv_type == 'hh':
        local_hap = seq_5p + inner_ins + seq_3p[::-1]
    elif sv_type == 'tt':
        local_hap = seq_5p[::-1] + inner_ins + seq_3p

    return local_hap


def simu_reads(ref_fasta,
               chrom_5p, bkpos_5p, strand_5p,
               chrom_3p, bkpos_3p, strand_3p, sv_type,
               inner_ins, junc_reads):
    span = 400
    local_hap = simu_local_hap(span, ref_fasta, chrom_5p, bkpos_5p, chrom_3p, bkpos_3p, sv_type, inner_ins)

    '''
    meta_info = 'Inner-ins:{},HM:{},BX:{}'.format(
        inner_ins if inner_ins else 'NONE',
        'NONE',
        'NONE'
    )
    res = '##{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        chrom_5p, bkpos_5p, strand_5p,
        chrom_3p, bkpos_3p, strand_3p,
        meta_info, junc_reads
    )
    '''
    meta_info = 'Inner-ins:{}'.format(
        inner_ins if inner_ins else 'NONE'
    )
    res = '##{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        chrom_5p, bkpos_5p, chrom_3p, bkpos_3p,
        meta_info, junc_reads
    )

    for i in range(junc_reads):
        res += simu_pe(local_hap, sv_type, span)
    return res


def main():
    ref_fasta = '/home/BIOINFO_DATABASE/reference/genome_DNA/Homo_sapiens/hg19/BWA_GATK_index/hg19.fa'

    sv_type = 'ht'
    chrom_5p = 'chr1'
    bkpos_5p = 100000
    strand_5p = '+'
    chrom_3p = 'chr1'
    bkpos_3p = 101000
    strand_3p = '+'
    inner_ins = 'TACCGATAT'
    junc_reads = random.randint(10, 20)
    res = simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                     chrom_3p, bkpos_3p, strand_3p, sv_type,
                     inner_ins, junc_reads)
    sv_type = 'th'
    chrom_5p = 'chr20'
    bkpos_5p = 200000
    strand_5p = '+'
    chrom_3p = 'chr20'
    bkpos_3p = 199000
    strand_3p = '+'
    inner_ins = ''
    junc_reads = random.randint(10, 20)
    res += simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                      chrom_3p, bkpos_3p, strand_3p, sv_type,
                      inner_ins, junc_reads)
    sv_type = 'hh'
    chrom_5p = 'chr3'
    bkpos_5p = 200000
    strand_5p = '+'
    chrom_3p = 'chr3'
    bkpos_3p = 200005
    strand_3p = '-'
    inner_ins = ''
    junc_reads = random.randint(10, 20)
    res += simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                      chrom_3p, bkpos_3p, strand_3p, sv_type,
                      inner_ins, junc_reads)
    sv_type = 'tt'
    chrom_5p = 'chr16'
    bkpos_5p = 10000000
    strand_5p = '-'
    chrom_3p = 'chr16'
    bkpos_3p = 10000005
    strand_3p = '+'
    inner_ins = ''
    junc_reads = random.randint(10, 20)
    res += simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                      chrom_3p, bkpos_3p, strand_3p, sv_type,
                      inner_ins, junc_reads)

    sv_type = 'ht'
    chrom_5p = 'chr2'
    bkpos_5p = 100000
    strand_5p = '+'
    chrom_3p = 'chr3'
    bkpos_3p = 101000
    strand_3p = '+'
    inner_ins = 'TACCGATAT'
    junc_reads = random.randint(10, 20)
    res += simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                      chrom_3p, bkpos_3p, strand_3p, sv_type,
                      inner_ins, junc_reads)
    sv_type = 'th'
    chrom_5p = 'chr2'
    bkpos_5p = 2000000
    strand_5p = '+'
    chrom_3p = 'chr3'
    bkpos_3p = 1990000
    strand_3p = '+'
    inner_ins = ''
    junc_reads = random.randint(10, 20)
    res += simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                      chrom_3p, bkpos_3p, strand_3p, sv_type,
                      inner_ins, junc_reads)
    sv_type = 'hh'
    chrom_5p = 'chr1'
    bkpos_5p = 2000000
    strand_5p = '+'
    chrom_3p = 'chr12'
    bkpos_3p = 200005
    inner_ins = ''
    strand_3p = '-'
    junc_reads = random.randint(10, 20)
    res += simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                      chrom_3p, bkpos_3p, strand_3p, sv_type,
                      inner_ins, junc_reads)
    sv_type = 'tt'
    chrom_5p = 'chr1'
    bkpos_5p = 10000000
    strand_5p = '-'
    chrom_3p = 'chr18'
    bkpos_3p = 1000005
    strand_3p = '+'
    inner_ins = ''
    junc_reads = random.randint(10, 20)
    res += simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                      chrom_3p, bkpos_3p, strand_3p, sv_type,
                      inner_ins, junc_reads)
    print(res)


if __name__ == "__main__":
    main()
