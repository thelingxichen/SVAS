import numpy as np
import pybedtools
import random


def simu_single_read(local_hap, read_id, read_flag, read_pos, span, barcode, is_read1):
    read_num = 127 if is_read1 else 150

    indels = []
    if np.random.randint(2):
        num_indels = np.random.randint(1, 5)
        for i in range(num_indels):
            indel_type = 'INS' if np.random.randint(2) else 'DEL'
            indel_len = np.random.randint(1, 4)
            indel_candidates = list(range(read_pos, -4)) + list(range(5, read_pos-1 + 50))
            indel_index = indel_candidates[np.random.randint(len(indel_candidates))]
            indels.append((indel_type, indel_index, indel_len))
    print(indels)

    read_seq = local_hap[span+read_pos-1: span+read_pos+read_num-1]
    read_seq = ''
    prev_pos = span + read_pos - 1
    for indel in indels:
        current_pos = 0
        read_seq += local_hap[prev_pos, current_pos]
        prev_pos = current_pos

    ins_str = 'NONE'
    del_str = 'NONE'

    read_qual = generate_read_qual(is_read1)

    # BX:NONE;INS:-12(3),15(2);DEL:-86(1),36(4);
    read_meta_info = 'BX:{};INS:{};DEL:{}'.format(
        barcode if barcode else 'NONE',
        ins_str,
        del_str
    )
    res = '{}/{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        read_id,
        1 if is_read1 else 2,
        read_flag,
        read_pos,
        read_seq,
        read_qual,
        read_meta_info
    )
    return res


def generate_read_qual(is_read1):
    seed = np.random.randint(17)
    fn = '/home/chenlingxi/mnt/chenlingxi/workspace/Bio_Projects/ComplexSV/demo_data/read{}_qual'.format(
        1 if is_read1 else 2)
    for i, line in enumerate(open(fn, 'r')):
        if i == seed:
            return line.strip()


def simu_pe(local_hap, sv_type, span, barcode):
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

    read1_str = simu_single_read(local_hap, read_id, read1_flag, read1_pos, span, barcode, True)
    read2_str = simu_single_read(local_hap, read_id, read2_flag, read2_pos, span, barcode, False)
    return read1_str + read2_str


def reverse_complementary(seq):
    maps = {'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C',
            'N': 'N'}
    return ''.join([maps[x] for x in seq][::-1])


def simu_local_hap(span, ref_fasta, chrom_5p, bkpos_5p, chrom_3p, bkpos_3p, sv_type, inner_ins):
    if sv_type in ['ht', 'th']:
        start_5p = bkpos_5p - span
        end_5p = bkpos_5p
        start_3p = bkpos_3p - 1  # bed start from 0
        end_3p = bkpos_3p + span - 1
    elif sv_type == 'hh':
        start_5p = bkpos_5p - span
        end_5p = bkpos_5p
        start_3p = bkpos_3p - span
        end_3p = bkpos_3p
    elif sv_type == 'tt':
        start_5p = bkpos_5p - 1
        end_5p = bkpos_5p + span - 1
        start_3p = bkpos_3p - 1
        end_3p = bkpos_3p + span - 1

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
        local_hap = seq_5p + inner_ins + reverse_complementary(seq_3p)
    elif sv_type == 'tt':
        local_hap = reverse_complementary(seq_5p) + inner_ins + seq_3p

    return local_hap


def simu_reads(ref_fasta,
               chrom_5p, bkpos_5p, strand_5p,
               chrom_3p, bkpos_3p, strand_3p, sv_type,
               inner_ins, hm, junc_reads):
    span = 400
    local_hap = simu_local_hap(span, ref_fasta, chrom_5p, bkpos_5p, chrom_3p, bkpos_3p, sv_type, inner_ins)

    meta_info = 'Inner-ins:{};HM:{}'.format(
        inner_ins if inner_ins else 'NONE',
        hm if hm else 'NONE'
    )
    res = '##{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        chrom_5p, bkpos_5p, strand_5p,
        chrom_3p, bkpos_3p, strand_3p,
        meta_info, junc_reads
    )

    barcodes = generate_barcodes(junc_reads)
    for i in range(junc_reads):
        res += simu_pe(local_hap, sv_type, span, barcodes[i])
    return res


def generate_barcodes(num):
    barcodes = []
    seeds = random.sample(range(1000), num)
    for i, line in enumerate(open('/home/chenlingxi/mnt/chenlingxi/workspace/Bio_Projects/ComplexSV/demo_data/barcode_list', 'r')):
        if i not in seeds:
            continue
        barcodes.append(line.strip()+'-1')
    return barcodes


def main():
    ref_fasta = '/home/BIOINFO_DATABASE/reference/genome_DNA/Homo_sapiens/hg19/BWA_GATK_index/hg19.fa'

    # KIF5B(-) ALK(-) gene fusion
    sv_type = 'ht'
    chrom_5p = 'chr2'  # ALK
    bkpos_5p = 29592774
    # chr2:29591774-29593774
    strand_5p = '+'
    chrom_3p = 'chr10'  # KIF5B
    bkpos_3p = 32307938
    # chr10:32306938-32308938
    strand_3p = '+'
    inner_ins = 'CCGA'
    hm = ''
    junc_reads = random.randint(10, 20)
    res = simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                     chrom_3p, bkpos_3p, strand_3p, sv_type,
                     inner_ins, hm, junc_reads)
    # FGFR3 (+) TACC3 (+) gene fusion
    sv_type = 'th'
    chrom_5p = 'chr4'  # FGFR3
    bkpos_5p = 1808661
    strand_5p = '+'
    chrom_3p = 'chr4'  # TACC3
    bkpos_3p = 1739325
    strand_3p = '+'
    inner_ins = ''
    hm = 'CTG,1'
    junc_reads = random.randint(10, 20)
    res += simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                      chrom_3p, bkpos_3p, strand_3p, sv_type,
                      inner_ins, hm, junc_reads)

    # NPM1 (+) ALK (-)
    sv_type = 'hh'
    chrom_5p = 'chr5'  # NPM1
    bkpos_5p = 170824083
    strand_5p = '+'
    chrom_3p = 'chr2'  # ALK
    bkpos_3p = 29792804
    strand_3p = '-'
    inner_ins = 'TTGCAA'
    hm = 'GGT,-2'
    junc_reads = random.randint(10, 20)
    res += simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                      chrom_3p, bkpos_3p, strand_3p, sv_type,
                      inner_ins, hm, junc_reads)
    # PARG (-) BMS1 (+)
    sv_type = 'tt'
    chrom_5p = 'chr10'  # PARG
    bkpos_5p = 51093249
    strand_5p = '-'
    chrom_3p = 'chr10'  # BMS1
    bkpos_3p = 43287075
    strand_3p = '+'
    inner_ins = ''
    hm = ''
    junc_reads = random.randint(10, 20)
    res += simu_reads(ref_fasta, chrom_5p, bkpos_5p, strand_5p,
                      chrom_3p, bkpos_3p, strand_3p, sv_type,
                      inner_ins, hm, junc_reads)

    print(res)


if __name__ == "__main__":
    main()
