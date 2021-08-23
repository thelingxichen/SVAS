import sv

def read_seg_fn(seg_fn):
    seg2cn = {} 
    if not seg_fn:
        return seg2cn

    for line in open(seg_fn, 'r'):
        splits = line.strip().split()
        avg_cn = float(splits[1])
        seg = splits[0]
        chrom = seg.split(':')[0]
        start = seg.split(':')[1].split('-')[0]
        end = seg.split(':')[1].split('-')[1]

        key = (chrom, start, end)
        seg2cn[key] = avg_cn

    return seg2cn

def get_segments_from_sv_group(group, contig2length, margin=100000):

    contig2bp = {}

    for r in group:
        if r.chrom_5p in contig2bp:
            contig2bp[r.chrom_5p].add(r.bkpos_5p)
        else:
            contig2bp[r.chrom_5p] = set([r.bkpos_5p])
        if r.chrom_3p in contig2bp:
            contig2bp[r.chrom_3p].add(r.bkpos_3p)
        else:
            contig2bp[r.chrom_3p] = set([r.bkpos_3p])

    contig2seg = {}
    contig_seg_list = []
    bp2seg = {}
    for contig in contig2bp.keys():

        bp_list = sorted(list(contig2bp[contig]))
        # append source segment
        if bp_list[0] != 0:
            if bp_list[0] - margin < 0:
                bp_list = [1] + bp_list
            else:
                bp_list = [bp_list[0] - margin] + bp_list

        length = contig2length.get(contig)
        if length and bp_list[-1] != length:
            if bp_list[-1] + margin > length:
                bp_list.append(length)
            else:
                bp_list.append(bp_list[-1] + margin)
        else:
            bp_list.append(bp_list[-1] + margin)

        seg_list = list(zip(bp_list[:-1], bp_list[1:], [None]*len(bp_list)))

        for i, bp in enumerate(bp_list):
            if i == 0:
                bp2seg[(contig, bp)] = (None, seg_list[i])
            elif i == len(bp_list)-1:
                bp2seg[(contig, bp)] = (seg_list[i-1], None)
            else:
                bp2seg[(contig, bp)] = (seg_list[i-1], seg_list[i])
        contig_seg_list.append([contig, bp_list, seg_list])
        contig2seg[contig] = seg_list

    return contig_seg_list, contig2seg, bp2seg


def encode_segments(contig2seg, seg2cn):
    res = ''
    segid2seg = {}
    seg2segid = {}
    count = 1
    for contig, _, seg_list in contig2seg:
        for start, end, avg_cn in seg_list:
            length = end - start + 1
            seg_id = 'H:{}'.format(count)
            seg = '{}_{}_{}'.format(contig, start, end)
            segid2seg[seg_id] = seg
            seg2segid[seg] = seg_id
            avg_cn = seg2cn.get((contig, start, end), avg_cn)
            res += 'SEG {}:{}:{}:{} {} {}\n'.format(
                seg_id, seg, 1, length, 30, avg_cn or -1
            )
            count += 1

    return res, seg2segid, segid2seg


def encode_sv(sv_records, bp2seg, seg2segid, contig2length, contig2seg):
    # JUNC H:144:+ H:145:+ 382 -1 U B
    res = ''
    for r in sv_records:
        if r.strand_5p == '-' and r.strand_3p == '-':
            tmp = r.chrom_5p
            r.chrom_5p = r.chrom_3p
            r.chrom_3p = tmp
            tmp = r.bkpos_5p
            r.bkpos_5p = r.bkpos_3p
            r.bkpos_3p = tmp
            r.strand_5p = '+'
            r.strand_3p = '+'

        bp_5p = (r.chrom_5p, r.bkpos_5p)
        bp_3p = (r.chrom_3p, r.bkpos_3p)

        # add sv
        if r.strand_5p == '+':
            seg_5p = '{}_{}_{}'.format(r.chrom_5p, bp2seg[bp_5p][0][0], bp2seg[bp_5p][0][1])
        else:
            seg_5p = '{}_{}_{}'.format(r.chrom_5p, bp2seg[bp_5p][1][0], bp2seg[bp_5p][1][1])

        if r.strand_3p == '+':
            seg_3p = '{}_{}_{}'.format(r.chrom_3p, bp2seg[bp_3p][1][0], bp2seg[bp_3p][1][1])
        else:
            seg_3p = '{}_{}_{}'.format(r.chrom_3p, bp2seg[bp_3p][0][0], bp2seg[bp_3p][0][1])

        res += 'JUNC {}:{} {}:{} 30 {} U B\n'.format(
            seg2segid[seg_5p], r.strand_5p,
            seg2segid[seg_3p], r.strand_3p,
            r.avg_cn or -1
        )
        if r.chrom_5p == r.chrom_3p:
            continue

        # add virtual edge for translocation
        join_type = sv.get_orientation_by_strand(r, list(contig2length.keys()))

        if join_type in ['ht', 'th']:
            # add 3pt+ -> 5ph+
            seg_5p = contig2seg[r.chrom_3p][-1]
            seg_3p = contig2seg[r.chrom_5p][0]
            seg_5p = '{}_{}_{}'.format(r.chrom_3p, seg_5p[0], seg_5p[1])
            seg_3p = '{}_{}_{}'.format(r.chrom_5p, seg_3p[0], seg_3p[1])
            strand_5p = '+'
            strand_3p = '+'

        if join_type == 'hh':
            # add 3ph- 5ph+
            seg_5p = contig2seg[r.chrom_3p][0]
            seg_3p = contig2seg[r.chrom_5p][0]
            seg_5p = '{}_{}_{}'.format(r.chrom_3p, seg_5p[0], seg_5p[1])
            seg_3p = '{}_{}_{}'.format(r.chrom_5p, seg_3p[0], seg_3p[1])
            strand_5p = '-'
            strand_3p = '+'

        if join_type == 'tt':
            # add 3pt+ 3pt-
            seg_5p = contig2seg[r.chrom_3p][-1]
            seg_3p = contig2seg[r.chrom_5p][-1]
            seg_5p = '{}_{}_{}'.format(r.chrom_3p, seg_5p[0], seg_5p[1])
            seg_3p = '{}_{}_{}'.format(r.chrom_5p, seg_3p[0], seg_3p[1])
            strand_5p = '+'
            strand_3p = '-'

        res += 'JUNC {}:{} {}:{} 30 {} I U\n'.format(
            seg2segid[seg_5p], r.strand_5p,
            seg2segid[seg_3p], r.strand_3p,
            r.avg_cn or -1
        )

    # add normal edges
    for contig in contig2seg:
        seg_list = contig2seg[contig]
        for seg1, seg2 in zip(seg_list[:-1], seg_list[1:]):
            seg_5p = '{}_{}_{}'.format(contig, seg1[0], seg1[1])
            seg_3p = '{}_{}_{}'.format(contig, seg2[0], seg2[1])
            res += 'JUNC {}:{} {}:{} 30 {} I U\n'.format(
                seg2segid[seg_5p], '+',
                seg2segid[seg_3p], '+',
                r.avg_cn or -1
            )

        seg_5p = '{}_{}_{}'.format(contig, seg_list[-1][0], seg_list[-1][1])
        seg_3p = '{}_{}_{}'.format(contig, seg_list[0][0], seg_list[0][1])
        res += 'JUNC {}:{} {}:{} 30 {} I U\n'.format(
            seg2segid[seg_5p], '+',
            seg2segid[seg_3p], '+',
            r.avg_cn or -1
        )

    return res

def encode_graph(group_i, group, contig_seg_list, bp2seg, contig2length, contig2seg, seg2cn):
    seg_num = len([s for _, _, seg_list in contig_seg_list for s in seg_list])
    res = '''
SAMPLE group{}
AVG_SEG_DP 30
AVG_JUNC_DP 30
PURITY 1
AVG_PLOIDY 1
PLOIDY 1m1
SOURCE H:1
SINK H:{}
'''.format(group_i, seg_num)
    seg_str, seg2segid, segid2seg = encode_segments(contig_seg_list, seg2cn)
    res += seg_str
    res += encode_sv(group, bp2seg, seg2segid, contig2length, contig2seg)

    return res


def encode_junc_db(group, contig2length):
    res = 'chrom_5p\tpos_5p\tstrand_5p\tchrom_3p\tpos_3p\tstrand_3p\tcount\ttype\tgene_5p\tgene_3p\tclone\n'
    for r in group:
        join_type = sv.get_orientation_by_strand(r, list(contig2length.keys()))
        res += '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            r.chrom_5p, r.bkpos_5p, r.strand_5p,
            r.chrom_3p, r.bkpos_3p, r.strand_3p,
            30, r.type,
            r.gene_5p, r.gene_3p, r.clone
        )
    return res
