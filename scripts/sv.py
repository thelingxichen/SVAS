import copy
import vcf
import collections

from bioutensil import myio


class SVRecord(myio.Record):
    fields = 'chrom_5p,bkpos_5p,strand_5p,chrom_3p,bkpos_3p,strand_3p,'
    fields += 'inner_ins,span_reads,junc_reads,id,qual,filter,'
    fields += 'group,meta_info,gene_5p,gene_3p,clone,genes,type,avg_cn'
    fields = fields.split(',')

    def __init__(self, parent=None, tool='svaba', *args, **kwargs):
        super(SVRecord, self).__init__(*args, **kwargs)
        self.parent = parent
        self.tool = tool
        if parent:
            self.from_parent()

    def __repr__(self):
        funcs = (self._format_info if 'info' in field else self._format_value
                 for field in self.fields)
        return '\t'.join(func(getattr(self, field))
                         for field, func in zip(self.fields, funcs))

    def _format_value(self, x):  # None, int, float, list
        if x == 0:
            return '0'
        if not x:
            return '.'   # notice not 0 is True
        if type(x) == list:
            return self._format_list(x)
        return str(x)

    def _format_list(self, x, sep=','):
        return sep.join(map(self._format_value, x))

    def _format_pair(self, pair):
        k, v = pair
        if v is True:
            return '{}'.format(self._format_value(k))
        formated_v = self._format_value(v)
        if formated_v == '.':
            return ''
        return '{}={}'.format(self._format_value(k), formated_v)

    def _format_info(self, info):
        if not info:
            return '.'
        xs = (self._format_pair(x) for x in info.items())
        return ';'.join(x for x in xs if x)

    def from_parent(self):
        args = {}
        args['meta_info'] = self._get_sv_meta_info
        if 'BND' in self.var_type:
            args.update(self._get_bnd_sv_fields)
        elif self.var_type == 'INS':
            args.update(self._get_ins_sv_fields)
        elif self.var_type == 'DEL':
            args.update(self._get_del_sv_fields)
        elif self.var_type == 'DUP':
            args.update(self._get_dup_sv_fields)

        args.update(self._get_general_sv_fields)
        self.set(**args)

    def set(self, group=None, chrom_5p=None, bkpos_5p=None, strand_5p=None,
            chrom_3p=None, bkpos_3p=None, strand_3p=None,
            inner_ins=None, span_reads=None, junc_reads=None,
            id=None, qual=None, filter=None,
            meta_info=None, anno_info=None,
            gene_5p=None, gene_3p=None, clone=None, genes=None, type=None,
            avg_cn=None,
            sample=None):
        self.gene_5p = self._validate(gene_5p)
        self.chrom_5p = self._validate(chrom_5p)
        self.bkpos_5p = self._validate(bkpos_5p, int)
        self.strand_5p = self._validate(strand_5p)
        self.gene_3p = self._validate(gene_3p)
        self.chrom_3p = self._validate(chrom_3p)
        self.bkpos_3p = self._validate(bkpos_3p, int)
        self.strand_3p = self._validate(strand_3p)
        self.inner_ins = self._validate(inner_ins)
        self.span_reads = self._validate(span_reads, int)
        self.junc_reads = self._validate(junc_reads, int)
        self.group = self._validate(group)
        self.id = self._validate(id)
        self.qual = self._validate(qual, float)
        self.filter = self._validate(filter, self._parse_list)
        self.meta_info = self._validate(meta_info, self._parse_info)
        self.anno_info = self._validate(anno_info, self._parse_info)
        self.clone = self._validate(clone)
        self.genes = self._validate(genes)
        self.type = self._validate(type)
        self.avg_cn = self._validate(avg_cn, float)

        self.sample = sample

    def _validate(self, x, func=None):
        if func in [None, int, float]:
            if not x:
                return x
            if x == '.':
                return None
        if func is str:
            if not x:
                return ''
        if func:
            return func(x)
        return x

    def parse(self, line):  # parse a stream_or_string to sv object
        args = dict(zip(self.fields, line.strip().split('\t')))
        self.set(**args)

    def _parse_info(self, info):
        if not info:
            return {}
        if type(info) == dict:
            return info
        # parse string
        if info == '.':
            return {}
        pairs = []
        for pair in info.split(';'):
            pair = self._parse_pair(pair)
            pairs.append(pair)

        return dict(pairs)

    def _parse_list(self, line, sep=','):
        if isinstance(line, list):
            return line
        if not line or line == '.':
            return []
        return line.split(sep)

    def _parse_value(self, x):  # None, int, float, list
        if x == '.':
            return None
        splits = x.split(',')
        if splits:
            return map(self._parse_value, splits)
        return float(x)

    @property
    def _get_general_sv_fields(self):
        args = {}

        args['id'] = self.parent.ID.replace(':', '_')
        args['qual'] = self.parent.QUAL
        args['filter'] = self.parent.FILTER
        args['inner_ins'] = self.parent.INFO.get('INSERTION', None)
        if self.tool == 'svaba':
            args['span_reads'] = [call['DR'] for call in self.parent.samples][0]
            args['junc_reads'] = [call['SR'] for call in self.parent.samples][0]
        if self.tool == 'manta':
            if self.parent.samples:
                ref, alt = [call['PR'] for call in self.parent.samples][0]
                args['span_reads'] = ref + alt
          
                ref, alt = [call.data._asdict().get('SR', (0, 0)) for call in self.parent.samples][0]
                args['junc_reads'] = ref + alt
            else:
                args['span_reads'] = self.parent.INFO.get('BND_PAIR_COUNT', None) 
                args['junc_reads'] = self.parent.INFO.get('BND_PAIR_COUNT', None) 

        return args

    @property
    def _get_sv_meta_info(self):
        meta_info = copy.deepcopy(self.parent.INFO)

        # add new tags
        meta_info['TOOL'] = self.tool 
        if self.tool == 'svaba':
            meta_info['MATEID'] = meta_info['MATEID'].replace(':', '_')
        if self.tool == 'manta':
            if 'MATEID' in meta_info.keys(): 
                meta_info['MATEID'] = meta_info['MATEID'][0].replace(':', '_')
            else:
                meta_info['MATEID'] = None

        genotype = []
        for call in self.parent.samples:
            g = call.sample.split('/')[-1] + '(' + ' '.join('{}:{}'.format(x, y) for x, y in zip(call.data._fields, call.data)) + ')'
            genotype.append(g)
        meta_info['GENOTYPE'] = genotype
        meta_info['ALT'] = self.parent.ALT

        meta_info['VARTYPE'] = self.var_type
        if 'BND' in self.var_type:
            meta_info['JOINTYPE'] = self.orientation
        elif self.var_type == 'INS':
            meta_info['JOINTYPE'] = 'ins'
        elif self.var_type == 'DEL':
            meta_info['JOINTYPE'] = 'ht'
        elif self.var_type == 'DUP':
            meta_info['JOINTYPE'] = 'th'
        else:
            print('vartype', self.var_type)

        return meta_info

    @property
    def _get_ins_sv_fields(self):
        args = {}
        args['strand_5p'] = '+'
        args['strand_3p'] = '+'
        args['chrom_5p'] = self.parent.CHROM 
        args['bkpos_5p'] = self.parent.POS
        args['chrom_3p'] = self.parent.CHROM
        args['bkpos_3p'] = self.parent.INFO['END']
        return args

    @property
    def _get_del_sv_fields(self):
        args = {}
        args['strand_5p'] = '+'
        args['strand_3p'] = '+'
        args['chrom_5p'] = self.parent.CHROM 
        args['bkpos_5p'] = self.parent.POS
        args['chrom_3p'] = self.parent.CHROM
        args['bkpos_3p'] = self.parent.INFO['END']
        return args

    @property
    def _get_dup_sv_fields(self):
        args = {}
        args['strand_5p'] = '+'
        args['strand_3p'] = '+'
        args['chrom_5p'] = self.parent.CHROM
        args['bkpos_5p'] = self.parent.INFO['END']
        args['chrom_3p'] = self.parent.CHROM 
        args['bkpos_3p'] = self.parent.POS
        return args

    @property
    def _get_bnd_sv_fields(self):
        args = {}
        alt = self.parent.ALT[0]

        # if strand1 == '-' and strand2 == '-':   # inversion
        if alt.remoteOrientation and alt.orientation:
            # [p[t tt-inversion
            args['strand_5p'] = '-'
            args['strand_3p'] = '+'
        # elif strand1 == '+' and strand2 == '+':  # inversion
        elif not alt.remoteOrientation and not alt.orientation:
            # t]p] hh-inversion
            args['strand_5p'] = '+'
            args['strand_3p'] = '-'
        else:                                   # not inversion
            args['strand_5p'] = '+'
            args['strand_3p'] = '+'

        if alt.chr == 'MT':
            alt.chr == 'chrM'

        number_chrs = list(map(str, range(1, 23))) + ['X', 'Y']
        if alt.chr in number_chrs:
            alt.chr = 'chr' + alt.chr

        if self.parent.CHROM == 'MT':
            self.parent.CHROM == 'chrM'
        if self.parent.CHROM in number_chrs:
            self.parent.CHROM = 'chr' + self.parent.CHROM

        if not alt.remoteOrientation and alt.orientation:
            # t[p[ ht
            args['chrom_3p'] = self.parent.CHROM
            args['bkpos_3p'] = self.parent.POS
            args['chrom_5p'] = alt.chr
            args['bkpos_5p'] = alt.pos

            args['inner_ins'] = alt.connectingSequence[::-1]

        else:
            # ]p]t  th
            args['chrom_5p'] = self.parent.CHROM
            args['bkpos_5p'] = self.parent.POS
            args['chrom_3p'] = alt.chr
            args['bkpos_3p'] = alt.pos

            args['inner_ins'] = alt.connectingSequence[1:]

        return args

    def __eq__(self, other, allowance=10):
        return self.var_type == other.var_type \
            and self.chrom_5p == other.chrom_5p \
            and self._pos_eq(self.bkpos_5p, other.bkpos_5p, allowance) \
            and self.strand_5p == other.strand_5p \
            and self.chrom_3p == other.chrom_3p \
            and self._pos_eq(self.bkpos_3p, other.bkpos_3p, allowance) \
            and self.strand_3p == other.strand_3p

    @staticmethod
    def _pos_eq(self, self_pos, other_pos, allowance=10):
        return self_pos >= other_pos - allowance \
            and self_pos <= other_pos + allowance

    @property
    def precise(self):
        if self.junc_reads >= 10:
            return True
        else:
            return False

    @property
    def orientation(self):
        '''
        update orientation in 2020.02.13,
        chrom bkpos strand result is correct, won't affect localhap
        '''
        '''
        # If the breakend is connected to sequence 3', extending to 5'
        breakend is head (based on + strand): alt.orientation == True
        # If the breakend is connected to sequence 5', extending to 3'
        breakend is tail (based on - strand): alt.orientation == False
        '''
        if not self.parent:
            return None

        alt = self.parent.ALT[0]
        if alt.remoteOrientation and alt.orientation:
            o = 'tt'
        elif not alt.remoteOrientation and not alt.orientation:
            o = 'hh'
        else:
            o = 'h' if alt.remoteOrientation else 't'   # self.parent.POS point
            o += 'h' if alt.orientation else 't'        # breakend point

        return o

    @property
    def bkpos_5p_orientation(self):
        return self.orientation[0]

    @property
    def bkpos_3p_orientation(self):
        return self.orientation[1]

    @property
    def sv_type(self):
        if self.parent.is_snp:
            return "snp"
        elif self.parent.is_indel:
            return "indel"
        elif self.parent.is_sv:
            return self.parent.INFO.get('SVTYPE', 'sv')
        else:
            return "unknown"

    @property
    def sv_type2(self):
        orientation = self.orientation
        if self.parent.ALT[0].chr != self.parent.CHROM:
            return 'TRX-'+orientation
        elif orientation == 'ht':
            return 'DEL-'+orientation
        elif orientation == 'th':
            return 'DUP-'+orientation
        else:   # 'hh' or 'tt'
            return 'INV-'+orientation

    @property
    def var_type(self):
        try:
            if self.sv_type == 'BND' and self.sv_type2:
                return 'BND:' + self.sv_type2
            else:
                return self.sv_type
        except Exception:
            return self.meta_info.get('VARTYPE', None)


def read_vcf(vcf_fn, tool='svaba', precise=True, id2genes={}, barcode2clone={}, **args):
    for vcf_record in vcf.Reader(filename=vcf_fn):
        sv_record = SVRecord(vcf_record, tool)
        if sv_record.id.endswith('_2'):
            continue
        if precise and (not sv_record.precise):
            continue

        id = sv_record.id.split('_')[0]
        if id in id2genes:
            sv_record.gene = list(id2genes[id])
        else:
            sv_record.gene = []

        if barcode2clone:
            add_barcode_clone(sv_record, barcode2clone)

        yield sv_record


def add_barcode_clone(sv_record, barcode2clone):
    bcs = sv_record.meta_info.get('BX', [])
    clones = [barcode2clone.get(bc.split('_')[0], '') for bc in bcs]
    sv_record.meta_info['CLONE'] = clones
    clone_prev = collections.Counter(c for c in clones if c)
    clone_prev = {c: f for c, f in clone_prev.items()}
    sv_record.meta_info['CLONE_PREV'] = clone_prev


def read_sv_group(sv_fn, group_header):
    reader = myio.Reader(in_fn=sv_fn, record_cls=SVRecord, sep='\t', has_header=True)
    return reader.read_chunks(group_header)


def get_orientation_by_strand(r, contig_list):
    if r.strand_5p == '+' and r.strand_3p == '-':
        return 'hh'
    if r.strand_5p == '-' and r.strand_3p == '+':
        return 'tt'
    if r.strand_5p == '+' and r.strand_3p == '+':
        if r.chrom_5p == r.chrom_3p:
            if r.bkpos_5p < r.bkpos_3p:
                return 'ht'
            else:
                return 'th'
        elif contig_list.index(r.chrom_5p) < contig_list.index(r.chrom_3p):
            return 'ht'
        else:
            return 'th'
    if r.strand_5p == '-' and r.strand_3p == '-':
        if r.chrom_5p == r.chrom_3p:
            if r.bkpos_5p < r.bkpos_3p:
                return 'th'
            else:
                return 'ht'
        elif contig_list.index(r.chrom_5p) < contig_list.index(r.chrom_3p):
            return 'th'
        else:
            return 'ht'

    return None
