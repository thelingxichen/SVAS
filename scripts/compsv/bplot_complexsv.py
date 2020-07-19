#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Usage:
    complexsv call --sample=STR --region_fn=IN_FILE --sv_fn=IN_FILE
    complexsv draw --in_fn=IN_FILE --out_dir=DIR
    complexsv -h | --help

Options:
    --out_dir=DIR

    -v | --verbose	Show log information on console.
    -h | --help		Show this screen.

Others:
    tenxtools.complexsv
    ~~~~~~~~~~~~~~~~~~~~~~~~

    @Input: DESCRIPTION
    @Return:

    @Copyright: (c) 2018-09 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""

from biotool.bplot import svgfigure
from biotool.bplot import svgobject
from biotool import myio


class SVJoin(svgobject.SVGObject):

    def __init__(self, insert, size,
                 cut_y=0,
                 join_type=None, arc_v=None,
                 precise=None,
                 large_sv=None,
                 enable_horizontal=None,
                 **extra):
        self.cut_y = cut_y
        self.join_type = join_type
        if size[0] == 0 and enable_horizontal:
            self.arc_h = arc_v or (size[1] - cut_y)/4 or abs(size[1]/4)
            self.arc_v = 0
        else:
            self.arc_v = arc_v or (size[1] - cut_y)/4 or abs(size[0]/8)
            self.arc_h = 0

        self.precise = precise
        self.large_sv = large_sv
        super(SVJoin, self).__init__(insert, size, **extra)

    def render(self):
        '''
            a       b
            ---------
            |       |
            |       |
            ---------
            c       d
                or
            a
            |       |   cut_y
            c       |
            |       b
            |       |
                    d
                or   horizontal presentation
            a
            |
            c
            |  cut_y
            |

            |  cut_y
            |
            b
            |
            d
        '''
        a = self.insert
        b = self.insert[0] + self.size[0], self.insert[1] + self.cut_y
        c = self.insert[0], self.insert[1] + self.size[1] - self.cut_y
        d = self.insert[0] + self.size[0], self.insert[1] + self.size[1]
        v = self.size[1] - self.cut_y
        if self.join_type == 'ht':
            o = c
            bp1 = a[0] + self.arc_h, a[1] + self.arc_v
            bp2 = b[0] + self.arc_h, b[1] + self.arc_v
            fp = b
            color = 'blue'
        elif self.join_type == 'th':
            o = c
            bp1 = a[0] - self.arc_h, a[1] - self.arc_v
            bp2 = b[0] - self.arc_h, b[1] - self.arc_v
            fp = b
            color = 'red'
        elif self.join_type == 'tt':
            o = c
            bp1 = a[0] - self.arc_h, a[1] - self.arc_v
            bp2 = b[0] - self.arc_h, b[1] - self.arc_v
            fp = b
            color = 'green'
        elif self.join_type == 'hh':
            o = c
            bp1 = a[0] + self.arc_h, a[1] + self.arc_v
            bp2 = b[0] + self.arc_h, b[1] + self.arc_v
            fp = b
            color = 'orange'
        else:
            o = c
            bp1 = a[0], a[1]
            bp2 = b[0], b[1]
            fp = b
            color = 'black'

        if self.large_sv:
            opacity = 1
        else:
            opacity = 0.2

        path_desc = "M {},{} v {} C {},{} {},{} {},{} v {}".format(o[0], o[1], -v,
                                                                   bp1[0], bp1[1],
                                                                   bp2[0], bp2[1],
                                                                   fp[0], fp[1], v)
        sv_join = self._dwg.path(d=path_desc,
                                 fill='none', stroke=color, opacity=opacity,
                                 stroke_dasharray=self.stroke_dasharray)

        if not self.precise:
            rect_5p = svgobject.Rect(c, (5, 5), fill=color, stroke='none', opacity=opacity,
                                     self_pos='middle:up').render()
            rect_3p = svgobject.Rect(d, (5, 5), fill=color, stroke='none', opacity=opacity,
                                     self_pos='middle:up').render()
        else:
            rect_5p = None
            rect_3p = None

        return sv_join, rect_5p, rect_3p


class RegionGroupsPlot(svgfigure.SVGFigure):

    def __init__(self, out_fn=None, sample=None, groups=None,
                 resolution=None,
                 sv_join_width=None,
                 **kwargs):
        self.sample = sample
        self.groups = groups
        self.figsize = kwargs['size']

        self.sv_join_width = sv_join_width or 100

        super().__init__(fn=out_fn, **kwargs)

    def _measure_group(self, group, i):
        length = 0
        region_num = 0
        for region in group.region_list:
            length += (region.end - region.start) + 1
            region_num += 1

        group.resolution = length/(self.figsize[0]*0.4)
        min_gap = 15
        for region in group.region_list:
            d = region.min_bkp_distance
            if d and (region.end - region.start)/float(d) < self.figsize[0]*0.02/float(region_num):
                region.resolution = d / min_gap
            else:
                region.resolution = group.resolution

        group_width = 200
        group_start = 150
        region_width = 150
        region_delta_v = 20
        prev_region = None
        for region in group.region_list:
            delta_x = (region.end - region.start) / float(region.resolution)
            if not prev_region:
                region.insert = (group_start,  group_width*i + region_width)
                region.size = (delta_x, region_width)
            else:
                region.insert = (prev_region.insert[0] + prev_region.size[0] + region_delta_v,
                                 group_width*i + region_width)
                region.size = (delta_x, region_width)
            prev_region = region

    def draw(self):
        for i, group in enumerate(self.groups):
            meta, group = group
            self._measure_group(group, i)
            self._draw_group(group)

        # default painting
        self._draw_info()
        self._draw_icon()

        self.dwg.add(self.canvas)
        self.dwg.save()

    def _draw_group(self, group):

        for i, region in enumerate(group.region_list):
            # print(region, region.insert, region.size)
            region_rect = self._draw_region(region)
            if i == 0:
                group_label = svgobject.Text(group.group_name, font_size=10, fill='black',
                                             anchor_object=region_rect, related_insert='left:center').render()
                self.canvas.add(group_label)

            for i, sv in enumerate(region.sv_list):
                bs, bs_5p, bs_3p = group.has_sv(sv)
                # print(sv.chrom_5p, sv.bkpos_5p/linkage_distance, sv.chrom_3p, sv.bkpos_3p/linkage_distance, sv.orientation)
                is_5p = [i for i, b in enumerate(bs_5p) if b]
                is_3p = [i for i, b in enumerate(bs_3p) if b]
                if len(is_5p) > 1 or len(is_3p) > 1 or len(is_5p) == 0 or len(is_3p) == 0:   # check this bug
                    print('Error is_5p > 1 or is_3p > 1')
                    print([group.region_list[i] for i in is_5p])
                    print([group.region_list[i] for i in is_3p])
                    print(sv.chrom_5p, sv.bkpos_5p, sv.strand_5p)
                    print(sv.chrom_3p, sv.bkpos_3p, sv.strand_3p)
                    continue

                region_5p = group.region_list[is_5p[0]]
                region_3p = group.region_list[is_3p[0]]
                self._draw_sv_join(region_5p, region_3p, sv)
                if sv.tran_5p and sv.tran_3p and sv.tran_5p.gene == sv.tran_3p.gene:
                    # draw one time
                    self._draw_gene(region_5p, sv.tran_5p, self.sv_join_width)
                elif sv.tran_5p:
                    self._draw_gene(region_5p, sv.tran_5p, self.sv_join_width)
                elif sv.tran_3p:
                    self._draw_gene(region_3p, sv.tran_3p, self.sv_join_width)

            for meta, enhancer in region.enhancer_list:
                self._draw_enhancer(region, enhancer, 'purple', self.sv_join_width + 25)

            for meta, enhancer in region.super_enhancer_list:
                self._draw_enhancer(region, enhancer, 'red', self.sv_join_width + 35)

            for i, cn_region in enumerate(region.cn_list):
                self._draw_cn(region, cn_region, self.sv_join_width)

    def _draw_region(self, region):
        insert, size = region.insert, region.size
        rect = svgobject.Rect(insert, size, stroke='gray', fill='none').render()
        label = svgobject.Text(region.chrom, font_size=10, fill='black',
                               anchor_object=rect, related_insert='middle:up').render()

        start = svgobject.Text(region, font_size=10, fill='black',
                               anchor_object=rect, related_insert='middle:bottom').render()

        self.canvas.add(rect)
        self.canvas.add(label)
        self.canvas.add(start)

        return rect

    def _draw_sv_join(self, region_5p, region_3p, sv):
        bkpos_5p = (sv.bkpos_5p - region_5p.start) / float(region_5p.resolution)
        bkpos_3p = (sv.bkpos_3p - region_3p.start) / float(region_3p.resolution)

        x_5p, y_5p = region_5p.insert
        x_3p, y_3p = region_3p.insert
        x_5p += bkpos_5p
        cut_y = y_3p - y_5p
        cut_x = x_3p - x_5p
        insert = (x_5p, y_5p)
        size = (cut_x + bkpos_3p, cut_y + self.sv_join_width)

        if sv.chrom_5p == sv.chrom_3p and abs(sv.bkpos_5p - sv.bkpos_3p) < 1000:
            large_sv = False
        else:
            large_sv = True
        sv_join, rect_5p, rect_3p = SVJoin(insert, size, cut_y=cut_y, join_type=sv.orientation,
                                           precise=sv.precise, large_sv=large_sv).render()
        self.canvas.add(sv_join)
        if rect_5p:
            self.canvas.add(rect_5p)
        if rect_3p:
            self.canvas.add(rect_3p)

    def _draw_gene(self, region, tran, down_y=0):
        if tran.group:
            color = 'red'
        else:
            color = 'black'
        x, y = region.insert
        y = y + down_y
        delta_x, delta_y = region.size
        gene_start = (tran.left - region.start) / float(region.resolution)
        gene_end = (tran.right - region.start) / float(region.resolution)

        if gene_start < 0:
            gene_start = 0
        if gene_end > delta_x:
            gene_end = delta_x

        insert = x + gene_start, y
        size = gene_end - gene_start, 5

        gene_object = svgobject.Rect(insert, size, stroke='none', fill=color).render()
        label = '{}{}:{}'.format(tran.strand, tran.gene, tran.func_loc)
        gene_label = svgobject.Text(label, font_size=10, font_style='italic', fill=color,
                                    anchor_object=gene_object, related_insert='middle:bottom').render()

        self.canvas.add(gene_object)
        self.canvas.add(gene_label)

    def _draw_enhancer(self, region, enhancer, color='purple', down_y=0):
        x, y = region.insert
        delta_x, delta_y = region.size
        y = y + down_y
        start = (enhancer.start - region.start) / float(region.resolution)
        end = (enhancer.end - region.start) / float(region.resolution)

        if start < 0:
            start = 0
        if end > delta_x:
            end = delta_x

        insert = x + start, y
        size = end - start, 5
        enhancer = svgobject.Rect(insert, size, stroke='none', fill=color).render()

        self.canvas.add(enhancer)

    def _draw_cn(self, region, cn_region, down_y=0):
        cn = cn_region.cn
        balance_cn = round(cn)
        x, y = region.insert
        delta_x, delta_y = region.size
        y = y + down_y
        start = (cn_region.start - region.start) / float(region.resolution)
        end = (cn_region.end - region.start) / float(region.resolution)

        if start < 0:
            start = 0
        if end > delta_x:
            end = delta_x

        size = end - start, 0

        cn_resolution = 5
        insert = x + start, y - cn * cn_resolution
        cn_line = svgobject.Line(insert, size, stroke='gray', fill='none').render()
        insert = x + start, y - balance_cn * cn_resolution
        balance_cn_line = svgobject.Line(insert, size, stroke='black', fill='none').render()
        cn_label = svgobject.Text(balance_cn, font_size=10, fill='black',
                                  anchor_object=balance_cn_line, related_insert='middle:bottom').render()
        self.canvas.add(cn_line)
        self.canvas.add(balance_cn_line)
        self.canvas.add(cn_label)


def run_call(sample=None, region_fn=None, sv_fn=None, **args):
    regions = []
    for line in myio.safe_open(region_fn, 'r'):
        region = line.strip().split('\t')[0]
        chrom = region.split(':')[0].replace('p', '').replace('q', '')
        try:
            start, end = map(int, region.split(':')[1].split('-'))
            regions.append((chrom, start, end))
        except Exception:
            continue

    # res = {}
    '''
    for i, record in enumerate(svaba.SVReader(in_fn=sv_fn, sample=sample)):
        meta, record = record
    '''
    # run_draw(region_fn, records)


def run_draw(groups=None, out_dir=None, **args):
    title = 'Join Type'
    groups = list(groups)
    group_width = 250
    size = (2500, group_width * len(groups) + 200)
    plot = RegionGroupsPlot(groups=groups, out_dir=out_dir, title=title, size=size,
                            group_width=group_width,
                            **args)
    plot.draw()


def run(call=None, draw=None, **args):
    if call:
        run_call(**args)
    if draw:
        run_draw(**args)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    run(**args)
