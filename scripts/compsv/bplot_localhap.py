#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Usage:
    localhap call --sample=STR --region_fn=IN_FILE --sv_fn=IN_FILE
    localhap draw --in_fn=IN_FILE --out_dir=DIR
    localhap -h | --help

Options:
    --out_dir=DIR

    -v | --verbose	Show log information on console.
    -h | --help		Show this screen.

Others:
    tenxtools.localhap
    ~~~~~~~~~~~~~~~~~~~~~~~~

    @Input: DESCRIPTION
    @Return:

    @Copyright: (c) 2018-09 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
import svgwrite

from tenxtools.mydocopt import docopt

from tenxtools.bplot.utils import svgfigure
from tenxtools.bplot.utils import svgobject
from tenxtools.bplot import complexsv


class GraphsPlot(svgfigure.SVGFigure):

    def __init__(self, out_fn=None, sample=None, graphs=None,
                 resolution=None,
                 **kwargs):
        self.sample = sample
        self.graphs = graphs
        self.figsize = kwargs['size']

        self.vertex_width = 10

        super().__init__(fn=out_fn, **kwargs)

    def _measure_graph(self, graph, prev_graph):

        width, height = 50, 40

        graph.graph_x = 50
        if prev_graph:
            prev_groups = len(prev_graph.mGroupPriority)
            graph.graph_y = prev_graph.graph_y + prev_groups * (100 + height)
        else:
            graph.graph_y = 100

        for seg in graph.mSegments:
            i = graph.mGroupPriority.index(seg.mGroup)
            j = seg.mIdx

            seg.insert = graph.graph_x + j*100, graph.graph_y + i*100
            seg.size = width, height

    def draw(self):
        prev_graph = None
        for i, graph in enumerate(self.graphs):
            meta, graph = graph
            self._measure_graph(graph, prev_graph)
            self._draw_graph(graph)
            prev_graph = graph

        # default painting
        self._draw_info()
        self._draw_icon()

        self.dwg.add(self.canvas)
        self.dwg.save()

    def _draw_graph(self, graph):
        for i, seg in enumerate(graph.mSegments):
            self._draw_seg(seg)
            self._draw_seg_strand(seg, strand='+')
            self._draw_seg_strand(seg, strand='-')

        for i, junc in enumerate(graph.mJunctions):
            self._draw_edge(junc.mEdgeA, junc)
            self._draw_edge(junc.mEdgeB, junc)

    def _draw_seg(self, seg):
        if seg.mCov.mWeight == 0:
            color = 'gray'
        else:
            color = 'black'
        rect = svgobject.Rect(seg.insert, seg.size, stroke='none', fill='none').render()
        seg_repr = '{}{}\n{}'.format(seg.mType, seg.mIdx, int(seg.mCov.mWeight))
        label = svgobject.Text(seg_repr, font_size=10, fill=color,
                               anchor_object=rect, related_insert='middle:center').render()

        self.canvas.add(rect)
        self.canvas.add(label)
        return rect

    def _draw_seg_strand(self, seg, strand='+'):
        if seg.mCov.mWeight == 0:
            color = 'gray'
        else:
            color = 'black'

        if strand == '+':
            insert = seg.insert
            size = seg.size[0], 0
        else:
            insert = seg.insert[0], seg.insert[1] + seg.size[1]
            size = seg.size[0], 0
        strand_line = svgobject.Line(insert, size, stroke=color, fill='none').render()

        c1 = svgobject.Circle(insert, (self.vertex_width, self.vertex_width), stroke=color, fill='none',
                              self_pos='left:center').render()
        c2 = svgobject.Circle((insert[0]+size[0], insert[1]), (10, 10), stroke=color, fill='none',
                              self_pos='right:center').render()

        c1_label = svgobject.Text(strand, font_size=10, fill=color,
                                  anchor_object=c1, related_insert='middle:center').render()
        c2_label = svgobject.Text(strand, font_size=10, fill=color,
                                  anchor_object=c2, related_insert='middle:center').render()

        self.canvas.add(strand_line)
        self.canvas.add(c1)
        self.canvas.add(c2)
        self.canvas.add(c1_label)
        self.canvas.add(c2_label)

    def _draw_edge(self, edge, junc):

        if edge.mSourceV.mType == edge.mTargetV.mType:
            pass
        seg = edge.mSourceV.mSeg
        if edge.mSourceV.mDir == '+':
            source_xy = seg.insert[0] + seg.size[0] + self.vertex_width, seg.insert[1]
        else:
            source_xy = seg.insert[0] - self.vertex_width, seg.insert[1] + seg.size[1]

        seg = edge.mTargetV.mSeg
        if edge.mTargetV.mDir == '+':
            target_xy = seg.insert[0] - self.vertex_width, seg.insert[1]
        else:
            target_xy = seg.insert[0] + seg.size[0] + self.vertex_width, seg.insert[1] + seg.size[1]

        if junc.mCov.mImaginary:
            stroke_dasharray = "4"
        else:
            stroke_dasharray = ""

        if edge.mCov.mWeight:
            large_sv = True
        else:
            large_sv = False

        insert = source_xy
        size = target_xy[0] - source_xy[0], target_xy[1] - source_xy[1]
        sv_join, _, _ = complexsv.SVJoin(insert, size, cut_y=size[1],
                                         join_type=junc.join_type, large_sv=large_sv,
                                         stroke_dasharray=stroke_dasharray,
                                         enable_horizontal=True,
                                         ).render()

        cn = int(edge.mCov.mWeight)
        text = svgwrite.text.Text("")
        label = svgwrite.text.TextPath(sv_join, cn, startOffset='50%', text_anchor='middle',
                                       font_size=8, font_family='Arial')

        text.add(label)

        sv_join = self.canvas.add(sv_join)
        self.canvas.add(text)


def run_call(sample=None, region_fn=None, sv_fn=None, **args):
    pass


def run_draw(graphs, out_dir=None, **args):
    title = 'local haplotype graph'
    group_width = 250
    size = (1000, 10000)

    plot = GraphsPlot(graphs=graphs, out_dir=out_dir, title=title, size=size,
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
