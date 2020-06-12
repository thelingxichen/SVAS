#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    tenxtools.fusesv
    ~~~~~~~~~~~~~~~~

    FuseSV output reader

    @Copyright: (c) 2018-10 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""


class SEG():

    def __init__(self, line):
        args = dict(x.split('=') for x in line.split(';'))
        self.set(**args)

    def set(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)


class Graph():

    def __init__(self, in_fn=None):
        self.seg_list = []
        self.read_fn(in_fn)

    def read_fn(self, in_fn):
        for line in open(fn, 'r'):
            if line.startswith('#'):
                continue
            splits = line.strip().split('\t')
            if 'SAMPLE' in line:
                self.sample = splits[1]
            if 'PLOIDY' in line:
                self.ploidy = splits[1]
            if 'SOURCE' in line:
                self.ploidy = splits[1]
            if 'SINK' in line:
                self.ploidy = splits[1]
            if line.startswith('SEG'):
                self.seg_list.append(SEG(splits[1]))


if __name__ == "__main__":

    fn = '/home/chenlingxi/mnt/chenlingxi/workspace/FS_Projects/WHOC.WGS_10X.batch01/ComplexSV/svaba/OC001T/group/group1/cn/out/OC001T_chr7_13499000_15001000.InputForLocalHapAlgorithm.txt'
    g = Graph(fn)
    for seg in g.seg_list:
        print(seg.CN_ORIG, seg.INTERVAL)
