#! /usr/bin/env python

import sys
from collections import Counter

import pysam
from telescope.utils import alignment
from telescope.utils.model import Assigner
from telescope.utils.model import process_overlap_frag
from telescope.utils.annotation import get_annotation_class

""" Constants """
NFKEY = "__no_feature" # "no feature" key
OMODE = "threshold"    # overlap mode
OTHRESH = 0.2          # overlap threshold     

COLORS = {
    'unique': '217,95,2',     # vermillion
    'best':   '230,171,2',    # yellow
    'lowconf': '209,236,228', # green    
    'white': '248,248,248',    
    # 'best':   '27,158,119',   # teal

    # 'green1': '164,216,201',
    # 'green2': '209,236,228',
    # 'green3': '232,245,241',
    # 'gray': '192,192,192',

}

def tag_multi(samin, samout, gtf, nfkey, omode, othresh, make_reports):
    # Parse annotation
    Annotation = get_annotation_class('intervaltree')
    annot = Annotation(gtf, 'locus')
    
    # Setup assigner
    assign = Assigner(annot, nfkey, omode, othresh).assign_func()
    
    # Setup input and output alignment files
    afin  = pysam.AlignmentFile(samin, "r")
    afout = pysam.AlignmentFile(samout, "w", template=afin)
    
    # Summary variables
    numalns_locus = Counter()
    report = []
    
    for ci, alns in alignment.fetch_fragments_seq(afin, until_eof=True):
        _mapped = [a for a in alns if not a.is_unmapped]
        _scores = [a.alnscore for a in _mapped]
        _feats = list(map(assign, _mapped))        
        _maps = process_overlap_frag(_mapped, _feats)
        
        _bestscore = max(_scores)
        _ambig = len(_mapped) > 1  
        _maps = process_overlap_frag(_mapped, _feats)
        for m in _maps:
            numalns_locus[m[1]] += 1
        
        read_id = alns[0].r1.query_name        
        num_best = sum(s == _bestscore for s in _scores)
        num_alns = len(_mapped)
        num_feats = len(set(f for f in _feats if f != nfkey))
        is_outside = any(f == nfkey for f in _feats)

        # Write report line
        report.append([read_id, _ambig, num_best, num_alns, num_feats, is_outside])
        
        if not _ambig:
            _mapped[0].set_tag('YC', COLORS['unique'])
        else:
            for a in _mapped:
                if a.r1.get_tag('ZT') == 'SEC':
                    a.set_tag('YC', COLORS['white'])
                else:
                    if a.alnscore == _bestscore:
                        a.set_tag('YC', COLORS['best'])
                    else:
                        a.set_tag('YC', COLORS['lowconf'])
                        # ratio = a.alnscore / _bestscore
                        # if ratio > 0.99:
                        #     a.set_tag('YC', COLORS['green0'])
                        # elif ratio > 0.975:
                        #     a.set_tag('YC', COLORS['green1'])
                        # elif ratio > 0.95:
                        #     a.set_tag('YC', COLORS['green2'])
                        # else:
                        #     a.set_tag('YC', COLORS['green3'])
        # Write all alignments
        for a in alns:
            _ = a.write(afout)
    
    afin.close()
    afout.close()
    
    if make_reports:
        with open('per_read.txt', 'w') as outh:
            print('id\tambig\tnum_best\tnum_alns\tnum_feats\toutside', file=outh)
            for l in report:
                print('\t'.join(map(str, l)), file=outh)
        
        with open('per_locus.txt', 'w') as outh:
            for t in numalns_locus.most_common():
                print('%s\t%d' % t, file=outh)
    

def console():
    import argparse
    parser = argparse.ArgumentParser(
        description='Tag reads in multi-alignment',
    )
    parser.add_argument('--overlap_mode', 
                        choices=['threshold', 'intersection-strict', 'union',],
                        default='threshold',
                        help='''Overlap mode. The method used to determine whether a 
                                fragment overlaps feature.'''
    )
    parser.add_argument('--overlap_threshold', 
                        default=0.2,
                        type=float,
                        help='''Fraction of fragment that must be contained within a 
                                feature to be assigned to that locus. Ignored if 
                                --overlap_method is not "threshold".'''
    )
    parser.add_argument('--no_feature_key',
                        default='__no_feature',
                        help='''Used internally to represent alignments. Must be different
                                from all other feature names.'''
    )
    parser.add_argument('--reports', action='store_true',
                        help='''Create per-read and per-locus reports (per_read.txt and 
                                per_locus.txt). Reports are output in working 
                                directory.'''
    )
    parser.add_argument('gtf')
    parser.add_argument('samin', nargs='?', default='-')    
    parser.add_argument('samout', nargs='?', default='-')
    args = parser.parse_args()
    tag_multi(
        args.samin, 
        args.samout, 
        args.gtf, 
        args.no_feature_key,
        args.overlap_mode,
        args.overlap_threshold,
        args.reports,
    )

if __name__ == '__main__':
    console()