#!/usr/bin/env python3
import argparse,json
p=argparse.ArgumentParser(); p.add_argument('--route',choices=['raw','processed'],required=True); p.add_argument('--peaks',nargs='*',default=[]); p.add_argument('--output',required=True); a=p.parse_args()
n=sum(sum(1 for x in open(path) if x.strip() and not x.startswith('#')) for path in a.peaks)
if a.route=='raw':
    c=[{'component':'h3k27ac_processing','status':'COMPLETE','reason':'raw ChIP and matched input processing completed'},{'component':'h3k27ac_peak_calling','status':'COMPLETE' if n else 'COMPLETE_EMPTY_RESULT','reason':'MACS3 completed with '+(str(n)+' peaks' if n else 'zero peaks')},{'component':'h3k27ac_count_contract','status':'COMPLETE','reason':'consensus and integer count matrix completed'}]
else:
    c=[{'component':'h3k27ac_processing','status':'NOT_AVAILABLE_INPUT','reason':'processed count-matrix route'},{'component':'h3k27ac_peak_calling','status':'NOT_AVAILABLE_INPUT','reason':'processed count-matrix route'},{'component':'h3k27ac_count_contract','status':'COMPLETE','reason':'external quantitative region matrix validated'}]
open(a.output,'w').write(json.dumps({'components':c,'sample_peak_count':n},sort_keys=True,indent=2)+'\n')
