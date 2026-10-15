#!/usr/bin/env python3
import argparse, csv, json

def main():
    p=argparse.ArgumentParser(); p.add_argument('--route',required=True,choices=['raw','processed']); p.add_argument('--validation-report',required=True); p.add_argument('--contracts',nargs='*',default=[]); p.add_argument('--output',required=True); a=p.parse_args()
    validation=json.load(open(a.validation_report)); loci=0
    for path in a.contracts:
        with open(path,encoding='utf-8') as handle: loci+=sum(1 for _ in csv.DictReader(handle,delimiter='\t'))
    if a.route=='raw':
        components=[{'component':'wgbs_processing','status':'COMPLETE','reason':'paired-end conventional WGBS processing completed'},{'component':'wgbs_methylation_extraction','status':'COMPLETE' if loci else 'COMPLETE_EMPTY_RESULT','reason':'CpG extraction completed' if loci else 'extraction completed with zero covered CpGs'},{'component':'wgbs_count_contract','status':'COMPLETE','reason':'count-preserving canonical CpG contract completed'}]
    else:
        quantitative=validation.get('quantitative_counts_available',False)
        components=[{'component':'wgbs_processing','status':'NOT_AVAILABLE_INPUT','reason':'externally processed WGBS route'},{'component':'wgbs_methylation_extraction','status':'NOT_AVAILABLE_INPUT','reason':'externally processed WGBS route'},{'component':'wgbs_count_contract','status':'COMPLETE' if quantitative else 'NOT_AVAILABLE_INPUT','reason':'external CpG counts validated' if quantitative else 'percentage-only descriptive input has no count contract'}]
    with open(a.output,'w',encoding='utf-8') as handle: json.dump({'components':components,'canonical_cpg_rows':loci},handle,sort_keys=True,indent=2); handle.write('\n')
if __name__=='__main__': main()
