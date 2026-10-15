#!/usr/bin/env python3
import argparse, csv, gzip, json, re, statistics

def number(pattern,text,default='NA'):
    match=re.search(pattern,text,re.I|re.M)
    return match.group(1).replace(',','') if match else default

def context_counts(path,control=None):
    methylated=unmethylated=0
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as handle:
        for line in handle:
            fields=line.rstrip('\n').split('\t')
            if len(fields)<5 or (control and fields[2]!=control): continue
            if fields[1]=='+': methylated+=1
            elif fields[1]=='-': unmethylated+=1
    return methylated,unmethylated

def main():
    p=argparse.ArgumentParser()
    for key in ('contract','rejected','mbias','split-report','chg','chh','dedup-report','alignment-report','sample-id','library-id','output','json-output'): p.add_argument('--'+key,required=True)
    p.add_argument('--trim-report', required=True, action='append')
    p.add_argument('--min-coverage',type=int,default=5); p.add_argument('--conversion-control-contig',default=''); a=p.parse_args()
    coverages=[]; meth=unmeth=0
    with open(a.contract,encoding='utf-8') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            cov=int(row['coverage']); coverages.append(cov); meth+=int(row['methylated_count']); unmeth+=int(row['unmethylated_count'])
    alignment=open(a.alignment_report,encoding='utf-8',errors='ignore').read(); dedup=open(a.dedup_report,encoding='utf-8',errors='ignore').read(); trim='\n'.join(open(path,encoding='utf-8',errors='ignore').read() for path in a.trim_report); split=open(a.split_report,encoding='utf-8',errors='ignore').read()
    chg_m,chg_u=context_counts(a.chg); chh_m,chh_u=context_counts(a.chh)
    if a.conversion_control_contig:
        cm1,cu1=context_counts(a.chg,a.conversion_control_contig); cm2,cu2=context_counts(a.chh,a.conversion_control_contig); total=cm1+cu1+cm2+cu2
        conversion_status='COMPLETE' if total else 'NOT_AVAILABLE_INPUT'; conversion=f'{(cu1+cu2)/total:.8f}' if total else 'NA'
    else: conversion_status='NOT_AVAILABLE_INPUT'; conversion='NA'
    total_cpg=meth+unmeth
    metrics={
      'sample_id':a.sample_id,'library_id':a.library_id,
      'input_pairs':number(r'(?:Total reads processed|Total number of read pairs analysed)[:\s]+([0-9,]+)',trim),
      'retained_pairs':number(r'Sequence pairs analysed in total[:\s]+([0-9,]+)',alignment),
      'unique_aligned_pairs':number(r'Number of paired-end alignments with a unique best hit[:\s]+([0-9,]+)',alignment),
      'mapping_efficiency_percent':number(r'Mapping efficiency[:\s]+([0-9.]+)%',alignment),
      'duplicates_removed':number(r'(?:Total number of alignments removed|Duplicated alignments removed)[:\s]+([0-9,]+)',dedup),
      'duplicate_fraction':'NA',
      'observed_cpg_loci':len(coverages),'nonreference_cpg_calls_excluded':max(0,sum(1 for _ in open(a.rejected,encoding='utf-8'))-1),'cpg_loci_min_coverage':sum(x>=a.min_coverage for x in coverages),
      'mean_cpg_coverage':f'{statistics.mean(coverages):.8f}' if coverages else 'NA','median_cpg_coverage':f'{statistics.median(coverages):.8f}' if coverages else 'NA',
      'global_cpg_methylation_fraction':f'{meth/total_cpg:.8f}' if total_cpg else 'NA',
      'chg_methylation_fraction':f'{chg_m/(chg_m+chg_u):.8f}' if chg_m+chg_u else 'NA','chh_methylation_fraction':f'{chh_m/(chh_m+chh_u):.8f}' if chh_m+chh_u else 'NA',
      'mbias_status':'COMPLETE','conversion_control_status':conversion_status,'conversion_efficiency_fraction':conversion
    }
    try:
        removed=int(metrics['duplicates_removed']); aligned=int(metrics['unique_aligned_pairs']); metrics['duplicate_fraction']=f'{removed/aligned:.8f}' if aligned else 'NA'
    except ValueError: pass
    with open(a.output,'w',newline='',encoding='utf-8') as out:
        writer=csv.DictWriter(out,fieldnames=list(metrics),delimiter='\t',lineterminator='\n'); writer.writeheader(); writer.writerow(metrics)
    with open(a.json_output,'w',encoding='utf-8') as out: json.dump(metrics,out,sort_keys=True,indent=2); out.write('\n')

if __name__=='__main__': main()
