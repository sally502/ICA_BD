#!/usr/bin/env python3
"""Retrieve QuickGO annotations from verified UniProt TSVs.

Usage:
  python3 02_quickgo_fixed.py uniprot_annotations_human.tsv
  python3 02_quickgo_fixed.py uniprot_annotations_mouse.tsv
  python3 02_quickgo_fixed.py input.tsv custom_output.tsv
  python3 02_quickgo_fixed.py input.tsv --overwrite
Requires: requests
"""
import argparse
import csv
import re
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE = 'https://www.ebi.ac.uk/QuickGO/services'
COLUMNS = ['Pair_ID','Species','Tax_ID','NCBI_Gene_ID','Gene_Symbol',
           'UniProt_Accession','GO_ID','GO_Name','GO_Aspect','Evidence_Code',
           'Qualifier','Reference','Assigned_By','Annotation_Taxon_ID','Status']
TAXA = {'Human':'9606', 'Mouse':'10090'}


def make_session():
    session = requests.Session()
    retry = Retry(total=4, backoff_factor=1, status_forcelist=[429,500,502,503,504],
                  respect_retry_after_header=True)
    session.mount('https://', HTTPAdapter(max_retries=retry))
    session.headers.update({'Accept':'application/json',
                            'User-Agent':'BioDB-ICA-student/1.0'})
    return session


def fetch_json(session, url, params=None):
    response = session.get(url, params=params, timeout=45)
    response.raise_for_status()
    return response.json()


def fetch_annotations(session, accession, taxid):
    """Fetch all pages; return (annotations, number of reported hits)."""
    results = []
    page = 1
    while True:
        data = fetch_json(session, BASE+'/annotation/search', {
            'geneProductId':'UniProtKB:'+accession,
            'taxonId':taxid, 'limit':100, 'page':page
        })
        batch = data.get('results', [])
        results.extend(batch)
        page_info = data.get('pageInfo') or {}
        total_pages = page_info.get('total')
        if total_pages is not None:
            if page >= int(total_pages):
                break
        elif not batch or len(batch) < 100:
            break
        if not batch:
            break
        page += 1
        if page > 1000:
            raise RuntimeError('Pagination safety limit exceeded')
        time.sleep(0.2)
    return results


def fetch_go_names(session, go_ids):
    names = {}
    for i in range(0, len(go_ids), 25):
        chunk = go_ids[i:i+25]
        try:
            data = fetch_json(session, BASE+'/ontology/go/terms/'+','.join(chunk))
            for term in data.get('results', []):
                if term.get('id'):
                    names[term['id']] = term.get('name','')
        except (requests.RequestException, ValueError) as exc:
            print('WARNING: GO name lookup failed:',str(exc)[:140])
        time.sleep(0.2)
    return names


def scalar(value):
    if value is None:
        return ''
    if isinstance(value, list):
        return '; '.join(str(x) for x in value)
    return str(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', nargs='?', type=Path)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    with args.input.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'Species','Tax_ID','NCBI_Gene_ID','Gene_Symbol',
                    'UniProt_Accession','Mapping_Status'}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError('Missing columns in UniProt TSV: '+', '.join(sorted(missing)))
        genes = list(reader)
    if not genes:
        raise ValueError('Input TSV has no rows')
    species_set = {g['Species'].strip().lower() for g in genes}
    suffix = next(iter(species_set)) if len(species_set)==1 else 'human_mouse'
    output = args.output or args.input.with_name('quickgo_annotations_'+suffix+'.tsv')
    if output.resolve() == args.input.resolve():
        raise ValueError('Output must not be the same as input')
    if output.exists() and not args.overwrite:
        raise FileExistsError(f'{output} already exists. Use another filename or --overwrite.')

    print('Input species:',sorted(species_set),'; output:',output,flush=True)
    session = make_session()
    rows = []
    go_ids = set()
    for gene in genes:
        species = gene['Species'].strip()
        tax = gene['Tax_ID'].strip()
        accession = gene['UniProt_Accession'].strip()
        status = gene['Mapping_Status'].strip()
        base = {key:gene.get(key,'') for key in
                ('Pair_ID','Species','Tax_ID','NCBI_Gene_ID','Gene_Symbol','UniProt_Accession')}
        if species not in TAXA or tax != TAXA[species]:
            rows.append({**base,'Status':'INVALID_SPECIES_OR_TAX_ID'})
            print(species, gene['Gene_Symbol'], 'INVALID_SPECIES_OR_TAX_ID',flush=True)
            continue
        if not accession or not status.startswith('VERIFIED_'):
            rows.append({**base,'Status':'SKIPPED_UNVERIFIED_UNIPROT'})
            print(species, gene['Gene_Symbol'], 'SKIPPED_UNVERIFIED_UNIPROT',flush=True)
            continue
        if not re.fullmatch(r'[A-Z0-9]{6,10}(?:-\d+)?', accession):
            rows.append({**base,'Status':'INVALID_ACCESSION'})
            continue
        try:
            annotations = fetch_annotations(session, accession, tax)
            count = 0
            for ann in annotations:
                # Query taxon should be verified against returned annotations.
                ann_tax = scalar(ann.get('taxonId'))
                if ann_tax and ann_tax != tax:
                    continue
                go_id = ann.get('goId','')
                if go_id:
                    go_ids.add(go_id)
                rows.append({**base,
                    'GO_ID':go_id,
                    'GO_Aspect':scalar(ann.get('goAspect')),
                    'Evidence_Code':scalar(ann.get('goEvidence') or ann.get('evidenceCode')),
                    'Qualifier':scalar(ann.get('qualifier')),
                    'Reference':scalar(ann.get('reference')),
                    'Assigned_By':scalar(ann.get('assignedBy')),
                    'Annotation_Taxon_ID':ann_tax,
                    'Status':'OK'})
                count += 1
            if not annotations:
                rows.append({**base,'Status':'NO_ANNOTATIONS'})
            elif not count:
                rows.append({**base,'Status':'TAXON_MISMATCH'})
            print(species, gene['Gene_Symbol'], count, 'annotations',flush=True)
        except (requests.RequestException, ValueError, RuntimeError) as exc:
            rows.append({**base,'Status':'API_ERROR: '+str(exc)[:160]})
            print(species, gene['Gene_Symbol'], 'API_ERROR:',str(exc)[:100],flush=True)
        time.sleep(0.2)

    print('Looking up names for',len(go_ids),'unique GO terms...',flush=True)
    names = fetch_go_names(session, sorted(go_ids))
    for row in rows:
        row['GO_Name'] = names.get(row.get('GO_ID',''),'')
    with output.open('x' if not args.overwrite else 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    print('Saved:',output,'(',len(rows),'rows )')


if __name__ == '__main__':
    main()
