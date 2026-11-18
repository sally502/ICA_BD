#!/usr/bin/env python3
"""Fetch verified UniProtKB annotations from NCBI JSON or human/mouse TSV.
Usage: python 01_uniprot.py input.json [output.tsv]
       python 01_uniprot.py orthologs.tsv [output.tsv]
Requires: requests
"""
import argparse
import csv
import json
import re
import time
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

API = 'https://rest.uniprot.org/uniprotkb/search'
HEADERS = ['Pair_ID','Species','Tax_ID','NCBI_Gene_ID','Gene_Symbol','Ensembl_Gene_ID',
           'UniProt_Accession','UniProt_Entry_Name','Protein_Name','Reviewed',
           'Protein_Length_aa','Molecular_Weight_Da','Protein_Function',
           'Subcellular_Location','Domains','Sequence','Mapping_Status','UniProt_URL']


def session_with_retries():
    s = requests.Session()
    retry = Retry(total=4, backoff_factor=1, status_forcelist=[429,500,502,503,504],
                  respect_retry_after_header=True)
    s.mount('https://', HTTPAdapter(max_retries=retry))
    s.headers.update({'User-Agent':'BioDB-ICA-student/1.0', 'Accept':'application/json'})
    return s


def load_genes(path):
    if path.suffix.lower() == '.json':
        with path.open(encoding='utf-8') as f:
            data = json.load(f)
        if isinstance(data, dict):
            data = data.get('genes', data.get('results', []))
        for idx, record in enumerate(data, 1):
            summary = record.get('summary') or {}
            if not isinstance(summary, dict):
                summary = {}
            ensembl = record.get('ensembl_gene_ids') or []
            yield dict(pair=str(idx), species='Human', tax='9606',
                       gene_id=str(record.get('gene_id') or summary.get('uid') or ''),
                       symbol=str(record.get('query_symbol') or summary.get('name') or ''),
                       ensembl=str(ensembl[0]) if ensembl else '')
    else:
        with path.open(encoding='utf-8-sig', newline='') as f:
            for idx, row in enumerate(csv.DictReader(f, delimiter='\t'), 1):
                for species, tax in [('Human','9606'),('Mouse','10090')]:
                    prefix = species + '_'
                    gene_id = (row.get(prefix+'Gene_ID') or '').strip()
                    symbol = (row.get(prefix+'Symbol') or '').strip()
                    ensembl = (row.get(prefix+'Ensembl_Gene_ID') or '').strip()
                    if gene_id or symbol or ensembl:
                        yield dict(pair=str(idx), species=species, tax=tax,
                                   gene_id=gene_id, symbol=symbol, ensembl=ensembl)


def query_uniprot(s, gene):
    # Search by symbol and organism; never trust a symbol match without xref verification.
    symbol, tax = gene['symbol'], gene['tax']
    if not symbol:
        return [], 'MISSING_SYMBOL'
    query = f'(gene_exact:{symbol}) AND (organism_id:{tax})'
    params = {'query':query, 'format':'json', 'size':100}
    response = s.get(API, params=params, timeout=35)
    response.raise_for_status()
    return response.json().get('results', []), None


def crossrefs(entry, db):
    return [x.get('id','') for x in entry.get('uniProtKBCrossReferences', [])
            if x.get('database') == db]


def verified(entry, gene):
    gene_id, ens = gene['gene_id'], gene['ensembl']
    if gene_id and gene_id in crossrefs(entry, 'GeneID'):
        return 'VERIFIED_GENE_ID'
    if ens:
        for ref in crossrefs(entry, 'Ensembl'):
            if ref.split('.')[0] == ens.split('.')[0]:
                return 'VERIFIED_ENSEMBL'
    return None


def choose_entry(entries, gene):
    matches = [(e, verified(e, gene)) for e in entries]
    matches = [(e, status) for e, status in matches if status]
    if not matches:
        return None, 'NO_VERIFIED_MATCH' if entries else 'NO_RESULTS'
    matches.sort(key=lambda item: (item[0].get('entryType') == 'UniProtKB reviewed (Swiss-Prot)',
                                   item[1] == 'VERIFIED_GENE_ID',
                                   item[0].get('primaryAccession','')), reverse=True)
    reviewed = [x for x in matches if x[0].get('entryType') == 'UniProtKB reviewed (Swiss-Prot)']
    if len(reviewed) > 1:
        return None, 'MULTIPLE_REVIEWED_MATCHES'
    if not reviewed and len(matches) > 1:
        return None, 'MULTIPLE_UNREVIEWED_MATCHES'
    return matches[0]


def comments(entry, comment_type):
    result = []
    for c in entry.get('comments', []):
        if c.get('commentType') != comment_type:
            continue
        if comment_type == 'SUBCELLULAR LOCATION':
            for loc in c.get('subcellularLocations', []):
                name = (loc.get('location') or {}).get('value')
                if name:
                    result.append(name)
        else:
            result += [t.get('value','') for t in c.get('texts', []) if t.get('value')]
    return '; '.join(dict.fromkeys(result))


def protein_name(entry):
    pd = entry.get('proteinDescription') or {}
    for key in ('recommendedName','submissionNames'):
        v = pd.get(key)
        if isinstance(v, list):
            v = v[0] if v else {}
        name = (v or {}).get('fullName') or {}
        if name.get('value'):
            return name['value']
    return ''


def flatten(entry, gene, status):
    seq = entry.get('sequence') or {}
    domains = []
    for feature in entry.get('features', []):
        if feature.get('type') in ('Domain','Region','Repeat','Motif'):
            label = feature.get('description') or feature.get('type')
            loc = feature.get('location') or {}
            start = (loc.get('start') or {}).get('value','')
            end = (loc.get('end') or {}).get('value','')
            domains.append(f'{label} ({start}-{end})' if start and end else label)
    return {'UniProt_Accession':entry.get('primaryAccession',''),
            'UniProt_Entry_Name':entry.get('uniProtkbId',''),
            'Protein_Name':protein_name(entry),
            'Reviewed':str(entry.get('entryType') == 'UniProtKB reviewed (Swiss-Prot)'),
            'Protein_Length_aa':seq.get('length',''),
            'Molecular_Weight_Da':seq.get('molWeight',''),
            'Protein_Function':comments(entry,'FUNCTION'),
            'Subcellular_Location':comments(entry,'SUBCELLULAR LOCATION'),
            'Domains':'; '.join(dict.fromkeys(domains)),
            'Sequence':seq.get('value',''),
            'Mapping_Status':status,
            'UniProt_URL':'https://www.uniprot.org/uniprotkb/'+entry.get('primaryAccession','')+'/entry'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path)
    p.add_argument('output', nargs='?', type=Path, default=Path('uniprot_annotations.tsv'))
    args = p.parse_args()
    genes = list(load_genes(args.input))
    s = session_with_retries()
    with args.output.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=HEADERS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for gene in genes:
            row = {'Pair_ID':gene['pair'],'Species':gene['species'],'Tax_ID':gene['tax'],
                   'NCBI_Gene_ID':gene['gene_id'],'Gene_Symbol':gene['symbol'],
                   'Ensembl_Gene_ID':gene['ensembl']}
            try:
                entries, error = query_uniprot(s, gene)
                entry, status = (None, error) if error else choose_entry(entries, gene)
                if entry:
                    row.update(flatten(entry, gene, status))
                else:
                    row['Mapping_Status'] = status
            except (requests.RequestException, ValueError) as exc:
                row['Mapping_Status'] = 'API_ERROR: '+str(exc)[:180]
            writer.writerow(row)
            print(f"{gene['species']:5} {gene['symbol']:10} {row['Mapping_Status']} {row.get('UniProt_Accession','')}")
            time.sleep(0.2)
    print('Saved:', args.output)

if __name__ == '__main__':
    main()
