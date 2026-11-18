#!/usr/bin/env python3
"""Retrieve Reactome pathways for verified UniProt entries, output TSV.

Usage:
  python3 03_reactome_fixed.py uniprot_annotations_human.tsv
  python3 03_reactome_fixed.py uniprot_annotations_mouse.tsv
  python3 03_reactome_fixed.py input.tsv custom_output.tsv
  python3 03_reactome_fixed.py input.tsv --overwrite
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

BASE = 'https://reactome.org/ContentService'
SPECIES = {'Human': ('9606', 'Homo sapiens'), 'Mouse': ('10090', 'Mus musculus')}
COLUMNS = ['Pair_ID', 'Species', 'Tax_ID', 'NCBI_Gene_ID', 'Gene_Symbol',
           'UniProt_Accession', 'Pathway_ID', 'Pathway_Name', 'Reactome_Species',
           'Is_Inferred', 'Reactome_URL', 'Status']


def make_session():
    session = requests.Session()
    retry = Retry(total=4, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504],
                  respect_retry_after_header=True)
    session.mount('https://', HTTPAdapter(max_retries=retry))
    session.headers.update({'Accept': 'application/json',
                            'User-Agent': 'BioDB-ICA-student/1.0'})
    return session


def fetch_pathways(session, accession):
    url = f'{BASE}/data/mapping/UniProt/{accession}/pathways'
    response = session.get(url, timeout=45)
    if response.status_code == 404:
        return []
    response.raise_for_status()
    data = response.json()
    if not isinstance(data, list):
        raise ValueError('Expected a list of Reactome pathways')
    return data


def species_name(pathway):
    """Reactome pathway JSON generally uses speciesName; support species object too."""
    name = pathway.get('speciesName')
    if isinstance(name, str) and name:
        return name
    species = pathway.get('species')
    if isinstance(species, dict):
        return species.get('displayName') or species.get('name') or ''
    return species if isinstance(species, str) else ''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', nargs='?', type=Path)
    parser.add_argument('--overwrite', action='store_true', help='Allow replacing an existing output file')
    args = parser.parse_args()

    with args.input.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'Species', 'Tax_ID', 'NCBI_Gene_ID', 'Gene_Symbol',
                    'UniProt_Accession', 'Mapping_Status'}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError('Missing columns in UniProt TSV: ' + ', '.join(sorted(missing)))
        genes = list(reader)
    if not genes:
        raise ValueError('Input TSV has no rows')

    species_set = {g['Species'].strip().lower() for g in genes}
    suffix = next(iter(species_set)) if len(species_set) == 1 else 'human_mouse'
    output = args.output or args.input.with_name('reactome_pathways_' + suffix + '.tsv')
    if output.resolve() == args.input.resolve():
        raise ValueError('Output file must differ from input file')
    if output.exists() and not args.overwrite:
        raise FileExistsError(f'{output} already exists; choose a new filename or use --overwrite')

    print('Input species:', sorted(species_set), '; output:', output, flush=True)
    session = make_session()
    rows = []
    for gene in genes:
        base = {key: gene.get(key, '') for key in
                ('Pair_ID', 'Species', 'Tax_ID', 'NCBI_Gene_ID', 'Gene_Symbol', 'UniProt_Accession')}
        species = gene['Species'].strip()
        taxid = gene['Tax_ID'].strip()
        accession = gene['UniProt_Accession'].strip()
        mapping = gene['Mapping_Status'].strip()
        label = f"{species} {gene['Gene_Symbol']}"

        if species not in SPECIES or taxid != SPECIES[species][0]:
            rows.append({**base, 'Status': 'INVALID_SPECIES_OR_TAX_ID'})
            print(label, 'INVALID_SPECIES_OR_TAX_ID', flush=True)
            continue
        if not accession or not mapping.startswith('VERIFIED_'):
            rows.append({**base, 'Status': 'SKIPPED_UNVERIFIED_UNIPROT'})
            print(label, 'SKIPPED_UNVERIFIED_UNIPROT', flush=True)
            continue
        if not re.fullmatch(r'[A-Z0-9]{6,10}(?:-\d+)?', accession):
            rows.append({**base, 'Status': 'INVALID_ACCESSION'})
            print(label, 'INVALID_ACCESSION', flush=True)
            continue

        try:
            pathways = fetch_pathways(session, accession)
            expected_species = SPECIES[species][1]
            seen = set()
            count = 0
            unknown_species = 0
            other_species = 0
            for pathway in pathways:
                if not isinstance(pathway, dict):
                    continue
                returned_species = species_name(pathway)
                # Do not silently mix species or treat an unknown species as verified.
                if not returned_species:
                    unknown_species += 1
                    continue
                if returned_species != expected_species:
                    other_species += 1
                    continue
                pathway_id = pathway.get('stId') or ''
                if not pathway_id or pathway_id in seen:
                    continue
                seen.add(pathway_id)
                inferred = pathway.get('isInferred')
                rows.append({**base,
                             'Pathway_ID': pathway_id,
                             'Pathway_Name': pathway.get('displayName') or pathway.get('name') or '',
                             'Reactome_Species': returned_species,
                             'Is_Inferred': '' if inferred is None else str(inferred),
                             'Reactome_URL': 'https://reactome.org/content/detail/' + pathway_id,
                             'Status': 'OK'})
                count += 1
            if count == 0:
                status = ('NO_PATHWAYS' if not pathways else
                          'SPECIES_NOT_VERIFIABLE' if unknown_species and not other_species else
                          'NO_PATHWAYS_FOR_SPECIES')
                rows.append({**base, 'Status': status})
                print(label, status, flush=True)
            else:
                print(label, count, 'pathways', flush=True)
            if unknown_species:
                print('  WARNING:', unknown_species, 'pathways missing species information', flush=True)
        except (requests.RequestException, ValueError, RuntimeError) as exc:
            rows.append({**base, 'Status': 'API_ERROR: ' + str(exc)[:160]})
            print(label, 'API_ERROR:', str(exc)[:140], flush=True)
        time.sleep(0.25)

    # Exclusive creation by default protects previous human/mouse results.
    with output.open('w' if args.overwrite else 'x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    print('Saved:', output, '(', len(rows), 'rows )', flush=True)


if __name__ == '__main__':
    main()
