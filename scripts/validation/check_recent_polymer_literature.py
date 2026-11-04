#!/usr/bin/env python3
"""Check the two recent polymer references against saved primary DOI metadata."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def normalize(value):
    return ' '.join(value.replace('{','').replace('}','').replace('‐','-').replace('–','-').replace('--','-').split()).casefold()


def closing_brace(text,start):
    depth=0;escaped=False
    for index in range(start,len(text)):
        char=text[index]
        if escaped:escaped=False;continue
        if char=='\\':escaped=True;continue
        if char=='{':depth+=1
        if char=='}':
            depth-=1
            if depth==0:return index
    raise ValueError('Unbalanced bibliographic braces')


def article_fields(text,key):
    starts=list(re.finditer(r'@article\s*\{\s*'+re.escape(key)+r'\s*,',text,re.I))
    if len(starts)!=1:
        raise ValueError('Exactly one cited bibliographic entry required')
    start=starts[0];end=closing_brace(text,text.index('{',start.start()))
    body=text[start.end():end];fields={}
    for match in re.finditer(r'(?m)^\s*(\w+)\s*=\s*',body):
        position=match.end()
        if position>=len(body) or body[position]!='{':
            raise ValueError('Explicit brace-delimited bibliographic fields required')
        close=closing_brace(body,position)
        name=match[1].lower()
        if name in fields:raise ValueError('Duplicated bibliographic field')
        fields[name]=body[position+1:close]
    return fields


def audit(bibliography, sources):
    registry = json.loads(sources.read_text())['records']
    text = bibliography.read_text()
    rows, used = [], set()
    for record in registry:
        path = sources.parent/record['metadata_file']
        if hashlib.sha256(path.read_bytes()).hexdigest() != record['metadata_sha256']:
            raise ValueError('Saved primary metadata changed')
        data = json.loads(path.read_text())['message']
        if normalize(data['DOI']) != normalize(record['doi']) or normalize(record['doi']) in used:
            raise ValueError('DOI identity is different or duplicated')
        used.add(normalize(record['doi']))
        fields = article_fields(text,record['key'])
        expected = dict(doi=data['DOI'], title=data['title'][0], journal=data['container-title'][0],
            author=' and '.join(a['family']+', '+a['given'] for a in data['author']),
            year=str(data['issued']['date-parts'][0][0]),volume=data['volume'],number=data['issue'],
            pages=data.get('page') or data['article-number'])
        if any(normalize(fields.get(key,'')) != normalize(value) for key,value in expected.items()):
            raise ValueError('Bibliographic fields differ from primary DOI metadata: '+record['key'])
        rows.append(dict(key=record['key'],doi=record['doi'],fields_checked=sorted(expected),metadata_sha256=record['metadata_sha256']))
    return dict(status='recent_polymer_references_match_primary_metadata',references=len(rows),records=rows,
                route_or_property_promotion_inferred=False,
                bibliography_sha256=hashlib.sha256(bibliography.read_bytes()).hexdigest())


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bibliography',type=Path,required=True)
    parser.add_argument('--sources',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();result=audit(args.bibliography,args.sources)
    with args.out.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result))
