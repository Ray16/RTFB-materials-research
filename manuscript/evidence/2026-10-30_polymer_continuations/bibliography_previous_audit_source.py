#!/usr/bin/env python
"""Audit BibTeX DOI identities against registered metadata without editing citations."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import difflib
import hashlib
import html
import json
from pathlib import Path
import re
import time
import unicodedata
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import Request, urlopen

PRODUCER_BYTES=Path(__file__).read_bytes()


def field(block, name):
    match=re.search(r'\b'+re.escape(name)+r'\s*=\s*([{"])',block,re.I)
    if not match:return None
    start=match.end();delimiter=match[1];depth=1
    for i in range(start,len(block)):
        escaped=i>0 and block[i-1]=='\\'
        if delimiter=='"' and block[i]=='"' and not escaped:return block[start:i]
        if delimiter=='{' and not escaped:
            if block[i]=='{':depth+=1
            if block[i]=='}':depth-=1
            if depth==0:return block[start:i]
    raise ValueError('Unclosed BibTeX field '+name)


def entries(text):
    rows=[]
    for block in re.split(r'(?=@\w+\{)',text):
        if not block.strip():continue
        match=re.match(r'@(\w+)\{([^,]+),',block)
        if not match:raise ValueError('Unsupported BibTeX record')
        rows.append(dict(key=match[2],kind=match[1],**{name:field(block,name) for name in
                        ('doi','title','author','year','journal','volume','pages','url','howpublished')}))
    if len({r['key'] for r in rows})!=len(rows):raise ValueError('Duplicate citation keys')
    return rows


def normalized(value):
    value=html.unescape(re.sub(r'<[^>]+>','',value or ''))
    for letter,name in [('ω','omega'),('Ω','omega'),('α','alpha'),('β','beta'),('γ','gamma')]:value=value.replace(letter,name)
    value=re.sub(r'\\(?:text|math)(?:rm|it|bf)?\s*','',value)
    value=re.sub(r'\\(?:v|H|c|k|r|b|u|t|d)(?![A-Za-z])\s*','',value)
    value=value.replace('\\omega','omega').replace('\\alpha','alpha').replace('\\beta','beta')
    value=unicodedata.normalize('NFKD',value).encode('ascii','ignore').decode()
    return ''.join(re.findall(r'[a-z0-9]+',value.lower()))


def get_metadata(doi):
    datacite=doi.lower().startswith('10.48550/')
    url=('https://api.datacite.org/dois/' if datacite else 'https://api.crossref.org/works/')+quote(doi,safe='')
    request=Request(url,headers={'User-Agent':'RedoxManuscriptBibliographyAudit/1.0 (read-only DOI metadata verification)'})
    for attempt in range(3):
        try:
            with urlopen(request,timeout=25) as response:raw=response.read()
            break
        except HTTPError as error:
            if error.code != 429 or attempt == 2:raise
            delay=float(error.headers.get('Retry-After',str(3*(attempt+1))))
            if not 0 <= delay <= 30:raise
            time.sleep(delay)
    parsed=json.loads(raw)
    if datacite:
        data=parsed['data']['attributes']
        metadata=dict(doi=data['doi'],title=data['titles'][0]['title'],year_candidates=[data['publicationYear']],
                      first_author=data['creators'][0].get('familyName',data['creators'][0].get('name','')),
                      authors=data['creators'])
    else:
        data=parsed['message'];years=[]
        for key in ('published','published-print','published-online','issued'):
            parts=data.get(key,{}).get('date-parts',[])
            if parts and parts[0]:years.append(parts[0][0])
        metadata=dict(doi=data['DOI'],title=data['title'][0],year_candidates=sorted(set(years)),
                      first_author=data.get('author',[{}])[0].get('family',''),authors=data.get('author',[]),
                      journal=data.get('container-title',[''])[0],volume=data.get('volume'),pages=data.get('page'),
                      article_number=data.get('article-number'))
    return dict(url=url,response_sha256=hashlib.sha256(raw).hexdigest(),registered=metadata,raw=parsed)


def compare(row,record):
    metadata=record['registered'];score=difflib.SequenceMatcher(None,normalized(row['title']),normalized(metadata['title'])).ratio()
    cited_family=(row['author'] or '').split(' and ')[0].split(',')[0]
    family_ok=normalized(cited_family)==normalized(metadata['first_author'])
    # DataCite sometimes supplies only a full personal name.
    if ',' not in metadata['first_author'] and normalized(cited_family) in normalized(metadata['first_author']):family_ok=True
    doi_ok=row['doi'].lower()==metadata['doi'].lower()
    year_ok=bool(row['year']) and int(row['year']) in metadata['year_candidates']
    return dict(status='metadata_matches' if doi_ok and score>=.98 and family_ok and year_ok else 'manual_review_required',
                DOI_matches=doi_ok,title_similarity=score,first_author_matches=family_ok,year_matches=year_ok)


def audit(bibliography,out,reuse=None):
    if out.exists():raise FileExistsError('Preserve previous bibliographic audits')
    source_bytes=bibliography.read_bytes();rows=entries(source_bytes.decode());out.mkdir(parents=True);(out/'metadata').mkdir()
    (out/'producer.py').write_bytes(PRODUCER_BYTES)
    cached={}
    if reuse is not None:
        prior=json.loads(reuse.read_text())
        for row in prior['checks']:
            if row.get('status') != 'metadata_matches':continue
            path=Path(row['metadata_file'])
            if hashlib.sha256(path.read_bytes()).hexdigest()!=row['metadata_sha256']:
                raise ValueError('Prior metadata cache changed')
            cached[row['citation']['doi'].lower()]=path
    def check(row):
        if not row['doi']:return dict(key=row['key'],status='non_DOI_source_requires_direct_review',citation=row)
        try:
            cached_path=cached.get(row['doi'].lower())
            if cached_path is None:
                time.sleep(.5)
                record=get_metadata(row['doi'])
            else:record=json.loads(cached_path.read_text())
            path=out/'metadata'/(row['key']+'.json')
            path.write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n')
            return dict(key=row['key'],citation=row,**compare(row,record),registered=record['registered'],
                        metadata_file=str(path),metadata_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        cached_from=str(cached_path) if cached_path is not None else None)
        except Exception as error:return dict(key=row['key'],status='metadata_unverified',citation=row,error=repr(error))
    with ThreadPoolExecutor(max_workers=1) as executor:checks=list(executor.map(check,rows))
    if bibliography.read_bytes()!=source_bytes:raise ValueError('Bibliography changed during audit')
    result=dict(status='bibliographic_identity_audit_not_content_validation',checked=time.time(),
                bibliography=str(bibliography),bibliography_sha256=hashlib.sha256(source_bytes).hexdigest(),
                producer_sha256=hashlib.sha256(PRODUCER_BYTES).hexdigest(),checks=checks,
                scope='DOI/title/first-author/year metadata only. Content and non-DOI sources require direct primary-source review.')
    (out/'summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--bib',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--reuse',type=Path)
    args=parser.parse_args();result=audit(args.bib.resolve(),args.out.resolve(),args.reuse.resolve() if args.reuse else None)
    from collections import Counter
    print(json.dumps(dict(counts=Counter(r['status'] for r in result['checks']))))
