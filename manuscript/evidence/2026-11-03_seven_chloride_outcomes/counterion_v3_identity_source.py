"""Stable publication identities for declared chloride starts and repair branches."""
import re


def publication_identity(start, branch):
    if not isinstance(start,str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',start):
        raise ValueError('Unknown declared chloride start')
    if branch not in (None,'minus','plus'):
        raise ValueError('Unknown signed repair branch')
    stem='N0' if start=='N0_above' and branch is not None else start
    key=stem+('_'+branch if branch else '')
    label=start.replace('_',' ')
    if branch:label+=r', $'+('-' if branch=='minus' else '+')+r'$ mode'
    return key,label


def validate_identities(records, starts):
    if not starts or len(set(starts))!=len(starts) or not records:
        raise ValueError('Distinct declared starts and nonempty checked outcomes required')
    for start in starts:publication_identity(start,None)
    identities=[]
    for row in records:
        if row['original_start'] not in starts:raise ValueError('Outcome is outside the declared starts')
        identities.append(publication_identity(row['original_start'],row['branch']))
    if len({key for key,label in identities})!=len(identities):
        raise ValueError('Duplicate outcome identity')
    return identities
