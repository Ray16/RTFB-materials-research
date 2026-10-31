"""Standard-library replay of electronic implementation differences only."""
import math


def reconstruct(states,endpoints,hartree_eV):
    if not states or not endpoints or not math.isfinite(hartree_eV) or hartree_eV<=0:
        raise ValueError('Nonempty declared controls and a finite energy conversion required')
    aliases={};keys=set();differences={}
    for row in states:
        if row['key'] in keys or type(row['charge']) is not int:raise ValueError('Distinct exact state identities required')
        keys.add(row['key'])
        values=[row[name] for name in ('stored_Ha','installed_Ha','corrected_Ha')]
        if any(type(v) not in (int,float) or not math.isfinite(v) for v in values):raise ValueError('Finite scalar electronic energies required')
        old,installed,corrected=values
        differences[row['key']]=dict(installed_minus_stored_Ha=installed-old,
            corrected_minus_stored_Ha=corrected-old,corrected_minus_installed_Ha=corrected-installed)
        for alias in row['aliases']:
            if alias in aliases:raise ValueError('Ambiguous state alias')
            aliases[alias]=row
    if not {'benchmark:ferrocene:ox','benchmark:ferrocene:neu'}<=set(aliases):raise ValueError('Exact matched reference pair missing')
    fo=aliases['benchmark:ferrocene:ox'];fr=aliases['benchmark:ferrocene:neu']
    if fo['charge']!=fr['charge']+1:raise ValueError('One-electron Fc reference required')
    rows=[];endpoint_keys=set()
    for endpoint in endpoints:
        key=(endpoint['id'],endpoint['event'],endpoint.get('doi',''))
        if key in endpoint_keys:raise ValueError('Duplicate endpoint')
        endpoint_keys.add(key)
        names=endpoint['event'].split('->')
        if len(names)!=2:raise ValueError('Explicit state pair required')
        names=['benchmark:'+endpoint['id']+':'+state for state in names]
        if not set(names)<=set(aliases):raise ValueError('Endpoint state missing')
        a,b=[aliases[n] for n in names]
        if a['charge']!=b['charge']+1:raise ValueError('One-electron reduction orientation required')
        shifts={name.replace('_Ha','_V'):hartree_eV*(differences[a['key']][name]-differences[b['key']][name]
            -differences[fo['key']][name]+differences[fr['key']][name]) for name in differences[a['key']]}
        entry=dict(id=key[0],event=key[1],electronic_contribution_to_potential_shift_V=shifts)
        if 'doi' in endpoint:entry['doi']=endpoint['doi']
        rows.append(entry)
    return dict(states=len(states),endpoints=len(rows),state_differences_Ha=differences,endpoint_shifts=rows,
        maximum_state_energy_differences_Ha={name:max(abs(r[name]) for r in differences.values()) for name in next(iter(differences.values()))},
        maximum_electronic_potential_shifts_V={name:max(abs(r['electronic_contribution_to_potential_shift_V'][name]) for r in rows) for name in rows[0]['electronic_contribution_to_potential_shift_V']},
        nuclear_terms_recomputed=False,new_redox_potentials_qualified=False,native_G_qualified=False)


def molecular_derivatives(cases):
    """Reconstruct every saved Cartesian column without NumPy or quantum code."""
    rows=[];identities=set()
    flatten=lambda value:[x for atom in value for x in atom]
    norm=lambda value:math.sqrt(sum(x*x for x in value))
    for case in cases:
        identity=(case['device'],case['control'])
        if identity in identities:raise ValueError('Duplicate device/control')
        identities.add(identity);center=case['coordinates_Bohr'];dimension=3*len(center)
        if any(len(atom)!=3 for atom in center):raise ValueError('Three Cartesian coordinates per atom required')
        for variant in ('installed','corrected'):
            data=case['variants'][variant];matrix=data['Hessian'];points=data['points']
            if len(matrix)!=dimension or any(len(row)!=dimension for row in matrix):raise ValueError('Full Cartesian Hessian required')
            if len(points)!=1+4*dimension:raise ValueError('Complete declared two-step Cartesian point set required')
            if any(not math.isfinite(x) for row in matrix for x in row):raise ValueError('Finite Hessian required')
            for point in points:
                if (not point['SCF_converged'] or len(point['gradient_Ha_Bohr'])!=len(center)
                        or any(len(atom)!=3 for atom in point['gradient_Ha_Bohr'])):raise ValueError('Converged matched gradients required')
                if not math.isfinite(point['energy_Ha']) or any(not math.isfinite(x) for x in flatten(point['gradient_Ha_Bohr'])):raise ValueError('Finite scalar/force data required')
            centers=[p for p in points if p['coordinates_Bohr']==center]
            if len(centers)!=1:raise ValueError('One exact center required')
            for atom in range(len(center)):
                for axis in range(3):
                    for step in (.0005,.001):
                        pair=[]
                        for sign in (-1,1):
                            target=[list(x) for x in center];target[atom][axis]+=sign*step
                            matches=[p for p in points if p['coordinates_Bohr']==target]
                            if len(matches)!=1:raise ValueError('One exact point for every Cartesian displacement required')
                            pair.append(matches[0])
                        minus,plus=pair
                        finite=[(a-b)/(2*step) for a,b in zip(flatten(plus['gradient_Ha_Bohr']),flatten(minus['gradient_Ha_Bohr']))]
                        delta=[v-matrix[i][3*atom+axis] for i,v in enumerate(finite)]
                        relative=norm(delta)/max(norm(finite),1e-14)
                        error=(plus['energy_Ha']-minus['energy_Ha'])/(2*step)-centers[0]['gradient_Ha_Bohr'][atom][axis]
                        if relative>.01 or abs(error)>2e-6:raise ValueError('Original derivative threshold failed')
                        rows.append(dict(device=case['device'],control=case['control'],variant=variant,
                            atom=atom,axis=axis,step_Bohr=step,relative_Hessian_error=relative,energy_gradient_error_Ha_Bohr=error))
    return dict(cases=len(cases),checks=len(rows),columns=rows,
        maximum_relative_Hessian_error=max(r['relative_Hessian_error'] for r in rows),
        maximum_energy_gradient_error_Ha_Bohr=max(abs(r['energy_gradient_error_Ha_Bohr']) for r in rows),
        native_G_qualified=False,quantum_calculations_repeated=False)
