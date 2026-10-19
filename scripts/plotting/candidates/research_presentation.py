"""Insight-led main talk, refreshed from the same audited feed as the backup tables."""
import json
import math
from pathlib import Path


def number(value,digits=1):
    if type(value) not in (int,float) or not math.isfinite(value):return 'pending'
    return f'{value:.{digits}f}'


def research_frames(payload,root):
    folder=Path(root)/'results/candidates/figures'
    manifest=json.loads((folder/'metadata/manifest.json').read_text())
    if manifest['snapshot']!=payload['snapshot_path']:
        raise RuntimeError('Main talk and comparison figures require the same audited snapshot')
    rows={r['id']:r for r in payload['original39']+payload['additional_candidates']+payload.get('comparison_controls',[])}
    highlights=json.loads((folder/'data/highlights.json').read_text())
    proposals=json.loads((Path(root)/'results/discovery/vinyl_monomers/proposals.json').read_text())['records']
    capacities=[rows[gid]['capacity_nominal_mAh_g'] for gid in highlights['highlighted_ids']]
    ndi=next((r for r in highlights['matched_contrasts'] if r['core']=='ndi'),None)
    topology_contrast=ndi is not None and ndi['id']=='hyp_cata_ndi_ethyl' and ndi['baseline']=='design_ndi_ethyl'
    # Conclusions are conditional on the actual maintained selection, never assumed.
    quinone_leaders=all(rows[gid]['family']=='quinone (n-type)' for gid in highlights['highlighted_ids'][:2])
    capacity_title='Quinone models lead nominal capacity' if quinone_leaders else 'Capacity highlights the lightest two-electron models'
    pieces=[r'\begin{frame}[plain]\titlepage',
            r'\begin{center}\large Molecular design $\rightarrow$ polymerizable redox centers\end{center}',r'\end{frame}',
            r'\begin{frame}{Start with the monomer synthesis route}',r'\centering',
            r'\begin{tikzpicture}[node distance=.75cm,>=Latex]',
            r'\node[draw=teal,fill=soft,rounded corners,text width=3.1cm,minimum height=1.5cm,align=center] (a) {Redox center\\+ \mbox{vinylbenzyl}\\handle};',
            r'\node[right=of a,draw=teal,fill=soft,rounded corners,text width=3.1cm,minimum height=1.5cm,align=center] (b) {Functionalized\\vinyl monomer};',
            r'\node[right=of b,draw=teal,fill=soft,rounded corners,text width=3.1cm,minimum height=1.5cm,align=center] (c) {Polymer bearing\\redox centers};',
            r'\draw[->,thick,teal] (a)--(b);\draw[->,thick,teal] (b)--(c);\end{tikzpicture}',
            r'\vspace{1.2em}\begin{block}{The design question}',
            r'Which center combines useful redox behavior with a workable route to the polymer?',r'\end{block}',
            r'\vspace{.5em}Potential and capacity guide screening; kinetics and synthesis guide the next tests.',r'\end{frame}']

    background=[r'\begin{frame}{The target is a redox-active polymer for a flow battery}',
        r'\begin{columns}[T]\begin{column}{.48\textwidth}',
        r'\begin{block}{Research objective}Choose pendant redox centers that remain useful after monomer functionalization and polymerization.\end{block}',
        r'\vspace{.8em}\large Which molecular advantages survive the route to the polymer?',
        r'\end{column}\begin{column}{.48\textwidth}',
        r'\begin{tikzpicture}[every node/.style={draw=teal,fill=soft,rounded corners,align=center,text width=2.5cm,minimum height=1.2cm,font=\small}]',
        r'\node (a) {Redox potential\\Operating window};',
        r'\node[right=.3cm of a] (b) {Active loading\\Capacity};',
        r'\node[below=.4cm of a] (c) {Charge transfer\\Kinetics};',
        r'\node[below=.4cm of b] (d) {Chemical route\\Redox retention};',r'\end{tikzpicture}',
        r'\end{column}\end{columns}',
        r'\vfill\footnotesize Molecular screening supplies design hypotheses; polymer accessibility and cycling must be tested.',r'\end{frame}']
    # The main talk opens with background, then the confirmed synthesis route.
    pieces[3:3]=background
    pieces.extend([r'\begin{frame}{Workflow: from candidate families to synthesis controls}',r'\centering',
        r'\begin{tikzpicture}[node distance=.45cm,>=Latex,every node/.style={font=\small}]',
        r'\node[draw=teal,fill=soft,rounded corners,text width=2.65cm,minimum height=1.6cm,align=center] (a) {Build exact\\molecular models};',
        r'\node[right=of a,draw=teal,fill=soft,rounded corners,text width=2.65cm,minimum height=1.6cm,align=center] (b) {Sample\\conformers};',
        r'\node[right=of b,draw=teal,fill=soft,rounded corners,text width=2.65cm,minimum height=1.6cm,align=center] (c) {Compute solution\\redox states};',
        r'\node[right=of c,draw=teal,fill=soft,rounded corners,text width=2.65cm,minimum height=1.6cm,align=center] (d) {Compare models\\and select\\controls};',
        r'\draw[->,thick,teal] (a)--(b);\draw[->,thick,teal] (b)--(c);\draw[->,thick,teal] (c)--(d);',r'\end{tikzpicture}',
        r'\vspace{1em}\begin{columns}[T]\begin{column}{.47\textwidth}',
        r'\textbf{Chemical space}\\Diimides; benzoquinones (BQ) and anthraquinones (AQ); viologens; positive-side controls.',
        r'\end{column}\begin{column}{.47\textwidth}',
        r'\textbf{Selection}\\Potential, nominal capacity, kinetic evidence and synthetic route.',
        r'\end{column}\end{columns}',
        r'\vfill\footnotesize Finished fields enter comparisons after their checks. Missing properties remain pending.',r'\end{frame}',
        r'\begin{frame}{Method: compare free energies of matched solution states}',
        r'\centering\begin{tikzpicture}[node distance=1.2cm,>=Latex]',
        r'\node[draw=teal,fill=soft,rounded corners,align=center,minimum width=2.3cm,minimum height=1.1cm] (a) {$M^{q}$\\$G_0$};',
        r'\node[right=of a,draw=teal,fill=soft,rounded corners,align=center,minimum width=2.3cm,minimum height=1.1cm] (b) {$M^{q-1}$\\$G_1$};',
        r'\node[right=of b,draw=teal,fill=soft,rounded corners,align=center,minimum width=2.3cm,minimum height=1.1cm] (c) {$M^{q-2}$\\$G_2$};',
        r'\draw[->,thick,teal] (a)--node[above] {$+e^-$}(b);\draw[->,thick,teal] (b)--node[above] {$+e^-$}(c);',r'\end{tikzpicture}',
        r'\vspace{.7em}\begin{columns}[T]\begin{column}{.55\textwidth}',
        r'\small\begin{tabular}{@{}ll@{}}\textbf{Geometry} & r2SCAN-D4 / SVP(D) \\ \textbf{Electronic energy} & $\omega$B97M-V / TZVPD \\ \textbf{Environment} & SMD acetonitrile \\ \textbf{Thermal correction} & Native Hessian, 298 K \\ \end{tabular}',
        r'\end{column}\begin{column}{.43\textwidth}',
        r'\small $G=E_{\rm electronic}+G_{\rm thermal}$\\[.6em]',
        r'Redox potentials from $\Delta G$, referenced consistently to Fc/Fc$^+$.\\[.6em]',
        r'Nominal capacity: $Q=\dfrac{nF}{3.6M}$.',r'\end{column}\end{columns}',
        r'\vfill\footnotesize Charge states preserve the declared molecular graph. Stationarity, derivatives, spin and reference consistency are checked separately.',r'\end{frame}'])

    def figure(title,name,insight,note):
        pieces.extend([r'\begin{frame}{'+title+'}',r'\centering',
            r'\includegraphics[width=.98\linewidth,height=.70\textheight,keepaspectratio]{../results/candidates/figures/'+name+'}',
            r'\par\smallskip\raggedright\textcolor{teal}{\textbf{'+insight+'}}',
            r'\par\smallskip\footnotesize '+note,r'\end{frame}'])

    figure(capacity_title,'plots/redox_capacity.pdf',
           ('The two BQ models reach '+number(capacities[0])+r' and '+number(capacities[1])+r' mAh g$^{-1}$.' if quinone_leaders else 'The highlighted models span '+number(min(capacities))+'--'+number(max(capacities))+r' mAh g$^{-1}$.'),
           str(manifest['candidates'])+' designs; '+str(manifest['admitted_potential_candidates'])+' checked ladders. Nominal capacity excludes polymer and counterions.')
    shift_insight='Matched contrasts separate core changes from substituent changes.'
    if ndi is not None:
        shift_insight=(r'NDI core topology shifts' if topology_contrast else r'The highlighted NDI contrast shifts')+r' $E_1$ by '+number(ndi['delta_E1_V'],3)+r' V and $E_2$ by '+number(ndi['delta_E2_V'],3)+' V.'
    figure('Core topology changes NDI reduction potentials' if topology_contrast else 'Matched contrasts reveal molecular design effects','plots/new_candidate_shifts.pdf',shift_insight,
           'Matched model differences in MeCN; voltage shifts alone do not establish battery performance.')
    figure('SA is a separate filter for candidate choice','plots/synthesis_tradeoff.pdf',
           'High capacity and low SA do not establish a polymerization route.',
           'SA is heuristic; polymerization compatibility and redox retention remain untested.')
    figure('Capacity leaders and the first monomer pilot differ','structures/candidate_structures.pdf',
           'The pilot tests NDI, PMDI and viologen with the actual vinyl handle.',
           'NDI: naphthalene diimide. PMDI: pyromellitic diimide. Vinyl potentials remain pending.')

    pieces.extend([r'\begin{frame}{The vinyl group carries a small nominal capacity penalty}',
                   r'\centering\renewcommand{\arraystretch}{1.65}',
                   r'\begin{tabular}{@{}lrr@{}}\toprule',
                   r'Core & Capped model & Vinyl monomer \\',
                   r'& \multicolumn{2}{c}{Nominal capacity / mAh g$^{-1}$} \\ \midrule'])
    for gid,label in [('vm_design_ndi_ethyl','NDI / ethyl'),('vm_design_pmdi_ethyl','PMDI / ethyl'),('vm_ethylviologen','Viologen / ethyl')]:
        row=next(r for r in proposals if r['id']==gid);source=rows[row['source_capped_id']]
        pieces.append(label+' & '+number(source['capacity_nominal_mAh_g'])+' & '+number(row['capacity_nominal_mAh_g'])+r' \\')
    pieces.extend([r'\bottomrule\end{tabular}',
                   r'\vspace{1em}\begin{block}{The larger open question is electronic}',
                   r'How much does the vinyl-bearing structure shift the redox potentials?',r'\end{block}',
                   r'\footnotesize Same declared two-electron count; one added carbon. Stoichiometric capacities exclude counterions and polymer. Vinyl potentials remain pending.',r'\end{frame}',
                   r'\begin{frame}{Three findings guide the next experiments}',
                   r'\begin{block}{Capacity}\textbf{'+('Substituted BQ models lead the current nominal-capacity screen.' if quinone_leaders else 'The highlighted designs expose the nominal capacity tradeoff.')+r'}\end{block}',
                   r'\begin{block}{Voltage tuning}\textbf{'+('NDI core topology provides a clear matched computational contrast.' if topology_contrast else 'Matched contrasts identify molecular changes for follow-up.')+r'}\end{block}',
                   r'\begin{block}{Synthesis direction}\textbf{Compare exact vinyl NDI, PMDI and viologen before committing to polymer synthesis.}\end{block}',
                   r'\vfill\footnotesize Next decision: combine monomer redox properties, qualified kinetics and polymerization evidence. Solubility and cycling performance remain unmeasured for these exact targets.',r'\end{frame}'])
    return '\n'.join(pieces)+'\n'
