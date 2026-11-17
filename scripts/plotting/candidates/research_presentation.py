"""Insight-led main talk, refreshed from the same audited feed as the backup tables."""
import json
import hashlib
import math
from pathlib import Path


def number(value,digits=1):
    if type(value) not in (int,float) or not math.isfinite(value):return 'pending'
    return f'{value:.{digits}f}'


def benchmark_frame(root):
    root = Path(root)
    include = root / 'slides/generated/retrospective_redox_benchmark.tex'
    metadata = root / 'results/figures/validation/retrospective_redox_benchmark.json'
    if include.exists() or metadata.exists():
        record = json.loads(metadata.read_text())
        if record['status'] != 'checked_retrospective_snapshot':
            raise ValueError('Expanded benchmark figure is not checked')
        for path, expected in record['output_sha256'].items():
            if hashlib.sha256((root / path).read_bytes()).hexdigest() != expected:
                raise ValueError('Expanded benchmark figure artifact changed')
        return include.read_text()
    return (root / 'slides/generated/matched_parent_validation_plot.tex').read_text()


def research_frames(payload,root):
    folder=Path(root)/'results/candidates/figures'
    manifest=json.loads((folder/'metadata/manifest.json').read_text())
    if manifest['snapshot']!=payload['snapshot_path']:
        raise RuntimeError('Main talk and comparison figures require the same audited snapshot')
    rows={r['id']:r for r in payload['original39']+payload['additional_candidates']+payload.get('comparison_controls',[])}
    highlights=json.loads((folder/'data/highlights.json').read_text())
    assessment_path=Path(root)/'results/grafting/selection_assessment.json'
    assessment=json.loads(assessment_path.read_text()) if assessment_path.exists() else {}
    bq_hold=any(v.get('status')=='on_hold' for key,v in assessment.get('candidates',{}).items() if key.startswith('exp_bq_'))
    bq_ids=('exp_bq_oxy_methoxy','exp_bq_oxy_trimethyl')
    bq_decisions=[assessment.get('candidates',{}).get(key,{}) for key in bq_ids]
    bq_base=all(v.get('status')=='base_required' and v.get('thermodynamic_status')=='favorable_base_assisted_model' for v in bq_decisions)
    base_note=('Base required; BQ retained. DBU: '+ '/'.join(number(v['base_assisted_deltaG_kJ_mol']['dbu']['free_ions_kJ_mol']) for v in bq_decisions)
               +r' kJ mol$^{-1}$ (MeCN model); selectivity and yield untested.') if bq_base else None
    capacities=[rows[gid]['capacity_nominal_mAh_g'] for gid in highlights['highlighted_ids']]
    ndi=next((r for r in highlights['matched_contrasts'] if r['core']=='ndi'),None)
    topology_contrast=ndi is not None and ndi['id']=='hyp_cata_ndi_ethyl' and ndi['baseline']=='design_ndi_ethyl'
    # Conclusions are conditional on the actual maintained selection, never assumed.
    quinone_leaders=all(rows[gid]['family']=='quinone (n-type)' for gid in highlights['highlighted_ids'][:2])
    capacity_title='Quinone models lead nominal capacity' if quinone_leaders else 'Capacity highlights the lightest two-electron models'
    pieces=[r'\begin{frame}[plain]\titlepage',
            r'\begin{center}\large Molecular design $\rightarrow$ polymer-grafted redox centers\end{center}',r'\end{frame}',
            r'\begin{frame}{Post-functionalize the existing polymer}',r'\centering',
            r'\begin{tikzpicture}[node distance=.75cm,>=Latex]',
            r'\node[draw=teal,fill=soft,rounded corners,text width=3.1cm,minimum height=1.5cm,align=center] (a) {Preformed polymer\\with reactive sites};',
            r'\node[right=of a,draw=teal,fill=soft,rounded corners,text width=3.1cm,minimum height=1.5cm,align=center] (b) {Attach \mbox{redox} groups\\at polymer sites};',
            r'\node[right=of b,draw=teal,fill=soft,rounded corners,text width=3.1cm,minimum height=1.5cm,align=center] (c) {Polymer bearing\\redox centers};',
            r'\draw[->,thick,teal] (a)--(b);\draw[->,thick,teal] (b)--(c);\end{tikzpicture}',
            r'\vspace{1.2em}\begin{block}{The design question}',
            r'Which center combines useful redox behavior with a workable route to the polymer?',r'\end{block}',
            r'\vspace{.5em}Potential and molecular capacity guide screening; attachment and accessibility guide downselection.',r'\end{frame}']

    background=[r'\begin{frame}{The target is a redox-active polymer for a flow battery}',
        r'\begin{columns}[T]\begin{column}{.48\textwidth}',
        r'\begin{block}{Research objective}Choose pendant redox centers that remain useful after grafting onto the polymer.\end{block}',
        r'\vspace{.8em}\large Which molecular advantages survive the grafted polymer environment?',
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
        r'\node[draw=teal,fill=soft,rounded corners,text width=2.65cm,minimum height=1.6cm,align=center] (a) {Screen precursor\\and attachment\\route};',
        r'\node[right=of a,draw=teal,fill=soft,rounded corners,text width=2.65cm,minimum height=1.6cm,align=center] (b) {Assess grafting\\under defined\\conditions};',
        r'\node[right=of b,draw=teal,fill=soft,rounded corners,text width=2.65cm,minimum height=1.6cm,align=center] (c) {Compute solution\\redox properties};',
        r'\node[right=of c,draw=teal,fill=soft,rounded corners,text width=2.65cm,minimum height=1.6cm,align=center] (d) {Compare dimers,\\trimers and\\stacking};',
        r'\draw[->,thick,teal] (a)--(b);\draw[->,thick,teal] (b)--(c);\draw[->,thick,teal] (c)--(d);',r'\end{tikzpicture}',
        r'\vspace{1em}\begin{columns}[T]\begin{column}{.47\textwidth}',
        r'\textbf{Chemical space}\\Diimides; benzoquinones (BQ) and anthraquinones (AQ); viologens; positive-side controls.',
        r'\end{column}\begin{column}{.47\textwidth}',
        r'\textbf{Selection}\\Potential, nominal capacity, kinetic evidence and grafting feasibility.',
        r'\end{column}\end{columns}',
        r'\vfill\footnotesize New candidates advance after a reviewed grafting route. Unknown routes stay on hold; method controls remain separate.',r'\end{frame}',
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
        r'\vfill\footnotesize E1/E2 identify successive one-electron steps. Diagram: reductions. Geometry, spin and reference checks are retained.',r'\end{frame}'])

    method_position=next(i for i,piece in enumerate(pieces) if piece.startswith(r'\begin{frame}{Method:'))
    pieces[method_position:method_position]=[r'\begin{frame}{Candidate expansion uses matched molecular changes}',r'\centering',r'\includegraphics[width=.98\linewidth,height=.77\textheight,keepaspectratio]{../results/candidates/figures/plots/candidate_selection_overview.pdf}',r'\par\smallskip\footnotesize Current route: attach redox groups to a preformed polymer.',r'\end{frame}']
    pieces.append(benchmark_frame(root))

    def figure(title,name,insight,note):
        pieces.extend([r'\begin{frame}{'+title+'}',r'\centering',
            r'\includegraphics[width=.98\linewidth,height=.70\textheight,keepaspectratio]{../results/candidates/figures/'+name+'}',
            r'\par\smallskip\raggedright\textcolor{teal}{\textbf{'+insight+'}}',
            r'\par\smallskip\footnotesize '+note,r'\end{frame}'])

    figure(capacity_title,'plots/redox_capacity.pdf',
           ('The two BQ models reach '+number(capacities[0])+r' and '+number(capacities[1])+r' mAh g$^{-1}$.' if quinone_leaders else 'The highlighted models span '+number(min(capacities))+'--'+number(max(capacities))+r' mAh g$^{-1}$.'),
           (base_note if bq_base else 'BQ/AQ synthesis selection is on hold pending a base-assisted grafting assessment.' if bq_hold else str(manifest['candidates'])+' designs; '+str(manifest['admitted_potential_candidates'])+' checked ladders. Nominal capacity excludes polymer and counterions.'))
    shift_insight='Matched contrasts separate core changes from substituent changes.'
    if ndi is not None:
        shift_insight=(r'NDI core topology shifts' if topology_contrast else r'The highlighted NDI contrast shifts')+r' $E_1$ by '+number(ndi['delta_E1_V'],3)+r' V and $E_2$ by '+number(ndi['delta_E2_V'],3)+' V.'
    figure('Core topology changes NDI reduction potentials' if topology_contrast else 'Matched contrasts reveal molecular design effects','plots/new_candidate_shifts.pdf',shift_insight,
           'Matched model differences in MeCN; voltage shifts alone do not establish battery performance.')
    core_path=folder/'data/core_topology_comparison.json'
    core=json.loads(core_path.read_text()) if core_path.exists() else {}
    if (core.get('comparison_published') and core.get('same_formula')
            and {r['id'] for r in core.get('structures',[])}=={'design_ndi_ethyl','hyp_cata_ndi_ethyl'}):
        figure('Cata-NDI shifts potential without a nominal capacity penalty','structures/core_topology_comparison.pdf',
               'A core-isomer comparison worth discussing alongside the capacity leaders.',
               'Local graft models; precursor availability, attachment and polymer behavior remain untested.')
    else:
        figure('SA is a separate filter for candidate choice','plots/synthesis_tradeoff.pdf',
               'High capacity and low SA do not establish a workable grafting route.',
               'SA is heuristic; attachment selectivity and site accessibility require testing.')
    figure('Capacity leaders and representative graft-site controls','structures/candidate_structures.pdf',
           ('Base required; BQ candidates retained for evaluation.' if bq_base else 'BQ/AQ remain descriptor leaders; synthesis selection is on hold.' if bq_hold else 'Quinones, diimides and viologens provide complementary graft-site hypotheses.'),
           'NDI: naphthalenediimide; PMDI: pyromellitic diimide. Polymer effects unvalidated.')

    pieces.extend([r'\begin{frame}{Downselection couples redox behavior with grafting feasibility}',
                   r'\begin{columns}[T]\begin{column}{.31\textwidth}\begin{block}{Attach}Precursor availability\\[.6em]Attachment selectivity\\[.6em]Redox-center retention\end{block}\end{column}',
                   r'\begin{column}{.31\textwidth}\begin{block}{Load}Degree of grafting\\[.6em]Accessible active sites\\[.6em]Counterion mass\end{block}\end{column}',
                   r'\begin{column}{.31\textwidth}\begin{block}{Operate}Swelling and electrolyte\\[.6em]Charge transport\\[.6em]Cycling stability\end{block}\end{column}\end{columns}',
                   r'\vspace{1em}\textcolor{teal}{\textbf{High nominal molecular capacity is a screening descriptor; usable polymer capacity depends on accessible grafted sites.}}',
                   r'\vfill\footnotesize Local benzyl-capped calculations do not include the full polymer environment.',r'\end{frame}',
                   r'\begin{frame}{Three findings guide the next experiments}',
                   r'\begin{block}{Capacity}\textbf{'+('Substituted BQ models lead the current nominal-capacity screen.' if quinone_leaders else 'The highlighted designs expose the nominal capacity tradeoff.')+r'}\end{block}',
                   r'\begin{block}{Voltage tuning}\textbf{'+('NDI core topology provides a clear matched computational contrast.' if topology_contrast else 'Matched contrasts identify molecular changes for follow-up.')+r'}\end{block}',
                   r'\begin{block}{Synthesis direction}\textbf{'+('Base required: retain BQ candidates and test selective attachment.' if bq_base else 'BQ/AQ are on hold pending a defined, credible grafting route.' if bq_hold else 'Review quinone, diimide and viologen precursors for post-grafting.')+r'}\end{block}',
                   r'\vfill\footnotesize Next decision: combine molecular redox properties with attachment feasibility, graft loading and accessible redox sites. Solubility and cycling performance remain unmeasured for these exact targets.',r'\end{frame}'])
    return '\n'.join(pieces)+'\n'
