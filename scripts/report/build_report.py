#!/usr/bin/env python3
"""Build the portable evidence-led report; requires dvisvgm. No new calculations.

Inputs are deliberately pinned. Advance the evidence and narrative together.
Figures are isolated SVG images, so repeated PDF-generated IDs cannot collide.
"""
from pathlib import Path
import base64
import datetime
import hashlib
import html
import json
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from enrich import visuals, INPUTS

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = ROOT / 'results/report'
FIGDIR = ROOT / 'results/figures/manuscript/2026-11-03_six_chloride_contacts'
SOURCES = [
 ('Design-axis benchmark record', 'docs/DESIGN_AXES.md'),
 ('Exact viologen: checked reduction steps', 'validation/viologen_native_recovery_2026-10-30/README.md'),
 ('Exact viologen: attachment reaction', 'validation/viologen_grafting_2026-10-31/README.md'),
 ('Exact phenothiazine and candidate registry', 'validation/phenothiazine_registry_2026-11-03/README.md'),
 ('Four finite polymer cycles', 'validation/ionic_cycle_completion_2026-11-02/README.md'),
 ('Sampling completeness', 'validation/corrected_sampling_2026-11-02/README.md'),
 ('Publication scope and qualifications', 'docs/PUBLICATION_READINESS.md'),
 ('Current results guide', 'results/README.md'),
]
FIGURES = {
 'benchmark': ('Calculated and measured voltages', 'Historical 39-potential comparison: MAE 0.25 V. The 35-point tier-A subset summarized above has MAE 0.22 V; these are different benchmark scopes.'),
 'pmdi_capacity': ('How attachment and loading change stored charge per gram', 'PMDI accounting scenarios, not measured polymer performance. Dry neutral H-capped fragments, with residual chloride retained; no added electrolyte or crosslinker mass. One or two electrons per graft are assumed.'),
 'grafting': ('Attachment depends on the reaction conditions', 'Calculated attachment free energies under DBU, DIPEA and HCl reference conventions. Contact-pair and separated-ion treatments are distinct scenarios; neither determines yield.'),
 'viologen_interaction_definitions_final': ('Two definitions of neighboring-site effects', 'Four retained configurations. Frozen capped-pair energies and whole-chain occupancy cycles answer different questions; neither is an electronic-coupling measurement.'),
 'conditional_grafting': ('Does an occupied neighbor change attachment?', 'Electronic and free-energy changes are labeled separately. These finite configurations do not establish an average over a real polymer.'),
 'decomposition': ('Contact energy versus structural distortion', 'Within this segment decomposition, frozen contact energies are −44 to −55 kJ/mol, while distortion costs about 0.4–0.6 kJ/mol. This is a separate calculation from the charged-viologen occupancy cycles.'),
 'finite_ensemble': ('How much does the choice of structure matter?', 'Finite-set comparisons only. The occupied-structure sampling has not passed its convergence requirement.'),
 'viologen_native_redox': ('The exact viologen model: two reduction steps', 'Calculated local-model potentials versus ferrocene in acetonitrile. Numerical sensitivity is not experimental uncertainty.'),
}

def digest(path):
 return hashlib.sha256(path.read_bytes()).hexdigest()

def figure(name):
 title, caption = FIGURES[name]
 pdf = FIGDIR / (name + '.pdf')
 if not pdf.is_file():
  raise FileNotFoundError(pdf)
 with tempfile.TemporaryDirectory(dir=OUT) as td:
  dest = Path(td) / 'figure.svg'
  subprocess.run(['dvisvgm', '--pdf', '--no-fonts', str(pdf), '-o', str(dest)], check=True, capture_output=True)
  data = dest.read_bytes()
  root = ET.fromstring(data)
  if 'viewBox' not in root.attrib:
   raise ValueError(f'No viewBox: {name}')
 uri = 'data:image/svg+xml;base64,' + base64.b64encode(data).decode()
 return (f'<figure id="fig-{name}"><button class="figure-open" aria-label="Enlarge: {html.escape(title)}">'
         f'<img src="{uri}" alt="{html.escape(title + ". " + caption)}" loading="lazy"></button>'
         f'<figcaption><strong>{title}</strong> · {caption} <span class="enlarge">Click to enlarge ↗</span></figcaption></figure>')

def build():
 OUT.mkdir(parents=True, exist_ok=True)
 template = (HERE / 'report.html').read_text()
 for key, value in visuals().items():
  template = template.replace('{{'+key+'}}', value)
 names = re.findall(r'\{\{figure:([\w_]+)\}\}', template)
 for name in names:
  template = template.replace('{{figure:' + name + '}}', figure(name))
 sources = ''.join(f'<li id="source-{i}"><a href="../../{p}">{html.escape(label)}</a><small>{p}</small></li>' for i,(label,p) in enumerate(SOURCES, 1))
 stamp = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
 template = template.replace('{{sources}}', sources).replace('{{built}}', stamp)
 if '{{' in template: raise ValueError('Unresolved placeholder')
 output = OUT / 'project_highlights.html'
 output.write_text(template)
 paths = INPUTS + [HERE / 'build_report.py', HERE / 'report.html'] + [ROOT / p for _,p in SOURCES] + [FIGDIR / (n+'.pdf') for n in names]
 manifest = {'built_utc': stamp, 'evidence_scope': 'Pinned local sources; historical benchmark summaries and checked finite-model results, not new physical validation.', 'figure_snapshot': str(FIGDIR.relative_to(ROOT)), 'git_head': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(), 'generator_sha256':digest(HERE/'build_report.py'), 'inputs': {str(p.relative_to(ROOT)): digest(p) for p in paths}, 'output_sha256':digest(output), 'figures': names}
 (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(f'Wrote {output} ({output.stat().st_size/1024:.0f} KiB), {len(names)} vector figures')

if __name__ == '__main__': build()
