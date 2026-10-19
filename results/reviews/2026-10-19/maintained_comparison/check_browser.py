"""Exercise the generated local comparison, data coverage and vector exports."""
import csv
import hashlib
import json
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright
from pypdf import PdfReader
from PIL import Image

ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
FIGURES=ROOT/'results/candidates/figures'
sys.path.insert(0,str(ROOT/'scripts/plotting/candidates'))
from collection_layout import artifact_path
def artifact(name):return FIGURES/artifact_path(name)
manifest=json.loads((artifact('manifest.json')).read_text())
records=json.loads((artifact('all_candidates.json')).read_text())
checks=[]


def check(name, passed):
    checks.append(dict(check=name,passed=bool(passed)))
    assert passed,name


check('Every registered candidate appears exactly once',len(records)==manifest['candidates']==len({r['id'] for r in records}))
check('Both cohorts retain their full counts',sum(r['cohort']=='Original' for r in records)==manifest['original'] and sum(r['cohort']=='Additional' for r in records)==manifest['additional'])
with (artifact('all_candidates.csv')).open() as stream:
    csv_rows=list(csv.DictReader(stream))
check('CSV and browser cover identical IDs',{r['id'] for r in csv_rows}=={r['id'] for r in records})
for name,digest in manifest['files'].items():
    check('Output checksum '+name,hashlib.sha256((FIGURES/name).read_bytes()).hexdigest()==digest)
for name,digest in manifest['producer_sha256'].items():
    check('Unchanged figure producer/input '+name,hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest)
for name in ('redox_capacity','synthesis_tradeoff','new_candidate_shifts','coverage','candidate_structures'):
    image=Image.open(artifact(f'{name}.png'))
    check(name+' PNG has 600-dpi metadata',min(image.info['dpi'])>599)
    document=PdfReader(artifact(f'{name}.pdf'))
    check(name+' vector PDF has double-column width',abs(float(document.pages[0].mediabox.width)-7.2*72)<.05)
    check(name+' PDF retains vector graphics',not list(document.pages[0].images))
check('Supplementary PDF paginates every candidate',len(PdfReader(artifact('redox_ladders.pdf')).pages)==(len(records)+17)//18)
highlights=json.loads((artifact('highlights.json')).read_text())
plate=json.loads((artifact('candidate_structures.json')).read_text())
proposal=json.loads((ROOT/'results/discovery/vinyl_monomers/proposals.json').read_text())
model_source={r['id']:r['exact_grafted_smiles'] for r in records}
vinyl_source={r['id']:r['smiles'] for r in proposal['records']}
check('The molecular plate has three models and three vinyl proposals',len(plate['structures'])==6)
for cell in plate['structures']:
    source=model_source if cell['scope']=='computed_capped_local_model' else vinyl_source
    check('Structure plate exact graph '+cell['id'],cell['smiles']==source[cell['id']])
lookup={r['id']:r for r in records}
eligible=[lookup[gid] for gid in highlights['eligible_ids']]
front=set(highlights['partial_front_ids'])
for row in eligible:
    vector=lambda r:(-r['E1_vs_Fc_V'],-r['E2_vs_Fc_V'],r['capacity_nominal_mAh_g'])
    dominated=any(all(a>=b for a,b in zip(vector(other),vector(row))) and
                  any(a>b for a,b in zip(vector(other),vector(row))) for other in eligible)
    check('Independent descriptor dominance '+row['id'],(row['id'] in front)==(not dominated))
check('Missing kinetics do not become a qualified primary front',not highlights['qualified_primary_Pareto_front_published'])
check('SA remains an annotation',highlights['SA_is_annotation_not_objective'])
with sync_playwright() as playwright:
    executable=Path('/nfs/lambda_stor_01/homes/rzhu/.cache/ms-playwright/chromium-1208/chrome-linux64/chrome')
    browser=playwright.chromium.launch(headless=True,executable_path=str(executable),args=['--no-sandbox'])
    page=browser.new_page(viewport=dict(width=1440,height=1050),device_scale_factor=1)
    errors=[];page.on('pageerror',lambda error:errors.append(str(error)))
    page.goto((FIGURES/'index.html').as_uri(),wait_until='networkidle')
    check('All candidates render in the table',page.locator('#body tr').count()==len(records))
    page.locator('#search').fill('hyp_cata_ndi_ethyl')
    check('Search finds the exact candidate',page.locator('#body tr').count()==1)
    page.locator('#body input').check()
    check('Structure comparison renders an SVG',page.locator('#cards .card svg:visible').count()==1)
    page.locator('#cards summary').filter(has_text='Proposed vinyl monomer').click()
    check('Vinyl structure can be inspected separately',page.locator('#cards .card svg:visible').count()==2)
    page.locator('#cards summary').filter(has_text='Proposed vinyl monomer').click()
    page.locator('#search').fill('')
    page.locator('#cohort').select_option('Additional')
    check('Additional cohort filter retains all new candidates',page.locator('#body tr').count()==manifest['additional'])
    page.locator('#cohort').select_option('')
    page.locator('#status').select_option('Pending / excluded / source changed')
    check('Pending candidates remain inspectable',page.locator('#body tr').count()==len(records)-manifest['admitted_potential_candidates'])
    page.locator('#status').select_option('')
    page.locator('th[data-key="capacity_nominal_mAh_g"]').click()
    actual=float(page.locator('#body tr').first.locator('td').nth(7).inner_text())
    check('Numeric sorting orders capacities correctly',abs(actual-min(r['capacity_nominal_mAh_g'] for r in records))<.051)
    page.locator('th[data-key="capacity_nominal_mAh_g"]').click()
    actual=float(page.locator('#body tr').first.locator('td').nth(7).inner_text())
    check('Numeric sorting reverses correctly',abs(actual-max(r['capacity_nominal_mAh_g'] for r in records))<.051)
    page.locator('#clear').click()
    if manifest['interactive_plotly']:
        check('Interactive plot has rendered points',page.locator('#redox-chart .scatterlayer .point').count()>0)
        plotted=page.evaluate("document.getElementById('redox-chart').data.map(t=>({axis:t.xaxis||'x',ids:t.customdata.map(r=>r[1]),x:t.x,y:t.y}))")
        for axis,key in [('x','E1_vs_Fc_V'),('x2','E2_vs_Fc_V')]:
            values={gid:(x,y) for t in plotted if t['axis']==axis for gid,x,y in zip(t['ids'],t['x'],t['y'])}
            expected={r['id']:(r[key],r['capacity_nominal_mAh_g']) for r in records if r['id'] in highlights['eligible_ids']}
            check('Every plotted endpoint matches the audited table '+key,values==expected)
        page.locator('#redox-chart .scatterlayer .point').first.click(force=True)
        check('Clicking a plot point opens the corresponding structure',page.locator('#cards .card svg:visible').count()==1)
        page.locator('#clear').click()
    page.evaluate('window.scrollTo(0,0)')
    page.screenshot(path=str(OUT/'desktop.png'),full_page=False)
    page.set_viewport_size(dict(width=390,height=844))
    page.wait_for_timeout(400)
    page.evaluate('window.scrollTo(0,0)')
    page.screenshot(path=str(OUT/'mobile.png'),full_page=False)
    check('Mobile page has no horizontal document overflow',page.evaluate('document.documentElement.scrollWidth <= window.innerWidth+2'))
    check('No browser JavaScript errors',not errors)
    browser.close()
result=dict(status='checked_data_figures_and_browser',checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),snapshot=manifest['snapshot'],records=len(records),details=checks)
(OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='details'}))
