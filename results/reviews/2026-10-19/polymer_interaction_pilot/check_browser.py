"""Exercise offline 3D model selection without using a scientific GPU."""
from pathlib import Path
import sys
import json
import os

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
sys.path.insert(0,'/tmp/redox-candidate-comparison-deps')
from playwright.sync_api import sync_playwright

out=Path(__file__).resolve().parent
checks=[]
with sync_playwright() as p:
    browser=p.chromium.launch(headless=True,args=['--no-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
    page=browser.new_page(viewport={'width':1440,'height':1000})
    errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto((ROOT/'results/polymer_models/index.html').as_uri(),wait_until='networkidle')
    page.wait_for_function('typeof Plotly!=="undefined" && document.querySelectorAll("#viewer canvas").length>0')
    assert page.locator('#model option').count()>=9
    checks.append('Monomer and distinct dimer previews available')
    assert 'stacked' in page.locator('#model option:checked').inner_text()
    checks.append('Default view shows a dimer stacking arrangement')
    page.locator('#viewer').screenshot(path=str(out/'stacked_dimer_3d.png'))
    choices=page.locator('#model option').evaluate_all('(opts)=>opts.map(o=>({value:o.value,label:o.textContent}))')
    open_dimer=next(o for o in choices if '2 centers' in o['label'] and o['label'].endswith('open'))
    page.select_option('#model',open_dimer['value'])
    assert 'normal projection' in page.locator('#metrics').inner_text()
    checks.append('Model selector updates three-dimensional coordinates and contact metrics')
    page.locator('#viewer').screenshot(path=str(out/'open_dimer_3d.png'))
    assert 'DFT pending' in page.locator('#stage').inner_text() or 'stationary_DFT' in page.locator('#stage').inner_text()
    checks.append('Preview versus DFT stage is visible')
    for width in (1440,390):
        page.set_viewport_size({'width':width,'height':1000})
        page.wait_for_function('document.documentElement.scrollWidth<=window.innerWidth+1')
        assert page.evaluate('document.documentElement.scrollWidth<=window.innerWidth+1')
        checks.append(f'No horizontal overflow at {width}px')
    assert not errors,errors
    checks.append('No browser JavaScript errors')
    browser.close()
(out/'browser_checks.json').write_text(json.dumps(dict(checks=checks,passed=len(checks),failed=0),indent=2)+'\n')
print(len(checks),'offline three-dimensional browser checks passed')
