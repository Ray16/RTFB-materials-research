"""Check the real local browser, filtering, SVG paths and vector reaction PDF."""
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():
    host_capacity(1)
    os.sched_setaffinity(0,cpu_allocation(1))
sys.path.insert(0,'/tmp/redox-candidate-comparison-deps')
from playwright.sync_api import sync_playwright

out=Path(__file__).resolve().parent
browser_path=ROOT/'results/grafting/index.html'
checks=[]
with sync_playwright() as p:
    browser=p.chromium.launch(headless=True,args=['--disable-gpu','--no-sandbox'])
    page=browser.new_page(viewport={'width':1440,'height':1100},device_scale_factor=1)
    errors=[]
    page.on('pageerror',lambda error:errors.append(str(error)))
    page.goto(browser_path.as_uri(),wait_until='networkidle')
    assert page.locator('article').count()==75
    checks.append('All 75 candidate cards present')
    targets=page.locator('img').evaluate_all('(imgs)=>imgs.map(i=>i.getAttribute("src"))')
    assert all((browser_path.parent/target).exists() for target in targets)
    checks.append('All structure and figure image paths resolve')
    page.locator('#search').fill('viologen')
    expected=page.locator('article').evaluate_all('(cards)=>cards.filter(c=>c.dataset.search.toLowerCase().includes("viologen")).length')
    assert page.locator('article:visible').count()==expected and expected>0
    checks.append('Search filters real candidate cards')
    page.locator('#search').fill('')
    assert page.locator('article:visible').count()==75
    page.evaluate('window.scrollTo(0,0)')
    page.screenshot(path=str(out/'browser_desktop.png'))
    for width in (1440,390):
        page.set_viewport_size({'width':width,'height':1100})
        page.evaluate('window.scrollTo(0,0)')
        assert page.evaluate('document.documentElement.scrollWidth<=window.innerWidth+1')
        checks.append(f'No horizontal overflow at {width}px')
    page.screenshot(path=str(out/'browser_mobile.png'))
    assert not errors,errors
    checks.append('No browser JavaScript errors')
    browser.close()
pdf=(ROOT/'results/grafting/figures/reaction_schemes.pdf').read_bytes()
scheme=(ROOT/'results/grafting/figures/reaction_schemes.svg').read_text()
assert pdf.startswith(b'%PDF-') and b'/Subtype /Image' not in pdf
assert scheme.count('<ns0:path')+scheme.count('<path')>20 and '<image' not in scheme
checks.append('Reaction PDF contains vector molecular paths and no raster image')
svg=list((ROOT/'results/grafting/structures').glob('*.svg'))
assert len(svg)==75 and all('<path' in p.read_text() for p in svg)
checks.append('75 exact-structure SVGs contain vector paths')
(out/'visual_checks.json').write_text(json.dumps(dict(passed=len(checks),failed=0,checks=checks),indent=2)+'\n')
print(len(checks),'visual checks passed')
