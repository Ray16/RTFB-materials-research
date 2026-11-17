"""Browser checks for the generated local report (Playwright required)."""
from pathlib import Path
import json
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/report'
with sync_playwright() as p:
 browser=p.chromium.launch(headless=True,args=['--no-sandbox'])
 page=browser.new_page(viewport={'width':1440,'height':1000},device_scale_factor=1)
 errors=[]
 page.on('pageerror',lambda e: errors.append(str(e)))
 page.goto((OUT/'project_highlights.html').as_uri())
 for d in page.locator('details').all(): d.evaluate('(e)=>e.open=true')
 for image in page.locator('figure img').all():
  image.scroll_into_view_if_needed()
  image.evaluate('(im)=>im.decode()')
 assert page.locator('figure img').count()==8
 for image in page.locator('.molecule').all(): image.evaluate('(im)=>im.decode()')
 assert page.locator('.molecule').count()==3
 page.get_by_role('button',name='TBAPF₆ mass scenario',exact=True).click()
 assert 'Viologen' in page.locator('#pareto-summary').inner_text()
 assert page.locator('#pareto-table tbody tr').count()==12
 page.get_by_role('button',name='LiPF₆ mass scenario',exact=True).click()
 assert 'Viologen' not in page.locator('#pareto-summary').inner_text()
 page.locator('#handle-count').select_option('5')
 page.locator('#graft-count').select_option('5')
 assert page.locator('#capacity-value').inner_text()=='148.58'
 page.locator('#electron-count').select_option('1')
 assert page.locator('#capacity-value').inner_text()=='74.29'
 page.locator('#handle-count').select_option('0')
 assert page.locator('#capacity-value').inner_text()=='0.00'
 page.locator('#handle-count').select_option('3')
 page.locator('#graft-count').select_option('3')
 page.locator('#electron-count').select_option('2')
 assert page.locator('#capacity-value').inner_text()=='124.52'
 assert page.evaluate('new Set([...document.querySelectorAll("[id]")].map(e=>e.id)).size===document.querySelectorAll("[id]").length')
 assert page.evaluate('[...document.querySelectorAll("a")].filter(a=>a.getAttribute("href").startsWith("#")).every(a=>document.getElementById(a.hash.slice(1)))')
 page.get_by_role('button',name='Electronic energy',exact=True).click()
 assert page.locator('.bar-row .value').all_text_contents()==['−3.46','+18.68','+1.20','−13.48']
 page.get_by_role('button',name='Free energy',exact=True).click()
 assert page.locator('.bar-row .value').all_text_contents()==['−3.94','+22.58','+8.28','−5.30']
 page.locator('.figure-open').first.click()
 assert page.locator('dialog').evaluate('(d)=>d.open')
 page.keyboard.press('Escape')
 assert not page.locator('dialog').evaluate('(d)=>d.open')
 for d in page.locator('details').all(): d.evaluate('(e)=>e.open=false')
 page.evaluate('window.scrollTo(0,0)');page.wait_for_timeout(350);page.screenshot(path=str(OUT/'desktop.png'),full_page=True);page.screenshot(path=str(OUT/'hero.png'))
 page.locator('#neighbors').scroll_into_view_if_needed();page.wait_for_timeout(350);page.screenshot(path=str(OUT/'neighbors.png'))
 page.locator('nav').evaluate("e=>e.style.visibility='hidden'")
 for target in ['shortlist','pareto','capacity-explorer']:
  page.locator('#'+target).screenshot(path=str(OUT/(target+'.png')))
 page.locator('nav').evaluate("e=>e.style.visibility=''")
 for width in [390,768,1440]:
  page.set_viewport_size({'width':width,'height':900})
  assert page.evaluate('document.documentElement.scrollWidth<=innerWidth'),f'Overflow at {width}'
  if width==390:
   page.evaluate('window.scrollTo(0,0)');page.screenshot(path=str(OUT/'mobile.png'),full_page=True)
 assert not errors,errors
 browser.close()
(OUT/'checks.json').write_text(json.dumps({'passed':True,'widths':[390,768,1440],'figures_decoded':8,'checks':['unique IDs','anchor targets','no page overflow','energy toggle values','figure dialog and Escape','no JavaScript exceptions','3 exact molecular images decoded','historical salt-scenario front flags','capacity zero/full/1e/2e scenarios']},indent=2)+'\n')
print('PASS: figures, anchors, chart, modal, desktop/tablet/mobile layout')
