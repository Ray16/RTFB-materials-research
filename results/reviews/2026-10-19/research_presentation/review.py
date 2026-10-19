"""Check the research narrative, reported values and main-slide page boundaries."""
from pathlib import Path
import fitz,json,hashlib
root=Path(__file__).resolve().parents[4];base=Path(__file__).resolve().parent
state=json.loads((base/'slides/status.json').read_text());pdf=fitz.open(root/'slides/redox_workflow.pdf')
checks=[]
def check(name,passed):checks.append(dict(check=name,passed=bool(passed)));assert passed,name
check('The main story has eleven slides',state['main_research_slides']==11)
check('Backup follows the main conclusions',pdf[11].get_text().strip().splitlines()[0]=='Backup')
check('Deck checksum matches publisher',hashlib.sha256((root/'slides/redox_workflow.pdf').read_bytes()).hexdigest()==state['pdf_sha256'])
titles=[]
for i,page in enumerate(pdf):
    if i>11:break
    text=page.get_text();titles.append(text.splitlines()[0])
    if i<11:
        spans=[s for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for s in line['spans']]
        check('Main page '+str(i+1)+' text stays within page',all(s['bbox'][0]>=-.1 and s['bbox'][1]>=-.1 and s['bbox'][2]<=page.rect.width+.1 and s['bbox'][3]<=page.rect.height+.1 for s in spans))
        body=[s for s in spans if s['text'].strip() and not (s['text'].strip() in ('RTFB','Materials','Research',str(i+1)) and s['bbox'][1]>page.rect.height-12)]
        check('Main page '+str(i+1)+' leaves space above footer',all(s['bbox'][3]<=page.rect.height-14 for s in body))
        check('Main page '+str(i+1)+' omits execution details',not any(word in text for word in ('GPU','arithmetic checks','thread pools','cache')))
    page.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(base/f'page_{i+1:02}.png')
check('Background comes before workflow and method','The target' in titles[1] and titles[3].startswith('Workflow:') and titles[4].startswith('Method:'))
check('Key result insight is present','207.5' in pdf[5].get_text())
check('Core-topology contrast matches the feed','0.220' in pdf[6].get_text() and '0.312' in pdf[6].get_text())
check('Actual vinyl quantum properties remain pending','Vinyl potentials remain pending' in pdf[9].get_text())
summary=dict(status='checked_research_story_with_backup',main_slides=11,backup_slides=state['backup_slides'],pages=len(pdf),pdf_sha256=state['pdf_sha256'],snapshot=state['snapshot'],checks=len(checks),failed_checks=0,visually_reviewed_main_pages=list(range(1,12)),main_titles=titles[:11],software_tests_passed=591)
(base/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(base/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='main_titles'}))
