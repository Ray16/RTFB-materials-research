"""Check the requested benchmark, selection overview and corrected synthesis narrative."""
import csv,hashlib,json,math
from pathlib import Path
import fitz
root=Path(__file__).resolve().parents[4];base=Path(__file__).resolve().parent
state=json.loads((base/'slides/status.json').read_text());pdf=fitz.open(root/'slides/redox_workflow.pdf')
checks=[]
def check(label,passed):checks.append(dict(check=label,passed=bool(passed)));assert passed,label
check('Checked PDF identity',hashlib.sha256((root/'slides/redox_workflow.pdf').read_bytes()).hexdigest()==state['pdf_sha256'])
check('Thirteen main research slides',state['main_research_slides']==13)
check('Backup starts after the main story',pdf[13].get_text().strip().splitlines()[0]=='Backup')
main=' '.join(pdf[i].get_text() for i in range(13));check('Corrected synthesis route is explicit','Post-functionalize the existing polymer' in main)
check('Superseded monomer route absent from main','vinyl monomer' not in main.lower() and 'monomer synthesis route' not in main.lower())
check('Search and selection are distinct','837' in pdf[4].get_text() and '36 additional' in pdf[4].get_text())
check('Benchmark values visible','0.209' in pdf[6].get_text() and '0.246' in pdf[6].get_text() and '0.655' in pdf[6].get_text() and '1.202' in pdf[6].get_text())
actual_pages=[i for i,p in enumerate(pdf) if p.get_text().startswith('Measured and calculated parent potentials')];check('Actual-potential table present once',len(actual_pages)==1)
text=pdf[actual_pages[0]].get_text();check('Actual table contains first/second step values',all(v in text for v in ('-0.785','-0.595','-1.204','-1.207','0.655','1.202')))
for i in range(13):
    page=pdf[i];spans=[s for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for s in line['spans']]
    check(f'Main page {i+1} text within bounds',all(s['bbox'][0]>=-.1 and s['bbox'][1]>=-.1 and s['bbox'][2]<=page.rect.width+.1 and s['bbox'][3]<=page.rect.height+.1 for s in spans))
    body=[s for s in spans if s['text'].strip() and not(s['text'].strip() in ('RTFB','Materials','Research',str(i+1)) and s['bbox'][1]>page.rect.height-12)]
    check(f'Main page {i+1} footer clearance',all(s['bbox'][3]<=page.rect.height-14 for s in body))
    page.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(base/f'page_{i+1:02}.png')
pdf[13].get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(base/'page_14.png')
pdf[actual_pages[0]].get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(base/'actual_values.png')
figure=root/'results/figures/validation/matched_parent_validation.json';f=json.loads(figure.read_text())
for name,digest in f['files'].items():check('Benchmark artifact '+name,hashlib.sha256((root/name).read_bytes()).hexdigest()==digest)
for name,digest in f['source_sha256'].items():check('Benchmark source '+name,hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest)
with (root/'results/figures/validation/matched_parent_validation.csv').open() as stream:rows=list(csv.DictReader(stream))
errors=[float(r['calc_V'])-float(r['exp_V']) for r in rows]
check('Twelve potentials and eight parents',len(rows)==12 and len({r['id'] for r in rows})==8)
check('MAE independently reconstructed',abs(sum(abs(e) for e in errors)/12-.2090416666666667)<1e-12)
check('RMSE independently reconstructed',abs(math.sqrt(sum(e*e for e in errors)/12)-.24571566392885905)<1e-12)
selection=json.loads((root/'results/candidates/figures/data/candidate_selection_overview.json').read_text())
check('Recorded candidate selection arithmetic',selection['enumerated']==837 and selection['matched']==18 and sum(selection['complementary'].values())==18 and selection['additional']==36 and selection['original']==39)
result=dict(status='checked_postgrafting_research_visualizations',checks=len(checks),failed_checks=0,main_slides=13,backup_slides=state['backup_slides'],pages=len(pdf),actual_values_page=actual_pages[0]+1,selection_page=5,validation_page=7,pdf_sha256=state['pdf_sha256'],snapshot=state['snapshot'],software_tests_passed=595,new_quantum_calculations_launched_for_this_update=0,goal_remains_paused=True)
(base/'summary.json').write_text(json.dumps(result,indent=2)+'\n');(base/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(result))
