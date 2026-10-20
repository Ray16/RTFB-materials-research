"""Validate and publish the grafting Backup checkpoint while preserving the main story."""
from pathlib import Path
import json
import hashlib
import shutil
import fitz

root=Path(__file__).resolve().parents[4]
out=Path(__file__).resolve().parent
old=fitz.open(root/'slides/redox_workflow.pdf')
new=fitz.open(out/'redox_workflow.pdf')
assert len(new)==69
assert all(old[i].get_text()==new[i].get_text() for i in range(13))
checks=['All 13 main slides preserve their prior text']
assert new[13].get_text().startswith('Backup')
assert 'define the reaction' in new[14].get_text()
assert 'first computed results' in new[15].get_text()
checks.append('Grafting checkpoint follows Backup on pages 15 and 16')
for i in (14,15):
    page=new[i]
    spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
    assert all(s['bbox'][0]>=-.2 and s['bbox'][1]>=-.2 and s['bbox'][2]<=page.rect.width+.2 and s['bbox'][3]<=page.rect.height+.2 for s in spans)
    page.get_pixmap(matrix=fitz.Matrix(2,2)).save(out/f'slide_{i+1}.png')
    checks.append(f'Page {i+1} text stays within page bounds')
count=json.loads((out/'first_results/current.json').read_text())['status']['completed_grafting_deltaG']
assert str(count)+'/75' in new[15].get_text()
checks.append('Slide completion count matches the immutable checkpoint')
inputs=[root/'slides/redox_workflow.tex',root/'slides/generated/grafting_reference_2026-10-19.tex',out/'first_results/grafting_slide.pdf',out/'first_results/current.json']
status=dict(status='checked_grafting_backup_checkpoint',pages=69,main_research_slides=13,backup_slides=56,grafting_pages=[15,16],qualified_reactions_at_checkpoint=count,checks=checks,failed_checks=0,pdf_sha256=hashlib.sha256((out/'redox_workflow.pdf').read_bytes()).hexdigest(),source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
(out/'slide_checks.json').write_text(json.dumps(status,indent=2)+'\n')
old.close();new.close()
shutil.copyfile(out/'redox_workflow.pdf',root/'slides/.redox_workflow.grafting.checked.pdf')
(root/'slides/.redox_workflow.grafting.checked.pdf').replace(root/'slides/redox_workflow.pdf')
p=root/'results/LATEST_REVIEW.json'
review=json.loads(p.read_text())
review['current']['slides'].update(status=status['status'],pdf_sha256=status['pdf_sha256'],pages=69,backup_slides=56,grafting_checkpoint='results/reviews/2026-10-19/grafting_submission/slide_checks.json',grafting_pages=[15,16],grafting_reference_reactions_at_checkpoint=count)
p.write_text(json.dumps(review,indent=2)+'\n')
print(len(checks),'slide checks passed; checked PDF published')
