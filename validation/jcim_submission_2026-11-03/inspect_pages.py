import hashlib,json,subprocess
from pathlib import Path
import fitz
from PIL import Image,ImageDraw
root=Path.cwd();out=root/'calcs/validation/jcim_submission_visual';out.mkdir(parents=True,exist_ok=True)
report={}
for name in ('manuscript','supporting_information'):
    p=root/'manuscript'/f'{name}.pdf';old=fitz.open(stream=subprocess.check_output(['git','show',f'HEAD:manuscript/{name}.pdf']),filetype='pdf');new=fitz.open(p)
    changed=[]
    for i,page in enumerate(new):
        pix=page.get_pixmap(matrix=fitz.Matrix(1,1),alpha=False)
        prev=old[i].get_pixmap(matrix=fitz.Matrix(1,1),alpha=False) if i<len(old) else None
        if prev is None or pix.samples!=prev.samples:
            changed.append(i+1)
            page.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False).save(out/f'{name}_{i+1:02d}.png')
    for start in range(0,len(changed),6):
        chunk=changed[start:start+6];canvas=Image.new('RGB',(1020,520*((len(chunk)+2)//3)),'#cccccc');draw=ImageDraw.Draw(canvas)
        for j,n in enumerate(chunk):
            im=Image.open(out/f'{name}_{n:02d}.png');im.thumbnail((330,490));x=(j%3)*340;y=(j//3)*520;canvas.paste(im,(x,y+25));draw.text((x+10,y+6),f'{name} {n}',fill='black')
        canvas.save(out/f'{name}_contact_{start//6+1}.png')
    report[name]={'pages':len(new),'previous_pages':len(old),'changed_pages':changed,'pdf_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(root/'validation/jcim_submission_2026-11-03/page_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
