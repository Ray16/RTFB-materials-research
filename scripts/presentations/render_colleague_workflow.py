#!/usr/bin/env python
"""Render a compiled scientific deck for a separate assistant visual review."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from PIL import Image,ImageDraw


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pdf',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    subprocess.run(['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dNumRenderingThreads=1',
        '-sDEVICE=png16m','-r140','-sOutputFile='+str(a.out/'page_%02d.png'),str(a.pdf)],
        stdout=subprocess.DEVNULL,check=True,timeout=180)
    pages=sorted(a.out.glob('page_*.png'));contacts=[]
    for first in range(0,len(pages),8):
        subset=pages[first:first+8]
        canvas=Image.new('RGB',(1240,370*((len(subset)+1)//2)),'#DCE5EB')
        draw=ImageDraw.Draw(canvas)
        for i,page in enumerate(subset):
            im=Image.open(page).convert('RGB');im.thumbnail((600,337))
            x=10+(i%2)*620;y=20+(i//2)*370;canvas.paste(im,(x,y))
            draw.text((x,y-14),f'Slide {first+i+1}',fill='#15304A')
        name=a.out/f'contact_{first+1:02d}_{first+len(subset):02d}.png';canvas.save(name);contacts.append(str(name))
    (a.out/'render.json').write_text(json.dumps(dict(status='rendered_for_visual_review_not_yet_reviewed',
        pdf=str(a.pdf),pdf_sha256=hashlib.sha256(a.pdf.read_bytes()).hexdigest(),pages=len(pages),
        contacts=contacts,render_command='Ghostscript png16m, 140 dpi, one rendering thread'),indent=2)+'\n')
    print(json.dumps(dict(pages=len(pages),contacts=contacts),indent=2))


if __name__=='__main__':main()
