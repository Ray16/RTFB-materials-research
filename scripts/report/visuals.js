const reportData=JSON.parse(document.getElementById('visual-data').textContent);
const escapeText=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const names={d3_05stji:'Naphthoquinone lead',aq_benzylamino:'AQ · amino tether',pmdi:'PMDI',aq_benzyloxy:'AQ · oxygen tether',nq_benzyloxy:'NQ · oxygen tether',mophquinone:'Mophquinone',dmophquinone:'Dmophquinone',ndi_ammonium:'NDI ammonium',viologen:'Viologen',ethylviologen:'Ethylviologen',bisviologen:'Bisviologen',d3_80izcq:'d3_80izcq'};
function drawPareto(salt){
 const rows=reportData.historical.filter(r=>r.scenario===salt&&r.pool==='anolyte');
 const x=v=>80+(v-.4)/1.3*680,y=v=>370-v/200*310;
 let svg='<svg viewBox="0 0 850 455" role="img" aria-label="Historical voltage versus capacity comparison. Exact values and flags in the table below."><rect width="850" height="455" fill="white"/>';
 for(const v of [0,50,100,150,200])svg+=`<line x1="80" x2="760" y1="${y(v)}" y2="${y(v)}" stroke="#d7e1dc"/><text x="65" y="${y(v)+5}" text-anchor="end" font-size="14">${v}</text>`;
 for(const v of [.4,.8,1.2,1.6])svg+=`<text x="${x(v)}" y="395" text-anchor="middle" font-size="14">−${v.toFixed(1)}</text>`;
 svg+='<text x="420" y="430" text-anchor="middle" font-size="16">Average reduction potential (V vs ferrocene)</text><text transform="translate(23 220) rotate(-90)" text-anchor="middle" font-size="16">Nominal capacity (mAh/g)</text>';
 rows.forEach((r,i)=>{const px=x(-Number(r.E_V)),py=y(Number(r.capacity)),front=r.pareto_optimal==='True';
 svg+=`<g><title>${escapeText(names[r.id]||r.id)}: ${r.E_V} V, ${r.capacity} mAh/g; λ QC: ${escapeText(r.lambda_qc)}</title>`;
 svg+=front?`<text x="${px}" y="${py+8}" text-anchor="middle" fill="#96551e" font-size="28">★</text>`:`<circle cx="${px}" cy="${py}" r="6" fill="white" stroke="#51656b" stroke-width="2"/>`;
 if(front)svg+=`<text x="${px+14}" y="${py-9}" font-size="14" fill="#824c1d">${escapeText(names[r.id]||r.id)}</text>`;
 svg+='</g>';
 });
 document.getElementById('pareto-plot').innerHTML=svg+'</svg>';
 const fronts=rows.filter(r=>r.pareto_optimal==='True').map(r=>names[r.id]||r.id);
 document.getElementById('pareto-summary').textContent='Archived front flags in this scenario: '+fronts.join(' and ')+'. These are historical screening results, not approved materials.';
 document.getElementById('pareto-table').innerHTML='<table><thead><tr><th>Candidate</th><th>Mean voltage / V</th><th>Capacity / mAh/g</th><th>Rearrangement / eV</th><th>Electron-exchange cost / kJ/mol</th><th>Historical checks</th></tr></thead><tbody>'+rows.map(r=>`<tr><td>${escapeText(names[r.id]||r.id)}<small>${escapeText(r.id)}${r.pareto_optimal==='True'?' · ★ archived front':''}</small></td><td>${r.E_V}</td><td>${r.capacity}</td><td>${r.lambda_i_eV||'—'}</td><td>${r.dG_disp_kJmol||'—'}</td><td>λ: ${escapeText(r.lambda_qc)}<small>${escapeText(r.lambda_flags||'No flags recorded in this field')}</small></td></tr>`).join('')+'</tbody></table>';
}
document.querySelectorAll('[data-salt]').forEach(b=>b.addEventListener('click',()=>{document.querySelectorAll('[data-salt]').forEach(x=>x.setAttribute('aria-pressed',String(x===b)));drawPareto(b.dataset.salt)}));drawPareto('LiPF6');
const handles=document.getElementById('handle-count'),grafts=document.getElementById('graft-count'),electrons=document.getElementById('electron-count');
function drawCapacity(){
 const h=Number(handles.value),g=Number(grafts.value),e=electrons.value;
 const record=reportData.capacity.find(r=>r.reactive_handles===h&&r.grafts===g)[e+'e'];
 document.getElementById('capacity-value').textContent=record.capacity_upper_bound_mAh_g.toFixed(2);
 document.getElementById('capacity-composition').textContent=`${h} handles / 5 repeat units · ${g} grafts · ${h?Math.round(g/h*100)+'% handle conversion':'conversion not applicable'} · ${g*Number(e)} electrons assumed per segment.`;
 document.getElementById('capacity-mass').textContent=`Counted fragment mass: ${record.mass_g_mol.toFixed(2)} g/mol. A bound of zero to ${record.capacity_upper_bound_mAh_g.toFixed(2)} mAh/g is compatible with the stated accounting; usable charge is not established.`;
 let svg='<svg viewBox="0 0 640 180" role="img" aria-label="Five repeat units: '+g+' PMDI grafts, '+(h-g)+' unused chloride handles and '+(5-h)+' styrene units"><path d="M30 110H610" stroke="#51656b" stroke-width="4"/>';
 for(let i=0;i<5;i++){const x=80+i*120,kind=i<g?'PMDI':i<h?'CH₂Cl':'Styrene',color=i<g?'#146b5c':i<h?'#96551e':'#71837d';svg+=`<circle cx="${x}" cy="110" r="9" fill="${color}"/>`;if(i<h)svg+=`<path d="M${x} 100V58" stroke="${color}" stroke-width="3"/><rect x="${x-38}" y="22" width="76" height="36" rx="8" fill="${i<g?'#e9f3ee':'#fbefdc'}"/>`;svg+=`<text x="${x}" y="${i<h?46:150}" text-anchor="middle" fill="${color}" font-size="17">${kind}</text>`;}
 document.getElementById('capacity-chain').innerHTML=svg+'</svg>';
}
function resetGrafts(){const old=Number(grafts.value||3),h=Number(handles.value);grafts.innerHTML=Array.from({length:h+1},(_,i)=>`<option value="${i}" ${i===Math.min(old,h)?'selected':''}>${i} of ${h}</option>`).join('');drawCapacity();}
handles.addEventListener('change',resetGrafts);grafts.addEventListener('change',drawCapacity);electrons.addEventListener('change',drawCapacity);resetGrafts();
