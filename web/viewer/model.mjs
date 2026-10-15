// 現行C出力の数値CSVを読み込む。物理モデル・保存CSVは変更しない。
export const LIMITS = {bytes: 32 * 1024 * 1024, rows: 500000, bonds: 200000};
function table(text, required, limit) {
  const lines = text.replace(/^\uFEFF/, '').trim().split(/\r?\n/);
  const fields = lines.shift().split(',').map(x => x.trim());
  if (required.some(x => !fields.includes(x))) throw Error(`必要な列: ${required.join(', ')}`);
  if (lines.length > limit) throw Error(`行数上限 ${limit.toLocaleString()} を超えています。小さい試料を選んでください。`);
  return lines.filter(x => x.trim()).map((line, i) => {
    const values = line.split(',').map(x => x.trim());
    if (values.length !== fields.length) throw Error(`${i + 2}行目の列数が不正です（生成プログラムのCSVを使用してください）。`);
    return Object.fromEntries(fields.map((f, k) => [f, values[k]]));
  });
}
function integer(v, name, min = 0) {
  if (v === undefined || !/^-?\d+$/.test(v)) throw Error(`${name}が整数ではありません。`);
  const n = Number(v);
  if (!Number.isSafeInteger(n) || n < min) throw Error(`${name}が範囲外です。`);
  return n;
}
export function parseModel(trajectory, bondText = null, {detailed = true, onProgress = () => {}} = {}) {
  onProgress(25, "軌跡CSVを解析中");
  const rows = table(trajectory, ['molecule','type','arm','step','x','y','z'], LIMITS.rows);
  if (!rows.length) throw Error('軌跡CSVが空です。');
  const mols = new Map(), occupied = new Set();
  for (const row of rows) {
    const id = integer(row.molecule, '分子ID'), arm = integer(row.arm, '腕ID', -1), step = integer(row.step, 'step');
    if (!['A','B'].includes(row.type) || arm > 3 || (arm === -1 ? step !== 0 : step < 1)) throw Error('分子型・腕ID・stepが不正です。');
    if (step > rows.length) throw Error('stepが軌跡行数を超えています。');
    const point = ['x','y','z'].map(k => integer(row[k], k));
    const key = point.join(',');
    if (detailed && occupied.has(key)) throw Error('軌跡に重複占有点があります。');
    if (detailed) occupied.add(key);
    if (!mols.has(id)) mols.set(id, {id, type: row.type, center:null, arms:[[],[],[],[]]});
    const m = mols.get(id);
    if (m.type !== row.type) throw Error('同じ分子のA/B型が一致しません。');
    if (arm === -1) { if (m.center) throw Error('分子中心が重複しています。'); m.center = point; }
    else { if (m.arms[arm][step-1]) throw Error('腕のstepが重複しています。'); m.arms[arm][step-1] = point; }
  }
  onProgress(50, "分子形状を構築中");
  let length = null;
  for (const m of mols.values()) {
    if (!m.center) throw Error('分子中心がありません。');
    for (const arm of m.arms) {
      if (!arm.length || arm.some(p => !p) || Object.keys(arm).length !== arm.length) throw Error('腕のstepが欠けています。');
      if (length === null) length = arm.length;
      if (length !== arm.length) throw Error('腕長が一致しません。');
      let previous = m.center;
      for (const p of arm) { if (detailed && p.reduce((s,v,k) => s + Math.abs(v-previous[k]),0) !== 1) throw Error('腕に最近接でない移動があります。'); previous = p; }
    }
  }
  onProgress(65, "結合と連結成分を集計中");
  const parent = new Map([...mols.keys()].map(id => [id,id]));
  const find = id => { let r=id; while(parent.get(r)!==r) r=parent.get(r); while(id!==r) {const n=parent.get(id);parent.set(id,r);id=n;} return r; };
  const bonds = [], used = new Set(), bondIds = new Set();
  if (bondText !== null) for (const r of table(bondText, ['bond_id','molecule_a','molecule_b','arm_a','arm_b','endpoint_a','endpoint_b'], LIMITS.bonds)) {
    const id=integer(r.bond_id,'結合ID'); if(bondIds.has(id)) throw Error('結合IDが重複しています。'); bondIds.add(id);
    const ends=['a','b'].map(side => {
      const m=integer(r[`molecule_${side}`],'分子ID'), a=integer(r[`arm_${side}`],'腕ID'), e=integer(r[`endpoint_${side}`],'末端ID');
      if(!mols.has(m)||a>3||e!==4*m+a||used.has(e)) throw Error('結合の分子・腕・末端IDが不正、または末端が重複しています。');
      used.add(e); const p=mols.get(m).arms[a].at(-1);
      for(const [k,axis] of ['x','y','z'].entries()) if(detailed && r[`${axis}_${side}`]!==undefined && integer(r[`${axis}_${side}`],axis)!==p[k]) throw Error('軌跡と結合の座標が一致しません。');
      return {m,p};
    });
    const [a,b]=ends, distance=a.p.reduce((s,v,k)=>s+Math.abs(v-b.p[k]),0);
    if(a.m===b.m||mols.get(a.m).type===mols.get(b.m).type||(detailed && (distance<1||distance>2))) throw Error('A–B結合または距離条件に違反しています。');
    const ra=find(a.m),rb=find(b.m); if(ra!==rb) parent.set(Math.max(ra,rb),Math.min(ra,rb));
    bonds.push({a:a.m,b:b.m,p:a.p,q:b.p});
  }
  onProgress(80, "表示データを準備中");
  const groups=new Map(); for(const id of mols.keys()) {const root=find(id);if(!groups.has(root))groups.set(root,[]);groups.get(root).push(id);}
  const sorted=[...groups.entries()].sort((a,b)=>b[1].length-a[1].length || a[0]-b[0]);
  const largest=sorted[0][1], tied=sorted.filter(g=>g[1].length===largest.length).length;
  const molecules=[...mols.values()].sort((a,b)=>a.id-b.id);
  return {validated:detailed,molecules,bonds,occupied:rows.length,armLength:length,hasBonds:bondText!==null,largest,tied,clusters:groups.size,groupById:Object.fromEntries([...mols.keys()].map(id=>[id,find(id)]))};
}
export function geometry(model, {largestOnly=false, showA=true, showB=true, arms=true, bonds=true, maxMolecules=5000, clusterColors=false}={}) {
  const selected=new Set(model.largest);
  const visible=model.molecules.filter(m=>(!largestOnly||selected.has(m.id))&&(m.type==='A'?showA:showB));
  const stride=Math.max(1,Math.ceil(visible.length/maxMolecules));
  const shown=visible.filter((_,i)=>i%stride===0), ids=new Set(shown.map(m=>m.id));
  const segments=[],colors=[],centers=[],centerColors=[],lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];
  const color=m=>{if(!clusterColors)return m.type==='A'?[.22,.55,.94]:[1,.43,.34];const h=(Number(model.groupById[m.id])*137.508)%360;const c=.68,x=c*(1-Math.abs((h/60)%2-1)),v=.22;return (h<60?[c,x,0]:h<120?[x,c,0]:h<180?[0,c,x]:h<240?[0,x,c]:h<300?[x,0,c]:[c,0,x]).map(z=>z+v);};
  function bounds(p){for(let k=0;k<3;k++){lo[k]=Math.min(lo[k],p[k]);hi[k]=Math.max(hi[k],p[k]);}}
  const segment=(p,q,c)=>{segments.push(...p,...q);colors.push(...c,...c);};
  for(const m of shown){const c=color(m);centers.push(...m.center);centerColors.push(...c);bounds(m.center);for(const arm of m.arms){let p=m.center;for(const q of arm){bounds(q);if(arms)segment(p,q,c);p=q;}}}
  let shownBonds=0;
  if(bonds)for(const b of model.bonds)if(ids.has(b.a)&&ids.has(b.b)){segment(b.p,b.q,[.34,.92,.65]);shownBonds++;}
  return {positions:new Float32Array(segments),colors:new Float32Array(colors),centers:new Float32Array(centers),centerColors:new Float32Array(centerColors),lo,hi,shown:shown.length,totalVisible:visible.length,shownBonds,stride};
}

export function compactModel(model) {
  const byId=new Map(model.molecules.map(m=>[m.id,m]));
  const arm=(id,p)=>byId.get(id).arms.findIndex(a=>a.at(-1).every((v,k)=>v===p[k]));
  return {format:'ab-saw-viewer-1',hasBonds:model.hasBonds,molecules:model.molecules.map(m=>[m.id,m.type,...m.center,m.arms.map(a=>a.flat())]),bonds:model.bonds.map(b=>[b.a,b.b,arm(b.a,b.p),arm(b.b,b.q)])};
}
export function parseCompact(data,{detailed=true,onProgress=()=>{}}={}) {
  onProgress(25,'表示専用データを展開中');
  const molecules=[],byId=new Map(),occupied=new Set();let count=0,length=null;
  function number(n){if(!Number.isSafeInteger(n)||n<0)throw Error('表示専用データの数値が不正です。');return n;}
  function point(p){p.forEach(number);count++;if(count>LIMITS.rows)throw Error('軌跡行数上限を超えています。');if(detailed){const key=p.join(',');if(occupied.has(key))throw Error('軌跡に重複占有点があります。');occupied.add(key);}return p;}
  if(data.format!=='ab-saw-viewer-1'||typeof data.hasBonds!=='boolean'||!Array.isArray(data.molecules)||!data.molecules.length||!Array.isArray(data.bonds)||data.bonds.length>LIMITS.bonds)throw Error('表示専用データの形式が不正です。');
  for(const row of data.molecules){
    if(!Array.isArray(row)||row.length!==6||!['A','B'].includes(row[1])||!Array.isArray(row[5])||row[5].length!==4)throw Error('分子データが不正です。');
    const id=number(row[0]);if(byId.has(id))throw Error('分子IDが重複しています。');
    const center=point(row.slice(2,5));const arms=row[5].map(flat=>{if(!Array.isArray(flat)||!flat.length||flat.length%3)throw Error('腕の座標が不正です。');if(length===null)length=flat.length/3;if(flat.length!==length*3)throw Error('腕長が一致しません。');let previous=center;const arm=[];for(let i=0;i<flat.length;i+=3){const p=point(flat.slice(i,i+3));if(detailed&&p.reduce((s,v,k)=>s+Math.abs(v-previous[k]),0)!==1)throw Error('腕に最近接でない移動があります。');arm.push(p);previous=p;}return arm;});
    const m={id,type:row[1],center,arms};molecules.push(m);byId.set(id,m);
  }
  onProgress(65,'結合と連結成分を集計中');
  const parent=new Map(molecules.map(m=>[m.id,m.id]));const find=id=>{let r=id;while(parent.get(r)!==r)r=parent.get(r);while(id!==r){const n=parent.get(id);parent.set(id,r);id=n;}return r;};
  const bonds=[],used=new Set();if(!data.hasBonds&&data.bonds.length)throw Error('結合の有無とデータが一致しません。');
  for(const row of data.bonds){if(!Array.isArray(row)||row.length!==4)throw Error('結合データが不正です。');row.forEach(number);const [a,b,aa,ab]=row,ma=byId.get(a),mb=byId.get(b);if(!ma||!mb||a===b||ma.type===mb.type||aa>3||ab>3)throw Error('結合先が不正です。');const ea=`${a}:${aa}`,eb=`${b}:${ab}`;if(used.has(ea)||used.has(eb))throw Error('末端が重複使用されています。');used.add(ea);used.add(eb);const p=ma.arms[aa].at(-1),q=mb.arms[ab].at(-1),distance=p.reduce((s,v,k)=>s+Math.abs(v-q[k]),0);if(detailed&&(distance<1||distance>2))throw Error('結合の距離条件に違反しています。');bonds.push({a,b,p,q});const ra=find(a),rb=find(b);if(ra!==rb)parent.set(Math.max(ra,rb),Math.min(ra,rb));}
  const groups=new Map();for(const m of molecules){const r=find(m.id);if(!groups.has(r))groups.set(r,[]);groups.get(r).push(m.id);}const sorted=[...groups.entries()].sort((a,b)=>b[1].length-a[1].length||a[0]-b[0]);const largest=sorted[0][1].sort((a,b)=>a-b);
  onProgress(80,'表示データを準備中');return {validated:detailed,molecules:molecules.sort((a,b)=>a.id-b.id),bonds,occupied:count,armLength:length,hasBonds:data.hasBonds,largest,tied:sorted.filter(g=>g[1].length===largest.length).length,clusters:groups.size,groupById:Object.fromEntries(molecules.map(m=>[m.id,find(m.id)]))};
}
