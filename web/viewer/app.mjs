import {createRenderer} from './renderer.mjs';
import {loadFiles} from './loader.mjs';
import {geometry,compactModel} from './model.mjs';
const $=id=>document.getElementById(id), canvas=$('scene');
let model=null, busy=false, sourceFiles=[], cancelLoad=null;
function error(message){$('error').textContent=message;$('error').hidden=false;$('status').textContent=message;$('progressPanel').hidden=true;}
try{start();}catch(e){error(`描画初期化に失敗しました: ${e.message}`);}
function start(){
  const renderer=createRenderer(canvas,{preferCanvas:new URLSearchParams(globalThis.location?.search||'').get('renderer')!=='webgl'});
  let yaw=.72,pitch=-.48,zoom=1,panX=0,panY=0,scheduled=false;
  let summary='データ未読込み';
  function refresh(){
    if(!model)return;
    const g=geometry(model,{largestOnly:$('scope').value==='largest',showA:$('showA').checked,showB:$('showB').checked,arms:$('arms').checked,bonds:$('bonds').checked,maxMolecules:Number($('limit').value),clusterColors:$('color').value==='cluster'});
    renderer.setGeometry(g);
    $('displayCount').textContent=`描画 ${g.shown.toLocaleString()} / 対象 ${g.totalVisible.toLocaleString()} 分子${g.stride>1?' · 間引き表示':''}`;
    $('validation').textContent=model.validated?'詳細検証済み':'表示優先・詳細未検証';
    $('validate').disabled=busy||model.validated;
    summary=`腕長 ${model.armLength} · ${model.hasBonds?`${model.clusters.toLocaleString()}成分 · 同率最大 ${model.tied}個（最小分子IDの成分を選択）`:'結合CSVなし：クラスタは未判定'}${g.stride>1?' · 表示は一部の分子と、その両端が表示される結合のみ':''}`;
    draw();
  }
  function draw(){
    if(scheduled)return;scheduled=true;
    requestAnimationFrame(()=>{scheduled=false;try{
      const mode=renderer.render({yaw,pitch,zoom,panX,panY});
      if(model){$('progressPanel').hidden=true;$('status').textContent=`表示成功 · ${mode} · ${summary}`;}
    }catch(e){error(`描画に失敗しました: ${e.message}`);}});
  }
  function reset(){yaw=.72;pitch=-.48;zoom=1;panX=panY=0;draw();}
  function load(files,detailed=false){
    if(!files.length)return;cancelLoad?.();busy=true;$('sample').disabled=true;$('validate').disabled=true;$('error').hidden=true;$('progressPanel').hidden=false;
    cancelLoad=loadFiles(files,{detailed,onProgress:(percent,message)=>{$('progress').value=percent;$('progressLabel').textContent=`${message} · ${percent}%`;$('status').textContent=message;},onError:message=>{busy=false;$('sample').disabled=false;$('validate').disabled=!model||model.validated;$('progressPanel').hidden=true;error(message);},onResult:data=>{
      if(!data||!Array.isArray(data.molecules)||!Array.isArray(data.bonds))throw Error('描画データの形式が不正です。');
      busy=false;$('sample').disabled=false;model=data;sourceFiles=files;
      $('statN').textContent=model.molecules.length.toLocaleString();$('statSites').textContent=model.occupied.toLocaleString();$('statBonds').textContent=model.hasBonds?model.bonds.length.toLocaleString():'—';$('statLargest').textContent=model.hasBonds?model.largest.length.toLocaleString():'—';
      $('scope').options[1].disabled=!model.hasBonds;$('color').options[1].disabled=!model.hasBonds;
      if(!detailed){$('scope').value='all';$('color').value='type';}
      $('bonds').disabled=!model.hasBonds;$('fileLabel').textContent=files.map(f=>f.name).join(' + ');refresh();if(!detailed)reset();
    }});
  }
  $('cancel').onclick=()=>{cancelLoad?.();busy=false;$('sample').disabled=false;$('validate').disabled=!model||model.validated;$('progressPanel').hidden=true;$('status').textContent='読込みを中止しました。';};
  $('validate').onclick=()=>{if(sourceFiles.length&&!busy)load(sourceFiles,true);};
  $('export').onclick=()=>{if(!model)return;const blob=new Blob([JSON.stringify(compactModel(model))],{type:'application/json'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='model.ab-saw.json';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000);};
  $('files').onchange=e=>load([...e.target.files]);const drop=$('drop');for(const event of ['dragover','dragleave','drop'])drop.addEventListener(event,e=>{e.preventDefault();drop.classList.toggle('drag',event==='dragover');if(event==='drop')load([...e.dataTransfer.files].filter(f=>(/\.(csv|json)$/i).test(f.name)));});
  $('sample').onclick=async()=>{if(busy)return;try{$('sample').disabled=true;const files=await Promise.all(['trajectories.csv','bonds.csv'].map(async name=>{const r=await fetch(`sample/${name}`);if(!r.ok)throw Error('サンプルを読み込めませんでした。');return new File([await r.text()],name,{type:'text/csv'});}));load(files);}catch(e){$('sample').disabled=false;error(e.message);}};
  for(const id of ['scope','color','limit','showA','showB','arms','bonds'])$(id).onchange=refresh;
  let dragging=null;canvas.addEventListener('pointerdown',e=>{dragging={x:e.clientX,y:e.clientY};canvas.setPointerCapture(e.pointerId);canvas.focus();});canvas.addEventListener('pointermove',e=>{if(!dragging)return;const dx=e.clientX-dragging.x,dy=e.clientY-dragging.y;dragging={x:e.clientX,y:e.clientY};if(e.shiftKey){panX+=2*dx/canvas.clientWidth;panY-=2*dy/canvas.clientHeight;}else{yaw+=dx*.008;pitch+=dy*.008;}draw();});for(const name of ['pointerup','pointercancel','lostpointercapture'])canvas.addEventListener(name,()=>{dragging=null;});canvas.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.max(.15,Math.min(20,zoom*Math.exp(-e.deltaY*.001)));draw();},{passive:false});
  canvas.addEventListener('keydown',e=>{const actions={ArrowLeft:()=>yaw-=.1,ArrowRight:()=>yaw+=.1,ArrowUp:()=>pitch-=.1,ArrowDown:()=>pitch+=.1,'+':()=>zoom=Math.min(20,zoom*1.1),'=':()=>zoom=Math.min(20,zoom*1.1),'-':()=>zoom=Math.max(.15,zoom/1.1),r:reset,R:reset};if(actions[e.key]){e.preventDefault();actions[e.key]();draw();}});
  $('reset').onclick=reset;for(const b of document.querySelectorAll('[data-view]'))b.onclick=()=>{const views={iso:[.72,-.48],xy:[0,0],xz:[0,Math.PI/2],yz:[Math.PI/2,0]};[yaw,pitch]=views[b.dataset.view];zoom=1;panX=panY=0;draw();};
  $('save').onclick=()=>{if(!model)return;const out=document.createElement('canvas');out.width=renderer.canvas.width;out.height=renderer.canvas.height;const ctx=out.getContext('2d');ctx.fillStyle='#0b1321';ctx.fillRect(0,0,out.width,out.height);ctx.drawImage(renderer.canvas,0,0);ctx.fillStyle='#c9d6e4';ctx.font=`${Math.max(12,out.width/90)}px sans-serif`;ctx.fillText(`AB-SAW | ${$('displayCount').textContent} | ${model.validated?'validated':'unverified'}`,18,out.height-18);out.toBlob(blob=>{if(!blob)return;const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='ab-saw-view.png';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000);});};
  canvas.addEventListener('webglcontextlost',e=>{e.preventDefault();draw();});new ResizeObserver(draw).observe(canvas);$('sample').click();
}
