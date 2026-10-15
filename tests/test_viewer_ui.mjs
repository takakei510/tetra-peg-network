// ブラウザを起動しないUI契約テスト。実GPUの動作確認とは区別する。
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
const base=new URL('../web/viewer/',import.meta.url),html=readFileSync(new URL('index.html',base),'utf8');
const elements=new Map(),draws=[],errors=[];
for(const m of html.matchAll(/<([a-z]+)\b([^>]*\bid="([^"]+)"[^>]*)>/g))elements.set(m[3],{textContent:'',value:'',checked:/\bchecked\b/.test(m[2]),hidden:/\bhidden\b/.test(m[2]),disabled:false,options:[{},{}],style:{},classList:{toggle(){}},addEventListener(){},click(){return this.onclick?.();}});
// 固定の疑似optionsではなく、HTMLに実在する開始タグから選択肢を作る。
for(const select of html.matchAll(/<select\b[^>]*id="([^"]+)"[^>]*>([\s\S]*?)<\/select>/g)){
  elements.get(select[1]).options=[...select[2].matchAll(/<option\b([^>]*)>/g)].map(m=>({value:m[1].match(/value="([^"]+)"/)?.[1],disabled:/\bdisabled\b/.test(m[1])}));
}
for(const id of ['scope','color','limit'])elements.get(id).value={scope:'all',color:'type',limit:'5000'}[id];
const context={clearRect(){},beginPath(){},stroke(){},moveTo(){},lineTo(){},arc(){},fill(){draws.push('center');}};
const fakeCanvas=()=>({width:800,height:600,style:{},setAttribute(){},getContext:()=>context});
globalThis.location={search:process.argv.includes('--canvas')?'':'?renderer=webgl'};
globalThis.document={getElementById:id=>elements.get(id)||null,querySelectorAll:()=>[],createElement:fakeCanvas};
const gl=new Proxy({NO_ERROR:0,getError:()=>0,isContextLost:()=>false,getShaderParameter:()=>true,getProgramParameter:()=>true,getAttribLocation:()=>1,getUniformLocation:()=>1,createShader:()=>({}),createProgram:()=>({}),createBuffer:()=>({}),drawArrays:(...a)=>draws.push(a)}, {get:(t,k)=>t[k]??(String(k).toUpperCase()===k?1:()=>{})});
const canvas=elements.get('scene');Object.assign(canvas,{getContext:()=>process.argv.includes('--canvas')?null:gl,clientWidth:800,clientHeight:600,ownerDocument:document,parentNode:{insertBefore(){}}});
globalThis.devicePixelRatio=1;globalThis.requestAnimationFrame=cb=>setTimeout(cb,0);globalThis.ResizeObserver=class{observe(){}};
globalThis.fetch=async url=>({ok:true,text:async()=>readFileSync(new URL(url,base),'utf8')});
globalThis.Worker=class{constructor(){throw Error('Worker blocked in this test');}};
const capture=e=>errors.push(e.stack);process.on('uncaughtException',capture);process.on('unhandledRejection',capture);
await import(new URL('app.mjs',base));await new Promise(r=>setTimeout(r,60));
process.off('uncaughtException',capture);process.off('unhandledRejection',capture);
assert.deepEqual(elements.get('scope').options.map(o=>o.value),['all','largest']);
assert.deepEqual(elements.get('color').options.map(o=>o.value),['type','cluster']);
assert.deepEqual(errors,[]);assert.equal(elements.get('error').textContent,'');
assert.equal(elements.get('statN').textContent,'20');assert.match(elements.get('fileLabel').textContent,/trajectories/);
assert.match(elements.get('displayCount').textContent,/描画 20/);assert.match(elements.get('status').textContent,/表示成功/);assert(draws.length>0);assert.equal(elements.get('progressPanel').hidden,true);
assert.equal(elements.get('scope').options[1].disabled,false);
assert.equal(elements.get('color').options[1].disabled,false);
elements.get('scope').value='largest';elements.get('scope').onchange();
await new Promise(r=>setTimeout(r,10));
assert.match(elements.get('displayCount').textContent,/描画 2 \/ 対象 2/);
elements.get('color').value='cluster';elements.get('color').onchange();
await new Promise(r=>setTimeout(r,10));
assert.match(elements.get('status').textContent,/表示成功/);
console.log(`UI: sample → direct read → geometry → ${process.argv.includes('--canvas')?'Canvas':'WebGL calls'} → display success: passed`);
