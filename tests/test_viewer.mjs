import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {parseModel,geometry} from '../web/viewer/model.mjs';
const path=new URL('../web/viewer/sample/',import.meta.url);
const t=readFileSync(new URL('trajectories.csv',path),'utf8'), b=readFileSync(new URL('bonds.csv',path),'utf8');
const model=parseModel(t,b);
assert.equal(model.molecules.length,20);assert.equal(model.occupied,420);assert.equal(model.armLength,5);assert.equal(model.bonds.length,4);assert.equal(model.clusters,16);assert.deepEqual(model.largest,[1,19]);assert.equal(model.tied,4);
const all=geometry(model);assert.equal(all.shown,20);assert.equal(all.positions.length,(400+4)*6);
const sparse=geometry(model,{maxMolecules:3});assert(sparse.shown<=3);assert.equal(sparse.totalVisible,20);
const largest=geometry(model,{largestOnly:true});assert.equal(largest.shown,model.largest.length);
assert.equal(geometry(model,{showA:false,showB:false}).shown,0);
assert.equal(parseModel(t).hasBonds,false);assert.equal(parseModel(t,'bond_id,molecule_a,molecule_b,arm_a,arm_b,endpoint_a,endpoint_b\n').largest.length,1);
assert.throws(()=>parseModel(t+'\n'+t.split('\n')[1]),/重複/);
assert.throws(()=>parseModel(t.replace(',A,',',C,')),/分子型/);
assert.throws(()=>parseModel(t,b+'\n'+b.split('\n')[1]),/重複/);
assert.throws(()=>parseModel('x,y,z\n1,2,3'),/必要な列/);
// 二重結合を保持し、連結成分は2分子の1個になる。
const header='molecule,type,arm,step,x,y,z\n';
const trajectory=header+'0,A,-1,0,2,2,2\n0,A,0,1,3,2,2\n0,A,1,1,1,2,2\n0,A,2,1,2,3,2\n0,A,3,1,2,2,3\n1,B,-1,0,4,4,2\n1,B,0,1,3,4,2\n1,B,1,1,4,3,2\n1,B,2,1,5,4,2\n1,B,3,1,4,4,3\n';
const bondHeader='bond_id,molecule_a,molecule_b,arm_a,arm_b,endpoint_a,endpoint_b\n';
const parallel=bondHeader+'0,0,1,0,1,0,5\n1,0,1,2,0,2,4\n';
const two=parseModel(trajectory,parallel);assert.equal(two.bonds.length,2);assert.equal(two.clusters,1);assert.deepEqual(two.largest,[0,1]);
const tied=parseModel(trajectory,bondHeader);assert.deepEqual(tied.largest,[0]);assert.equal(tied.tied,2);
assert.throws(()=>parseModel(trajectory,parallel.replace('1,0,1,2,0,2,4','1,0,1,0,0,0,4')),/末端/);
assert.throws(()=>parseModel(trajectory.replace('0,A,0,1','0,A,0,999999999999'),parallel),/step/);
console.log('Viewer parser, multibonds, geometry, sampling, ties, invalid data: passed');

// 表示専用形式も全体とクラスタを再計算。保存済みの集計値は信用しない。
const {compactModel,parseCompact}=await import('../web/viewer/model.mjs');
const compact=compactModel(model),roundtrip=parseCompact(compact);
assert.deepEqual(roundtrip,model);
assert.equal(parseModel(t,b,{detailed:false}).validated,false);
const corrupt=structuredClone(compact);corrupt.molecules[1][2]=corrupt.molecules[0][2];corrupt.molecules[1][3]=corrupt.molecules[0][3];corrupt.molecules[1][4]=corrupt.molecules[0][4];
assert.throws(()=>parseCompact(corrupt),/重複/);
assert.equal(parseCompact(compact,{detailed:false}).validated,false);
const {readFiles}=await import('../web/viewer/io.mjs');
const light=new File([JSON.stringify(compact)],'model.ab-saw.json');
assert.deepEqual(await readFiles([light],{detailed:true}),model);
const {loadFiles}=await import('../web/viewer/loader.mjs');
await new Promise((resolve,reject)=>loadFiles([light],{workerFactory:()=>{throw Error('Worker blocked');},onResult:m=>{assert.equal(m.molecules.length,20);resolve();},onError:reject}));
await new Promise((resolve,reject)=>loadFiles([light],{startupMs:5,workerFactory:()=>({terminate(){}}),onResult:()=>resolve(),onError:reject}));
await new Promise((resolve,reject)=>loadFiles([{size:3*1024*1024}],{workerFactory:()=>{throw Error('blocked');},onResult:()=>reject(Error('large fallback must be rejected')),onError:e=>{assert.match(e,/2 MiB/);resolve();}}));
await new Promise((resolve,reject)=>{const cancel=loadFiles([light],{startupMs:5,workerFactory:()=>({terminate(){}}),onResult:()=>reject(Error('cancel failed')),onError:reject});cancel();setTimeout(resolve,15);});
await new Promise((resolve,reject)=>loadFiles([light],{workerFactory:()=>{const w={terminate(){},postMessage(){queueMicrotask(()=>{w.onmessage({data:{type:'progress',percent:50,message:'test'}});w.onmessage({data:{type:'result',model}});});}};queueMicrotask(()=>w.onmessage({data:{type:'ready'}}));return w;},onResult:m=>{assert.equal(m,model);resolve();},onError:reject}));
console.log('Compact data, staged loading, Worker handshake/fallback/cancel: passed');

const {execFileSync}=await import('node:child_process');
const converted=JSON.parse(execFileSync('python3',['-c',"import json;from scripts.export_viewer import convert;print(json.dumps(convert('web/viewer/sample/trajectories.csv','web/viewer/sample/bonds.csv')))"],{encoding:'utf8'}));
assert.deepEqual(parseCompact(converted),model);
console.log('Python exporter and JavaScript importer agree');
// 読込み成功後の表示例外を握りつぶさない。
await new Promise((resolve,reject)=>loadFiles([light],{workerFactory:()=>{throw Error('blocked');},onResult:()=>{throw Error('render callback failed');},onError:e=>{assert.match(e,/表示処理でエラー: render callback failed/);resolve();}}));
const {cameraMatrix,normalizeGeometry}=await import('../web/viewer/renderer.mjs');
const projected=normalizeGeometry(geometry(model));assert(projected.positions.every(Number.isFinite));
for(const [yaw,pitch] of [[.72,-.48],[0,0],[0,Math.PI/2],[Math.PI/2,0]]){const matrix=cameraMatrix({yaw,pitch},800,600);assert(matrix.every(Number.isFinite));assert(projected.positions.length>0);for(let i=0;i<projected.positions.length;i+=3){const a=projected.positions;const x=matrix[0]*a[i]+matrix[4]*a[i+1]+matrix[8]*a[i+2],y=matrix[1]*a[i]+matrix[5]*a[i+1]+matrix[9]*a[i+2];assert(Math.abs(x)<1&&Math.abs(y)<1);}}
assert(normalizeGeometry(geometry(model,{showA:false,showB:false})).positions.every(Number.isFinite));
console.log('Display exceptions and finite visible camera projection: passed');
