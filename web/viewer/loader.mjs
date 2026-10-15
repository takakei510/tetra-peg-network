import {readFiles} from './io.mjs';
// Worker起動失敗とデータ処理のエラーを区別。小さいデータは同じパーサで直接読める。
export function loadFiles(files,{detailed=false,onProgress=()=>{},onResult=()=>{},onError=()=>{},workerFactory=()=>new Worker(new URL('./worker.mjs?v=2',import.meta.url),{type:'module'}),directRead=readFiles,startupMs=2500,timeoutMs=60000}={}) {
  let worker=null,done=false,started=false,fallback=false,startup,total;
  const cleanup=()=>{clearTimeout(startup);clearTimeout(total);worker?.terminate();worker=null;};
  const fail=e=>{if(done)return;done=true;cleanup();onError(e.message||String(e));};
  const success=model=>{if(done)return;cleanup();try{onProgress(90,'描画データを受渡し中');onResult(model);done=true;}catch(e){done=true;onError(`表示処理でエラー: ${e.message||String(e)}`);}};
  const direct=async()=>{if(done||fallback)return;fallback=true;cleanup();if(files.reduce((n,f)=>n+f.size,0)>2*1024*1024){fail(Error('Workerを起動できませんでした。再読込みするか、READMEのローカル起動をお試しください。直接読込みは2 MiBまでです。'));return;}onProgress(5,'直接読込みに切替中');total=setTimeout(()=>fail(Error('読込みが60秒を超えたため停止しました。')),timeoutMs);await new Promise(resolve=>setTimeout(resolve,0));if(done)return;try{success(await directRead(files,{detailed,onProgress:(p,m)=>{if(!done)onProgress(p,m);}}));}catch(e){fail(e);}};
  try{worker=workerFactory();worker.onmessage=({data})=>{if(done||fallback)return;if(data.type==='ready'){started=true;clearTimeout(startup);onProgress(5,'読込み処理を開始');worker.postMessage({files,detailed});}else if(data.type==='progress'){onProgress(data.percent,data.message);}else if(data.type==='result'){success(data.model);}else if(data.type==='error'){fail(Error(data.error));}};worker.onerror=e=>{e.preventDefault?.();if(!started)direct();else fail(Error('読込み処理が中断されました。ブラウザの空きメモリを確認してください。'));};startup=setTimeout(direct,startupMs);total=setTimeout(()=>fail(Error('読込みが60秒を超えたため停止しました。')),timeoutMs);}catch{direct();}
  return ()=>{done=true;cleanup();};
}
