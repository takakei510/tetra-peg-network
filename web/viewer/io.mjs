import {parseModel,parseCompact,LIMITS} from './model.mjs';
export async function readFiles(files,{detailed=false,onProgress=()=>{}}={}) {
  if(!files.length)throw Error('ファイルを選択してください。');
  if(files.reduce((n,f)=>n+f.size,0)>LIMITS.bytes)throw Error('ファイル合計は32 MiBまでです。');
  const light=files.find(f=>f.name.endsWith('.ab-saw.json'));
  onProgress(10,'ファイルを読み込み中');
  if(light){
    const data=JSON.parse(await light.text());
    if(data.format!=='ab-saw-viewer-1'||!Array.isArray(data.molecules)||!Array.isArray(data.bonds))throw Error('表示専用データの形式が違います。');
    if(data.molecules.length>LIMITS.rows/5||data.bonds.length>LIMITS.bonds)throw Error('表示専用データの件数上限を超えています。');
    return parseCompact(data,{detailed,onProgress});
  }
  let t=files.find(f=>/trajectories\.csv$/i.test(f.name));if(!t&&files.length===1)t=files[0];if(!t)throw Error('trajectories.csvを選択してください。');
  const b=files.find(f=>/^bonds\.csv$/i.test(f.name));
  const trajectory=await t.text(),bonds=b?await b.text():null;
  return parseModel(trajectory,bonds,{detailed,onProgress});
}
