import {readFiles} from './io.mjs';
self.postMessage({type:'ready'});
self.onmessage=async ({data})=>{try{const model=await readFiles(data.files,{detailed:data.detailed,onProgress:(percent,message)=>self.postMessage({type:'progress',percent,message})});self.postMessage({type:'result',model});}catch(e){self.postMessage({type:'error',error:e.message});}};
