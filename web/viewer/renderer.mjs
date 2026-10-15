// 表示データとカメラを受け取る。WebGL失敗時は同じ座標をCanvas 2Dで投影する。
export function cameraMatrix({yaw,pitch,zoom=1,panX=0,panY=0},width,height){
  const cy=Math.cos(yaw),sy=Math.sin(yaw),cp=Math.cos(pitch),sp=Math.sin(pitch),s=.8*zoom;
  const ax=s/Math.max(1,width/height),ay=s*Math.min(1,width/height);
  return new Float32Array([ax*cy,ay*sp*sy,cp*sy/4,0,0,ay*cp,-sp/4,0,ax*sy,-ay*sp*cy,-cp*cy/4,0,panX,panY,0,1]);
}
export function normalizeGeometry(g){
  const center=g.shown?g.lo.map((x,i)=>(x+g.hi[i])/2):[0,0,0];
  const radius=g.shown?Math.max(1,Math.hypot(...g.hi.map((x,i)=>(x-g.lo[i])/2))):1;
  const normalize=values=>Float32Array.from(values,(value,i)=>(value-center[i%3])/radius);
  return {positions:normalize(g.positions),colors:g.colors,centers:normalize(g.centers),centerColors:g.centerColors};
}
export function paintCanvas(ctx,g,matrix,width,height){
  ctx.clearRect(0,0,width,height);
  function point(a,i){const x=a[i],y=a[i+1],z=a[i+2];return [(matrix[0]*x+matrix[4]*y+matrix[8]*z+matrix[12]+1)*width/2,(1-matrix[1]*x-matrix[5]*y-matrix[9]*z-matrix[13])*height/2];}
  const css=(a,i)=>`rgb(${Math.round(a[i]*255)},${Math.round(a[i+1]*255)},${Math.round(a[i+2]*255)})`;
  // 同じ色の線をまとめ、線の本数だけstrokeを呼ばない。
  const paths=new Map();
  for(let i=0;i<g.positions.length;i+=6){const key=css(g.colors,i);if(!paths.has(key))paths.set(key,[]);paths.get(key).push(i);}
  ctx.lineWidth=Math.max(1,width/850);ctx.globalAlpha=.9;
  for(const [color,indices] of paths){ctx.strokeStyle=color;ctx.beginPath();for(const i of indices){const p=point(g.positions,i),q=point(g.positions,i+3);ctx.moveTo(...p);ctx.lineTo(...q);}ctx.stroke();}
  ctx.globalAlpha=1;const r=Math.max(1.5,width/380);
  for(let i=0;i<g.centers.length;i+=3){ctx.fillStyle=css(g.centerColors,i);ctx.beginPath();const p=point(g.centers,i);ctx.arc(p[0],p[1],r,0,Math.PI*2);ctx.fill();}
}
export function createRenderer(canvas,{preferCanvas=false}={}){
  let gl=null,program,position,color,matrix,size,buffers,geometry=null,context=null,active=canvas,mode='WebGL';
  function useCanvas(){
    if(context)return;
    // WebGLコンテキストを取得したcanvasは2Dに切り替えられないので別面を作る。
    active=canvas.ownerDocument.createElement('canvas');active.className='render-fallback';active.setAttribute('aria-hidden','true');
    canvas.parentNode.insertBefore(active,canvas);canvas.style.opacity='0';context=active.getContext('2d');
    if(!context)throw Error('WebGLとCanvasの両方で描画を開始できませんでした。');mode='Canvas';
  }
  try{
    if(preferCanvas)throw Error('Compatibility renderer selected');
    gl=canvas.getContext('webgl',{antialias:true,preserveDrawingBuffer:true});if(!gl)throw Error('WebGL unavailable');
    function shader(type,source){const s=gl.createShader(type);gl.shaderSource(s,source);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(s));return s;}
    program=gl.createProgram();gl.attachShader(program,shader(gl.VERTEX_SHADER,'attribute vec3 position;attribute vec3 color;uniform mat4 matrix;uniform float size;varying vec3 tint;void main(){gl_Position=matrix*vec4(position,1.0);gl_PointSize=size;tint=color;}'));
    gl.attachShader(program,shader(gl.FRAGMENT_SHADER,'precision mediump float;varying vec3 tint;void main(){gl_FragColor=vec4(tint,1.0);}'));gl.linkProgram(program);if(!gl.getProgramParameter(program,gl.LINK_STATUS))throw Error('WebGL link failed');
    gl.useProgram(program);position=gl.getAttribLocation(program,'position');color=gl.getAttribLocation(program,'color');matrix=gl.getUniformLocation(program,'matrix');size=gl.getUniformLocation(program,'size');
    buffers=Array.from({length:4},()=>gl.createBuffer());gl.enable(gl.DEPTH_TEST);gl.depthFunc(gl.LEQUAL);gl.clearColor(0,0,0,0);
  }catch{useCanvas();}
  function setGeometry(g){
    geometry=normalizeGeometry(g);
    if(mode==='WebGL')try{[geometry.positions,geometry.colors,geometry.centers,geometry.centerColors].forEach((a,i)=>{gl.bindBuffer(gl.ARRAY_BUFFER,buffers[i]);gl.bufferData(gl.ARRAY_BUFFER,a,gl.STATIC_DRAW);});if(gl.getError()!==gl.NO_ERROR)throw Error('GPU buffer error');}catch{useCanvas();}
  }
  function render(view){
    const dpr=Math.min(globalThis.devicePixelRatio||1,2),width=Math.max(1,Math.round(canvas.clientWidth*dpr)),height=Math.max(1,Math.round(canvas.clientHeight*dpr));
    if(!canvas.clientWidth||!canvas.clientHeight)throw Error('描画領域のサイズを取得できません。');
    if(active.width!==width||active.height!==height){active.width=width;active.height=height;}
    const m=cameraMatrix(view,width,height);
    if(mode==='WebGL')try{
      if(gl.isContextLost())throw Error('GPU context lost');gl.viewport(0,0,width,height);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);
      if(geometry){gl.uniformMatrix4fv(matrix,false,m);gl.uniform1f(size,Math.min(7,3*dpr));const bind=(index,attribute)=>{gl.bindBuffer(gl.ARRAY_BUFFER,buffers[index]);gl.enableVertexAttribArray(attribute);gl.vertexAttribPointer(attribute,3,gl.FLOAT,false,0,0);};bind(0,position);bind(1,color);gl.drawArrays(gl.LINES,0,geometry.positions.length/3);bind(2,position);bind(3,color);gl.drawArrays(gl.POINTS,0,geometry.centers.length/3);}
      if(gl.getError()!==gl.NO_ERROR)throw Error('GPU draw error');
    }catch{useCanvas();active.width=width;active.height=height;}
    if(mode==='Canvas'&&geometry)paintCanvas(context,geometry,m,width,height);
    return mode;
  }
  return {setGeometry,render,get canvas(){return active;},get mode(){return mode;}};
}
