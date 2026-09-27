"use strict";

let resultCharts=[];
window.addEventListener("resize",()=>{for(const chart of resultCharts)if(chart.canvas.isConnected)chart.redraw();});

// Results are requested on demand; polling the run never re-reads large artifacts.
async function loadResults(force=false){
  const id=state.selected;
  if(!id||(!force&&(state.resultsId===id||state.resultsLoading===id)))return;
  const request=(state.resultsRequest||0)+1;
  state.resultsRequest=request;
  state.resultsLoading=id;
  const root=$("results-content");
  root.replaceChildren(notice("正在校验并读取结果…"));
  $("refresh-results").disabled=true;
  try{
    const data=await api("/api/results?id="+encodeURIComponent(id));
    if(state.selected!==id||state.resultsRequest!==request)return;
    state.resultsId=id;
    renderResults(data);
  }catch(error){
    if(state.selected===id&&state.resultsRequest===request){root.replaceChildren(notice("结果暂时无法读取："+error.message,"error"));presentError(error);}
  }finally{
    if(state.resultsRequest===request){state.resultsLoading=null;$("refresh-results").disabled=false;}
  }
}
$("refresh-results").onclick=()=>loadResults(true);

function resultCard(title,description){
  const card=el("section",null,"result-card");
  card.append(el("h3",title),el("p",description,"hint"));
  return card;
}
function resultSource(card,uri){if(uri)card.append(el("p","产物来源 · "+uri,"result-source"));}
const resultCategories=[
  {key:"quality",title:"数据质量与实验设计",description:"输入数据与 QC 检查。",capabilities:["dataset_audit","qc"]},
  {key:"composition",title:"细胞组成与聚类",description:"展示已审计的细胞嵌入预览；差异丰度等尚无专属图表时，任务状态和证据仍可在这里、报告及证据包中查看。",capabilities:["clustering","subset_cells","differential_abundance"]},
  {key:"deg",title:"基因表达比较",description:"展示差异表达预览和已审计完整结果表的分页查询。",capabilities:["deg"]},
  {key:"trajectory",title:"轨迹与其他分析任务",description:"当前没有专属轨迹或其他扩展任务图表；请查看保存报告及完整证据包中的审计记录。",capabilities:["trajectory_inference","fate_prediction","cellrank_fate","donor_sensitivity","functional_activity","doublet_detection","cell_annotation","background_removal","cell_cell_communication","knowledge_retrieval","spatial_domain","spatial_deg"]}
];
function resultGroup(root,category){const section=el("section",null,"research-question-group");section.append(el("h2",category.title),el("p",category.description,"hint"));const rows=(state.resultTasks||[]).filter(task=>category.capabilities.includes(task.capability));if(rows.length){const list=el("div",null,"task-coverage");for(const task of rows){const row=el("div",null,"task-coverage-row"),description=el("div");description.append(el("strong",capabilities[task.capability]||task.capability),el("small",[task.task_id,task.method].filter(Boolean).join(" · ")||"快照未记录方法"));row.append(description,badge(task.status));if(task.audited)row.append(el("small","通过快照回执校验"));else row.append(el("small","产物未作为已审计结果展示","task-not-audited"));list.append(row);}section.append(list);}else section.append(el("p","快照未记录此类任务。","task-empty"));const cards=el("div",null,"results-grid");section.append(cards);root.append(section);return cards;}
function numberText(value){return typeof value==="number"&&Number.isFinite(value)?value.toLocaleString(undefined,{maximumFractionDigits:3}):"—";}
function renderResults(data){
  const root=$("results-content");root.replaceChildren();resultCharts=[];state.resultTasks=data.tasks||[];state.resultWarnings=data.warnings||[];
  if(data.notice)root.append(notice(data.notice));
  for(const warning of data.warnings||[])root.append(notice(typeof warning==="string"?warning:JSON.stringify(warning),"warning"));
  const groups=Object.fromEntries(resultCategories.map(category=>[category.key,resultGroup(root,category)]));
  if(data.qc)groups.quality.append(renderQc(data.qc));
  if(data.embedding)groups.composition.append(renderEmbedding(data.embedding));
  if(data.deg?.length)groups.deg.append(renderDeg(data.deg));
  groups.deg.append(renderDegTableQuery());
  if(!data.qc&&!data.embedding&&!data.deg?.length)root.prepend(el("div","当前快照没有可展示的 QC、聚类或 DEG 预览。未通过校验的产物不会用于绘图；仍可查看下方任务状态，并下载通过校验的完整证据包。","results-empty"));
  const reportLink=el("button","查看保存的分析报告","button secondary");reportLink.onclick=()=>selectDetailTab("report");root.append(reportLink);
  requestAnimationFrame(()=>{for(const chart of resultCharts)if(chart.canvas.isConnected)chart.redraw();});
}
function renderQc(qc){
  const card=resultCard("质量控制","查看本次过滤前后的细胞数量。");
  const summary=el("div",null,"result-summary");
  for(const [label,value] of [["输入细胞",qc.initial_cells],["保留细胞",qc.retained_cells],["过滤细胞",qc.filtered_cells]]){
    const item=el("div");item.append(el("strong",numberText(value)),el("small",label));summary.append(item);
  }
  card.append(summary);
  if(Number.isFinite(qc.initial_cells)&&qc.initial_cells>0&&Number.isFinite(qc.retained_cells)){
    const percent=Math.max(0,Math.min(100,qc.retained_cells/qc.initial_cells*100));
    const svg=document.createElementNS("http://www.w3.org/2000/svg","svg");
    svg.setAttribute("viewBox","0 0 100 4");svg.setAttribute("role","img");svg.setAttribute("aria-label","细胞保留率 "+percent.toFixed(1)+"%");
    for(const [width,color] of [[100,"#e8eeea"],[percent,"#368068"]]){const rect=document.createElementNS(svg.namespaceURI,"rect");rect.setAttribute("width",width);rect.setAttribute("height","4");rect.setAttribute("fill",color);svg.append(rect);}
    card.append(svg,el("p","保留率 "+percent.toFixed(1)+"%","hint"));
  }
  resultSource(card,qc.source_uri);return card;
}
const chartColors=["#267a63","#526cb3","#be783a","#976ba0","#368fa6","#a4555a","#818c39","#775d47","#5e8693","#a5758e","#49684b","#977d2d"];
function chartArea(card,label){
  const shell=el("div",null,"chart-shell"),canvas=el("canvas",null,"result-chart");
  canvas.width=1000;canvas.height=610;canvas.setAttribute("role","img");canvas.setAttribute("aria-label",label);
  const readout=el("div","悬停查看点位信息。","chart-readout");
  shell.append(canvas);card.append(shell,readout);
  return {canvas,readout};
}
function addPngExport(card,canvas,label){const button=el("button","下载带范围说明的 PNG","button ghost");button.onclick=()=>{const output=document.createElement("canvas");output.width=canvas.width;const first=output.getContext("2d");first.font="22px Segoe UI, sans-serif";const maxWidth=output.width-48,lines=[];for(const value of label()){let line="";for(const char of String(value)){if(line&&first.measureText(line+char).width>maxWidth){lines.push(line);line=char;}else line+=char;}if(line)lines.push(line);}const lineHeight=32;output.height=canvas.height+32+lines.length*lineHeight+16;const context=output.getContext("2d");context.fillStyle="#ffffff";context.fillRect(0,0,output.width,output.height);context.drawImage(canvas,0,0);context.fillStyle="#324b3d";context.font="22px Segoe UI, sans-serif";context.textAlign="left";lines.forEach((line,index)=>context.fillText(line,24,canvas.height+30+index*lineHeight));output.toBlob(blob=>{if(blob)downloadBlob(blob,(state.selected||"run").replaceAll("/","-")+"-chart.png","image/png");else presentError(Error("浏览器无法生成 PNG；请缩小页面窗口后重新绘图并下载。"));},"image/png");};card.append(button);}
function simulationMarker(){return (state.resultWarnings||[]).find(value=>/演示|模拟/.test(typeof value==="string"?value:JSON.stringify(value)));}
function drawScatter(canvas,points,{xLabel,yLabel,color,reference=null}){
  const displayWidth=canvas.getBoundingClientRect().width||500;
  canvas.width=Math.round(displayWidth*2);canvas.height=Math.round(Math.max(230,displayWidth/1.65)*2);
  const ctx=canvas.getContext("2d"),w=canvas.width,h=canvas.height;
  const left=88,right=28,top=28,bottom=80,plotW=w-left-right,plotH=h-top-bottom;
  const xs=points.map(p=>p.x),ys=points.map(p=>p.y);
  let xmin=points.length?Math.min(...xs):0,xmax=points.length?Math.max(...xs):1;
  let ymin=points.length?Math.min(...ys):0,ymax=points.length?Math.max(...ys):1;
  if(reference!==null){ymin=Math.min(ymin,reference);ymax=Math.max(ymax,reference);}
  const dx=(xmax-xmin)||1,dy=(ymax-ymin)||1;xmin-=dx*.06;xmax+=dx*.06;ymin-=dy*.06;ymax+=dy*.06;
  const px=x=>left+(x-xmin)/(xmax-xmin)*plotW,py=y=>top+plotH-(y-ymin)/(ymax-ymin)*plotH;
  ctx.clearRect(0,0,w,h);ctx.fillStyle="#f8faf7";ctx.fillRect(0,0,w,h);ctx.font="20px Segoe UI, sans-serif";
  ctx.lineWidth=1;
  const ticks=displayWidth<400?3:4;
  for(let i=0;i<=ticks;i++){
    const x=xmin+(xmax-xmin)*i/ticks,y=ymin+(ymax-ymin)*i/ticks;
    ctx.strokeStyle="#e0e8e1";ctx.beginPath();ctx.moveTo(px(x),top);ctx.lineTo(px(x),top+plotH);ctx.moveTo(left,py(y));ctx.lineTo(left+plotW,py(y));ctx.stroke();
    ctx.fillStyle="#60756a";ctx.textAlign="center";ctx.fillText(x.toFixed(1),px(x),h-bottom+32);ctx.textAlign="right";ctx.fillText(y.toFixed(1),left-13,py(y)+7);
  }
  if(reference!==null){ctx.strokeStyle="#97aa9c";ctx.setLineDash([7,6]);ctx.beginPath();ctx.moveTo(left,py(reference));ctx.lineTo(left+plotW,py(reference));ctx.stroke();ctx.setLineDash([]);}
  const hit=[];
  for(const point of points){const x=px(point.x),y=py(point.y);ctx.beginPath();ctx.fillStyle=color(point);ctx.globalAlpha=.78;ctx.arc(x,y,points.length>1500?3.3:5,0,Math.PI*2);ctx.fill();hit.push({x,y,point});}
  ctx.globalAlpha=1;ctx.fillStyle="#405c4e";ctx.textAlign="center";ctx.font="22px Segoe UI, sans-serif";ctx.fillText(xLabel,left+plotW/2,h-15);ctx.save();ctx.translate(24,top+plotH/2);ctx.rotate(-Math.PI/2);ctx.fillText(yLabel,0,0);ctx.restore();
  if(!points.length){ctx.fillText("没有匹配的数据点",left+plotW/2,top+plotH/2);}
  return hit;
}
function bindReadout(canvas,readout,getPoints,describe){
  canvas.onpointermove=event=>{
    const box=canvas.getBoundingClientRect();const x=(event.clientX-box.left)*canvas.width/box.width,y=(event.clientY-box.top)*canvas.height/box.height;
    let nearest=null,distance=400;
    for(const item of getPoints()){const d=(x-item.x)**2+(y-item.y)**2;if(d<distance){nearest=item;distance=d;}}
    readout.textContent=nearest?describe(nearest.point):"悬停查看点位信息。";
  };
  canvas.onpointerleave=()=>{readout.textContent="悬停查看点位信息。";};
}
function renderEmbedding(data){
  const card=resultCard(data.label||"细胞二维分布","展示 "+numberText(data.shown_cells)+" / "+numberText(data.total_cells)+" 个细胞"+(data.sampled?" · 已抽样；图例计数对应展示样本":"")+"。");
  const toolbar=el("div",null,"chart-toolbar"),label=el("label","颜色分组"),select=el("select");select.id="embedding-group";label.htmlFor=select.id;
  for(const [value,text] of [["cluster","聚类"],["cell_type","细胞类型"]]){const option=el("option",text);option.value=value;select.append(option);}
  toolbar.append(label,select);card.append(toolbar);
  const {canvas,readout}=chartArea(card,(data.label||"二维分布")+"，下方图例提供分组计数并可筛选。");
  const legend=el("div",null,"chart-legend");card.append(legend);
  let hit=[],hidden=new Set();
  const group=p=>String(p[select.value]??"未标注");
  const valid=data.points.filter(p=>Number.isFinite(p.x)&&Number.isFinite(p.y));
  let groups=[],colors=new Map();
  function redraw(){
    hit=drawScatter(canvas,valid.filter(p=>!hidden.has(group(p))),{xLabel:data.method?.includes("umap")?"UMAP 1":"维度 1",yLabel:data.method?.includes("umap")?"UMAP 2":"维度 2",color:p=>colors.get(group(p))});
  }
  function regroup(){
    hidden=new Set();groups=[...new Set(valid.map(group))];colors=new Map(groups.map((g,i)=>[g,chartColors[i%chartColors.length]]));legend.replaceChildren();
    const counts=new Map();for(const point of valid)counts.set(group(point),(counts.get(group(point))||0)+1);
    for(const g of groups){
      const button=el("button",g+" · "+counts.get(g),"legend-button");button.setAttribute("aria-pressed","true");button.title="显示或隐藏分组："+g;
      const dot=document.createElementNS("http://www.w3.org/2000/svg","svg");dot.setAttribute("width","10");dot.setAttribute("height","10");dot.setAttribute("aria-hidden","true");const circle=document.createElementNS(dot.namespaceURI,"circle");circle.setAttribute("cx","5");circle.setAttribute("cy","5");circle.setAttribute("r","4");circle.setAttribute("fill",colors.get(g));dot.append(circle);button.prepend(dot);
      button.onclick=()=>{if(hidden.has(g))hidden.delete(g);else hidden.add(g);button.setAttribute("aria-pressed",String(!hidden.has(g)));redraw();};legend.append(button);
    }
    redraw();
  }
  select.onchange=regroup;regroup();resultCharts.push({canvas,redraw});
  bindReadout(canvas,readout,()=>hit,p=>"聚类 "+(p.cluster??"未标注")+" · "+(p.cell_type??"未标注")+" · ("+p.x.toFixed(2)+", "+p.y.toFixed(2)+")");
  addPngExport(card,canvas,()=>[`来源 run=${state.selected} · task=${data.task_id||"见产物 URI"} · 方法=${data.method||"未记录"}`,`范围：${data.shown_cells} / ${data.total_cells} 个细胞${data.sampled?"（等距抽样）":"（全部有限坐标）"}；省略非有限坐标 ${data.omitted_nonfinite||0} 个`,`当前颜色分组=${select.value}；隐藏分组=${hidden.size}${simulationMarker()?"；含模拟 / 合成标记":""}`,...(simulationMarker()?[simulationMarker()]:[])]);
  card.append(el("p","点击图例可显示或隐藏分组；二维距离用于探索，不代表经过验证的生物学关系。","hint"));
  if(data.omitted_nonfinite)card.append(el("p","已省略 "+data.omitted_nonfinite+" 个含非有限坐标的细胞。","hint"));resultSource(card,data.source_uri);return card;
}
function scientificNumber(value){if(typeof value!=="number"||!Number.isFinite(value))return "—";return value!==0&&Math.abs(value)<.001?value.toExponential(2):value.toLocaleString(undefined,{maximumFractionDigits:4});}
function renderDeg(datasets){
  const card=resultCard("差异表达","按 FDR 从小到大预览最多 200 个基因；搜索和筛选仅作用于当前预览，完整结果保留在产物目录。");card.classList.add("result-wide");
  const toolbar=el("div",null,"chart-toolbar"),contrast=el("select"),threshold=el("select"),search=el("input"),only=el("input");
  contrast.id="deg-contrast";contrast.setAttribute("aria-label","选择差异表达结果");
  for(const [i,d] of datasets.entries()){const option=el("option",(d.target_cell_type||"全部细胞")+" · "+(d.contrast_label||((d.condition_a||"A")+" / "+(d.condition_b||"B")))+" · "+d.task_id);option.value=String(i);contrast.append(option);}
  threshold.id="deg-threshold";threshold.setAttribute("aria-label","FDR 显示阈值");
  for(const v of [.01,.05,.1]){const option=el("option","FDR < "+v);option.value=String(v);threshold.append(option);}threshold.value="0.05";
  search.type="search";search.placeholder="搜索基因…";search.setAttribute("aria-label","搜索差异表达基因");
  only.type="checkbox";only.id="deg-only";const onlyLabel=el("label",null,"check-label");onlyLabel.append(only,document.createTextNode("仅显示低于阈值的基因"));
  toolbar.append(contrast,threshold,onlyLabel);card.append(toolbar,search);
  const caption=el("p",null,"hint");card.append(caption);
  const {canvas,readout}=chartArea(card,"差异表达火山图：横轴 log2 fold change，纵轴负 log10 FDR；下方提供基因表。");
  const previewDetails=el("details",null,"preview-table-details"),previewSummary=el("summary","展开前 100 条预览基因表"),tableWrap=el("div",null,"results-table-wrap"),table=el("table",null,"results-table");tableWrap.append(table);
  const tableNote=el("p",null,"hint"),source=el("p",null,"result-source");previewDetails.append(previewSummary,tableWrap,tableNote,source);card.append(previewDetails);
  let hit=[];
  function redraw(){
    const d=datasets[Number(contrast.value)],limit=Number(threshold.value),query=search.value.trim().toLocaleLowerCase();
    const rows=d.rows.filter(r=>(!query||String(r.gene).toLocaleLowerCase().includes(query))&&(!only.checked||(Number.isFinite(r.fdr_q_value)&&r.fdr_q_value<limit)));
    const valid=rows.filter(r=>Number.isFinite(r.log2_fold_change)&&Number.isFinite(r.fdr_q_value)&&r.fdr_q_value>=0&&r.fdr_q_value<=1);
    const points=valid.map(r=>({...r,x:r.log2_fold_change,y:-Math.log10(Math.max(r.fdr_q_value,1e-300))}));
    hit=drawScatter(canvas,points,{xLabel:"log2 fold change",yLabel:"−log10 FDR",reference:-Math.log10(limit),color:p=>p.fdr_q_value<limit?(p.x>=0?"#b86d46":"#337ba0"):"#9aa99f"});
    caption.textContent="方法："+d.method+" · 统计单位："+(d.statistical_unit||"见报告")+" · 返回 "+numberText(d.shown_genes)+" / "+numberText(d.total_genes)+" 个基因 · 当前匹配 "+rows.length+" 个。灰色为未低于阈值，橙 / 蓝色区分正 / 负变化。";
    caption.textContent+=" 分析 alpha："+(d.alpha??0.05)+"；显示筛选不会改变分析结论。"+(d.effect_definition?" 效应定义："+d.effect_definition:"");
    table.replaceChildren();const head=el("thead"),heading=el("tr");for(const text of ["基因","log2 fold change","FDR q 值","p 值"]){const th=el("th",text);th.scope="col";heading.append(th);}head.append(heading);table.append(head);
    const body=el("tbody");for(const r of rows.slice(0,100)){const tr=el("tr");for(const value of [r.gene,scientificNumber(r.log2_fold_change),scientificNumber(r.fdr_q_value),scientificNumber(r.p_value)])tr.append(el("td",value));body.append(tr);}table.append(body);
    tableNote.textContent=(rows.length?"表格展示当前匹配的前 "+Math.min(rows.length,100)+" 项。":"没有匹配的基因。")+"缺失或非法统计值不绘图；FDR 为 0 时按 10⁻³⁰⁰ 显示。完整结果保留在产物目录。";
    source.textContent="产物来源 · "+d.source_uri;
  }
  contrast.onchange=threshold.onchange=only.onchange=search.oninput=redraw;redraw();resultCharts.push({canvas,redraw});
  bindReadout(canvas,readout,()=>hit,p=>p.gene+" · log2FC "+scientificNumber(p.log2_fold_change)+" · FDR "+scientificNumber(p.fdr_q_value));
  addPngExport(card,canvas,()=>{const d=datasets[Number(contrast.value)],filters="当前筛选：搜索="+(search.value.trim()||"无")+"；FDR 阈值="+threshold.value+"；仅显著="+(only.checked?"是":"否");return ["来源 run="+state.selected+" · task="+d.task_id+" · 方法="+(d.method||"未记录"),"火山图范围：DEG 预览 "+d.shown_genes+" / "+d.total_genes+" 个基因；不代表完整表",filters+"；完整表查询请使用下方后端分页表",...(simulationMarker()?[simulationMarker()]:[])]});
  return card;
}
function renderDegTableQuery(){
  const card=resultCard("已审计完整 DEG 表","独立于最多 200 行的图表预览。目录与每次分页查询均检查当前快照、成功任务和唯一通过审计回执；未重算差异分析。单表上限 512 MB；每次请求的完整来源哈希读取预算为 2 GB。点“查询”后才扫描完整表。");card.classList.add("result-wide");
  const toolbar=el("div",null,"chart-toolbar deg-query-toolbar"),discover=el("button","发现已审计 DEG 表","button secondary"),task=el("select"),query=el("input"),size=el("select"),fdr=el("input"),sig=el("input"),fc=el("input"),run=el("button","查询第 1 页","button primary"),download=el("button","下载完整 CSV","button ghost");
  task.setAttribute("aria-label","完整 DEG 结果表");task.add(new Option("先发现已审计结果表",""));query.type="search";query.maxLength=200;query.placeholder="搜索完整 gene ID（最多 200 字符）";query.setAttribute("aria-label","完整结果表中的基因 ID 搜索");
  size.setAttribute("aria-label","每页行数");size.add(new Option("每页 25 行","25"));size.add(new Option("每页 50 行","50"));size.add(new Option("每页 100 行","100"));size.value="50";
  fdr.type="number";fdr.min="0";fdr.max="1";fdr.step="0.01";fdr.value="0.05";fdr.setAttribute("aria-label","FDR 最大值");sig.type="checkbox";const fdrLabel=el("label",null,"check-label deg-fdr-filter");fdrLabel.append(sig,document.createTextNode(" 仅 FDR ≤ "),fdr);
  fc.type="number";fc.min="0";fc.step="0.1";fc.placeholder="可选最小 |log2FC|";fc.setAttribute("aria-label","最小绝对 log2 fold change");const tableLabel=el("label","结果表");tableLabel.append(task);toolbar.append(discover,tableLabel,query,size,fdrLabel,fc,run,download);card.append(toolbar);
  const status=el("p","尚未发现完整 DEG 表。","hint"),warnings=el("div"),pageNav=el("div",null,"table-pagination"),previous=el("button","← 上一页","button ghost"),pageText=el("strong","—"),next=el("button","下一页 →","button ghost"),tableWrap=el("div",null,"results-table-wrap"),table=el("table",null,"results-table"),foot=el("p","","hint");tableWrap.append(table);pageNav.append(previous,pageText,next);card.append(status,warnings,tableWrap,pageNav,foot);run.disabled=true;download.disabled=true;previous.disabled=true;next.disabled=true;
  let tables=[],page=1;state.degQuerySequence=state.degQuerySequence||0;
  function drawRows(rows){table.replaceChildren();const head=el("thead"),heading=el("tr");for(const name of ["完整 gene ID","log2 fold change","FDR q 值","p 值"]){const th=el("th",name);th.scope="col";heading.append(th);}head.append(heading);table.append(head);const body=el("tbody");for(const row of rows){const tr=el("tr");for(const value of [row.gene,scientificNumber(row.log2_fold_change),scientificNumber(row.fdr_q_value),scientificNumber(row.p_value)])tr.append(el("td",value));body.append(tr);}table.append(body);}
  function setOptions(){task.replaceChildren(new Option(tables.length?"选择完整 DEG 表":"未发现可查询表",""));for(const item of tables){const label=[item.target_cell_type||"全部细胞",(item.contrast_label||((item.condition_a||"A")+" / "+(item.condition_b||"B"))),item.method||"方法未记录",item.task_id].join(" · ");task.add(new Option(label,item.task_id));}if(tables.length)task.value=tables[0].task_id;run.disabled=!task.value;download.disabled=!task.value;}
  async function discoverTables(){const id=state.selected;if(!id)return;const sequence=++state.degQuerySequence;discover.disabled=true;status.textContent="正在校验快照与已审计 DEG 表目录…";warnings.replaceChildren();try{const data=await api("/api/results/deg-tables?id="+encodeURIComponent(id));if(state.selected!==id||state.degQuerySequence!==sequence)return;tables=data.tables||[];setOptions();status.textContent="发现 "+tables.length+" 个成功且通过持久审计的完整 DEG 表；单表上限 "+numberText(data.max_table_bytes)+" 字节。 "+(data.notice||"");for(const warning of data.warnings||[])warnings.append(notice(warning,"warning"));if(!tables.length)status.textContent+=" 当前运行没有符合条件的可查询表。";page=1;drawRows([]);pageText.textContent="—";}catch(error){if(state.selected===id&&state.degQuerySequence===sequence){status.textContent="目录校验失败："+error.message;presentError(error);}}finally{discover.disabled=false;}}
  function tableFilter(){return {q:query.value.trim(),size:size.value,fdr:sig.checked?fdr.value:null,abs_log2fc_min:fc.value.trim()||null};}
  function invalidateTableQuery(){state.degQuerySequence++;page=1;drawRows([]);pageText.textContent="待重新查询";previous.disabled=true;next.disabled=true;run.disabled=!task.value;run.textContent="查询第 1 页";status.textContent="筛选条件已修改；点击“查询第 1 页”应用新条件。";}
  async function queryPage(targetPage=1){
    if(!task.value)throw Error("请先发现并选择一份完整 DEG 表");
    const id=state.selected,taskId=task.value,filters=tableFilter(),filtersKey=JSON.stringify(filters),sequence=++state.degQuerySequence;
    const params=new URLSearchParams({id:id,task_id:taskId,page:String(targetPage),size:filters.size,q:filters.q,significant_only:String(sig.checked)});
    if(filters.fdr!==null){params.set("fdr_max",filters.fdr);params.set("significant_only","true");}
    if(filters.abs_log2fc_min!==null)params.set("abs_log2fc_min",filters.abs_log2fc_min);
    run.disabled=true;status.textContent="正在校验运行快照并扫描 "+taskId+" 的完整 CSV 表…";
    try{
      const result=await api("/api/results/deg?"+params.toString());
      if(state.selected!==id||task.value!==taskId||state.degQuerySequence!==sequence||JSON.stringify(tableFilter())!==filtersKey)return;
      page=result.page;drawRows(result.rows||[]);
      pageText.textContent=result.total_items?"第 "+result.page+" / "+result.total_pages+" 页 · "+numberText(result.total_items)+" 条匹配结果":"没有匹配结果 · 0 条";
      previous.disabled=!result.total_items||page<=1;next.disabled=!result.total_items||page>=result.total_pages;
      status.textContent=[result.target_cell_type||"全部细胞",(result.contrast_label||((result.condition_a||"A")+" / "+(result.condition_b||"B"))),result.method,result.statistical_unit||"统计单位见报告"].join(" · ")+"。"+result.notice;
      foot.textContent="列：gene、log2_fold_change、fdr_q_value、p_value；显示完整 gene ID。当前筛选：搜索="+(filters.q||"无")+"，仅 FDR≤"+(filters.fdr===null?"关闭":filters.fdr)+"，最小 |log2FC|="+(filters.abs_log2fc_min||"关闭")+"。";
    }catch(error){if(state.selected===id&&state.degQuerySequence===sequence){status.textContent="结果表查询失败："+error.message;presentError(error);}}
    finally{if(state.selected===id&&task.value===taskId&&state.degQuerySequence===sequence)run.disabled=false;}
  }
  discover.onclick=handle(discoverTables);run.onclick=handle(()=>queryPage(1));previous.onclick=handle(()=>queryPage(Math.max(1,page-1)));next.onclick=handle(()=>queryPage(page+1));task.onchange=()=>{state.degQuerySequence++;drawRows([]);page=1;pageText.textContent="待查询";previous.disabled=true;next.disabled=true;status.textContent="已切换完整结果表；点击“查询”执行新条件。";run.disabled=!task.value;download.disabled=!task.value;};
  query.addEventListener("input",invalidateTableQuery);size.addEventListener("change",invalidateTableQuery);sig.addEventListener("change",invalidateTableQuery);fdr.addEventListener("input",invalidateTableQuery);fc.addEventListener("input",invalidateTableQuery);
  download.onclick=handle(()=>downloadEndpoint("/api/results/deg.csv?id="+encodeURIComponent(state.selected)+"&task_id="+encodeURIComponent(task.value),task.value+"-full-deg.csv"));
  return card;
}
