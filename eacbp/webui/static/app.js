"use strict";
const $ = id => document.getElementById(id);
const token = document.querySelector('meta[name="eacbp-token"]').content;
const state = {view:"new", runs:[], jobs:[], selected:null, preview:null, previewKey:null, designConfirmed:false, busy:false, folder:null, detail:null, fileTarget:null, importPreview:null, importRequest:0, reuseOriginal:null, reuseBaseline:null};
const statuses = {success:"已完成", failed:"失败", partial_failure:"部分失败", running:"运行中", queued:"准备中", interrupted:"已中断", unknown:"状态待检查", created:"已创建", ready:"待运行", unreadable:"记录不可读", blocked:"已阻断"};
const capabilities = {dataset_audit:"数据与实验设计检查", qc:"质量控制", normalization:"归一化", integration:"批次校正", clustering:"聚类与细胞注释", subset_cells:"选择目标细胞", deg:"差异表达", differential_abundance:"细胞丰度比较", trajectory_inference:"拟时序分析", fate_prediction:"细胞命运推断", cellrank_fate:"细胞命运推断", donor_sensitivity:"供体留一敏感性", functional_activity:"通路与转录因子活性", doublet_detection:"双细胞检测", cell_annotation:"参考模型注释", background_removal:"背景去除", cell_cell_communication:"细胞通讯", knowledge_retrieval:"知识检索", spatial_domain:"空间分区", spatial_deg:"空间差异分析"};
const eventLabels = {run_started:"开始运行", plan_created:"执行计划已生成", task_started:"正在执行", method_resolved:"分析方法已选择", attempt_started:"开始计算", attempt_finished:"计算完成，等待审计", artifacts_committed:"分析产物已保存", audit_started:"正在独立审计", audit_finished:"审计已通过", audit_rejected:"审计未通过", task_finished:"任务完成", task_failed:"任务失败", task_blocked:"任务被阻断", task_resumed:"已复用保存的计算结果", plan_adapted:"计划已根据数据调整", run_finished:"本轮执行结束", run_failed:"本轮执行失败", retry_scheduled:"准备重试"};
function el(tag, text, cls){const node=document.createElement(tag);if(text!==undefined&&text!==null)node.textContent=String(text);if(cls)node.className=cls;return node;}
function notice(text, type=""){return el("div",text,"notice "+type);}
function badge(value){return el("span",statuses[value]||value||"未知","status "+(Object.hasOwn(statuses,value)?value:""));}
function toast(message,error=false){$("toast").textContent=message;$("toast").className="toast"+(error?" error":"");clearTimeout(toast.timer);toast.timer=setTimeout(()=>$("toast").classList.add("hidden"),error?12000:5000);}
function errorText(error){return [error.message||"操作失败",error.detail?"技术详情："+error.detail:"",error.nextSteps?.length?"下一步："+error.nextSteps.join("；"):""].filter(Boolean).join("\n");}
function presentError(error){const message=error.message||String(error);$("operation-error-title").textContent=message;$("operation-error-detail").textContent=error.detail?"技术详情："+error.detail:"";$("operation-error-next").textContent=error.nextSteps?.length?"下一步："+error.nextSteps.join("；"):"";$("operation-error").classList.remove("hidden");toast(errorText(error),true);}
$("close-operation-error").onclick=()=>$("operation-error").classList.add("hidden");
async function api(path, payload){const response=await fetch(path,{method:payload===undefined?"GET":"POST",headers:{"X-EACBP-Token":token,...(payload===undefined?{}:{"Content-Type":"application/json"})},body:payload===undefined?undefined:JSON.stringify(payload)});const data=await response.json();if(!response.ok){const error=Error(typeof data.error==="string"?data.error:data.error?.message||"请求失败");error.detail=data.detail;error.nextSteps=data.next_steps||[];throw error;}return data;}
function handle(fn){return async event=>{try{await fn(event);}catch(error){presentError(error);}};}
async function withButton(button,fn){button.disabled=true;try{return await fn();}finally{button.disabled=false;if(state.view==="detail"&&state.detail)renderDetail(state.detail);syncStart();}}
function syncStart(){$("start").disabled=state.busy||!state.preview?.can_run||!state.designConfirmed||state.jobs.some(j=>["running","queued","unknown"].includes(j.status));}
function invalidate(){const hadPreview=Boolean(state.preview||state.previewKey);state.preview=null;state.previewKey=null;state.designConfirmed=false;updateGuidance();if(!hadPreview&&$("plan-count").textContent==="待生成"){syncStart();return;}$("plan-count").textContent="待更新";$("plan-content").replaceChildren(notice("配置已修改，请重新预览执行计划。"));$("plan-content").classList.remove("hidden");$("plan-empty").classList.add("hidden");syncStart();}
function navigate(view){state.view=view;for(const name of ["new","history","detail"])$(name+"-view").classList.toggle("hidden",name!==view);$("nav-new").classList.toggle("active",view==="new");$("nav-history").classList.toggle("active",view!=="new");$("breadcrumb").textContent=view==="new"?"新建分析":view==="history"?"运行记录":"运行详情";window.scrollTo({top:0});}
function downloadBlob(text,name,type="text/plain;charset=utf-8"){const url=URL.createObjectURL(text instanceof Blob?text:new Blob([text],{type}));const a=el("a");a.href=url;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);}
async function downloadEndpoint(path,name){const response=await fetch(path,{headers:{"X-EACBP-Token":token}});if(!response.ok){const data=await response.json(),error=Error(data.error||"下载失败");error.detail=data.detail;error.nextSteps=data.next_steps||[];throw error;}downloadBlob(await response.blob(),name);}
async function downloadRun(name){if(!state.selected)return;await downloadEndpoint("/api/download?id="+encodeURIComponent(state.selected)+"&name="+encodeURIComponent(name),name);}
async function browse(path="."){const listing=await api("/api/files?path="+encodeURIComponent(path));state.folder=listing;$("folder-path").textContent=listing.path;$("parent-folder").disabled=!listing.parent;$("file-list").replaceChildren();if(!listing.entries.length)$("file-list").append(el("div","这里没有可选择的数据文件（h5ad、CSV 或 TSV）。","empty-state"));for(const file of listing.entries){const button=el("button",null,"file-entry");button.append(el("span",file.directory?"▣":"▤"),el("span",file.name));if(!file.directory)button.append(el("small",(file.size/1024/1024).toFixed(1)+" MB"));button.addEventListener("click",handle(async()=>{if(file.directory)await browse(file.path);else if(state.fileTarget){$(state.fileTarget).value=file.path;state.fileTarget=null;clearImportPreview();$("file-dialog").close();}else if(file.name.toLocaleLowerCase().endsWith(".h5ad")){$("data").value=file.path;clearDataset({preserveMappings:Boolean(state.reuseOriginal)});invalidate();$("file-dialog").close();}else toast("CSV/TSV 输入请使用“导入 CSV/TSV 为 h5ad”向导。",true);}));$("file-list").append(button);}if(listing.truncated)$("file-list").append(el("small","仅显示前 300 项；也可以直接输入完整路径。"));}
function clearImportPreview(){state.importPreview=null;$("confirm-import").disabled=true;$("import-preview").replaceChildren();$("import-preview").classList.add("hidden");$("import-error").replaceChildren();$("import-error").classList.add("hidden");}
function importPayload(){return {counts:$("import-counts").value.trim(),metadata:$("import-metadata").value.trim(),orientation:document.querySelector('input[name="import-orientation"]:checked')?.value};}
function showImportPreview(data){const root=$("import-preview");root.replaceChildren();root.classList.remove("hidden");root.append(el("strong","确认导入预览"),el("p",`${data.n_cells.toLocaleString()} 个细胞 × ${data.n_genes.toLocaleString()} 个基因 · ${data.orientation==="cells_by_genes"?"行=细胞，列=基因":"行=基因，列=细胞"}`),el("p",data.alignment+(data.metadata_rows_reordered?" 元数据行顺序不同，已按细胞 ID 对齐。":" 元数据顺序与计数矩阵一致。")),el("p",data.count_check),el("p",`计数范围 ${data.count_range.minimum}–${data.count_range.maximum}，非零值 ${data.count_range.nonzero.toLocaleString()}。`));const missing=Object.entries(data.metadata_missing||{}).filter(([,count])=>count>0);root.append(el("p",missing.length?"元数据空值（确认后转换为缺失值）：":"元数据列没有空值。"));if(missing.length){const list=el("ul");for(const [column,count] of missing)list.append(el("li",`${column}: ${count.toLocaleString()} 项`));root.append(list);}root.append(el("small",`SHA-256 · counts ${data.counts_sha256} · metadata ${data.metadata_sha256}`),el("p","确认后会在当前工作目录的 .eacbp/imports 下创建新 h5ad；源 CSV/TSV 不会修改。请核对维度、方向、ID 对齐和计数检查后再确认。"));}
async function previewImport(){clearImportPreview();$("preview-import").disabled=true;const sequence=++state.importRequest,payload=importPayload(),key=JSON.stringify(payload);try{const data=await api("/api/import/preview",payload);if(sequence!==state.importRequest||JSON.stringify(importPayload())!==key)return;state.importPreview=data;showImportPreview(data);$("confirm-import").disabled=false;}finally{$("preview-import").disabled=false;}}
async function confirmImport(){if(!state.importPreview?.preview_token)throw Error("请先预检并确认导入预览");const button=$("confirm-import");button.disabled=true;button.textContent="正在创建 h5ad…";try{const result=await api("/api/import/confirm",{preview_token:state.importPreview.preview_token});state.importPreview=null;$("data").value=result.data;clearDataset();invalidate();$("import-dialog").close();toast(result.notice);await inspectDataset();}catch(error){clearImportPreview();presentError(error);$("import-error").textContent=errorText(error);$("import-error").classList.remove("hidden");}finally{button.textContent="确认并新建 h5ad";button.disabled=!state.importPreview;}}
$("csv-import-open").onclick=()=>$("import-dialog").showModal();
for(const [buttonId,target] of [["browse-import-counts","import-counts"],["browse-import-metadata","import-metadata"]])$(buttonId).onclick=handle(async()=>{state.fileTarget=target;await browse();$("file-dialog").showModal();});
for(const id of ["import-counts","import-metadata"])$(id).addEventListener("input",()=>{state.importRequest++;clearImportPreview();});
for(const node of document.querySelectorAll('input[name="import-orientation"]'))node.addEventListener("change",()=>{state.importRequest++;clearImportPreview();});
$("preview-import").onclick=handle(previewImport);$("confirm-import").onclick=handle(confirmImport);
for(const id of ["close-import","cancel-import"])$(id).onclick=()=>$("import-dialog").close();
$("close-files").onclick=()=>{$("file-dialog").close();state.fileTarget=null;};
function stat(label,value){const node=el("div",null,"stat");node.append(el("small",label),typeof value==="object"?value:el("strong",value));return node;}
function renderHistory(){const allRows=state.runs;const query=$("run-search").value.trim().toLocaleLowerCase();const filter=$("run-filter").value;const rows=allRows.filter(r=>(!query||[r.title,r.study_id,r.id].filter(Boolean).join(" ").toLocaleLowerCase().includes(query))&&(filter==="all"||(filter==="success"&&r.status==="success")||(filter==="active"&&["running","queued"].includes(r.status))||(filter==="attention"&&["failed","partial_failure","interrupted","unknown","unreadable","blocked"].includes(r.status))));$("run-filter-count").textContent=allRows.length?"显示 "+rows.length+" / "+allRows.length+" 次运行":"";$("run-count").textContent=allRows.length;$("history-stats").replaceChildren(stat("全部运行",allRows.length),stat("已完成",allRows.filter(r=>r.status==="success").length),stat("正在运行",allRows.filter(r=>["running","queued"].includes(r.status)).length),stat("需要检查",allRows.filter(r=>["failed","partial_failure","interrupted","unknown","unreadable"].includes(r.status)).length));const root=$("runs-list");root.replaceChildren();if(!rows.length){const empty=el("div",null,"empty-state");empty.append(el("strong",allRows.length?"没有匹配的运行记录":"你的第一份分析记录，将从这里开始"),el("p",allRows.length?"试试其他关键词或状态筛选。":"创建分析后，可在这里跟踪进度并查看结果。"));root.append(empty);return;}const table=el("table",null,"runs-table");const head=el("thead"),tr=el("tr");["研究 / 运行目录","状态","最近更新",""].forEach(t=>tr.append(el("th",t)));head.append(tr);table.append(head);const body=el("tbody");for(const run of rows){const row=el("tr"),title=el("td"),button=el("button",null,"run-open");button.append(el("strong",run.title||run.study_id),el("small",run.id));button.addEventListener("click",handle(()=>openRun(run.id)));title.append(button);const status=el("td");status.append(badge(run.status));const date=el("td",run.updated_at?new Date(run.updated_at).toLocaleString():"—");const action=el("td"),go=el("button","查看 →","text-button");go.addEventListener("click",handle(()=>openRun(run.id)));action.append(go);row.append(title,status,date,action);body.append(row);}table.append(body);root.append(table);}
async function refreshRuns(){const data=await api("/api/runs");state.runs=data.runs;state.jobs=data.jobs;const key=JSON.stringify(data.runs);if(state.historyKey!==key){state.historyKey=key;renderHistory();}syncStart();syncDirectoryControls();}
async function openRun(id){state.selected=id;state.detail=null;state.resultsId=null;state.resultsLoading=null;$("refresh-results").disabled=false;state.resultsRequest=(state.resultsRequest||0)+1;navigate("detail");selectDetailTab("events");for(const id of ["detail-notice","detail-stats","report-output","config-output","log-output","results-content"])$(id).replaceChildren();for(const id of ["resume","rebuild-report","download-report","download-config","download-summary"])$(id).disabled=true;$("detail-title").textContent="正在读取运行…";$("detail-path").textContent=id;$("events-panel").replaceChildren(el("p","正在加载…"));await refreshDetail();}
function inline(parent,text){const regex=/(\*\*([^*]+)\*\*|`([^`]+)`)/g;let last=0;for(const match of text.matchAll(regex)){parent.append(document.createTextNode(text.slice(last,match.index)));parent.append(el(match[2]?"strong":"code",match[2]||match[3]));last=match.index+match[0].length;}parent.append(document.createTextNode(text.slice(last)));}
function renderMarkdown(text){const root=el("div");let code=null,list=null,table=null;const lines=text.split("\n");for(let i=0;i<lines.length;i++){const line=lines[i];if(line.trim().startsWith("```")){if(code){root.append(el("pre",code.join("\n")));code=null;}else code=[];list=null;table=null;continue;}if(code){code.push(line);continue;}if(!line.trim()){list=null;table=null;continue;}const heading=line.match(/^(#{1,6})\s+(.+)/);if(heading){const node=el("h"+heading[1].length);inline(node,heading[2]);root.append(node);list=null;table=null;continue;}if(/^\s*\|.*\|\s*$/.test(line)){const cells=line.trim().slice(1,-1).split("|").map(c=>c.trim());if(cells.every(c=>/^:?-+:?$/.test(c)))continue;if(!table){table=el("table");root.append(table);}const row=el("tr");for(const cell of cells){const td=el(i+1<lines.length&&/^\s*\|[\s:|\-]+\|\s*$/.test(lines[i+1])?"th":"td");inline(td,cell);row.append(td);}table.append(row);continue;}table=null;if(/^\s*[-*]\s+/.test(line)){if(!list){list=el("ul");root.append(list);}const li=el("li");inline(li,line.replace(/^\s*[-*]\s+/,""));list.append(li);continue;}list=null;if(/^\s*[-_]{3,}\s*$/.test(line)){root.append(el("hr"));continue;}const p=el("p");inline(p,line);root.append(p);}if(code)root.append(el("pre",code.join("\n")));return root;}
function renderDetail(data){const prior=state.detail;state.detail=data;$("detail-title").textContent=data.title||data.study_id;$("detail-path").textContent=data.path;const active=state.jobs.some(j=>["running","queued","unknown"].includes(j.status));$("resume").disabled=active||!data.can_resume;$("rebuild-report").disabled=active||!data.can_report;$("download-report").disabled=!data.has_report;$("download-config").disabled=false;$("download-summary").disabled=!data.summary;const notifications=$("detail-notice");notifications.replaceChildren();if(data.job?.error)notifications.append(notice(data.job.error,"error"));if(["running","queued"].includes(data.status))notifications.append(notice("后台作业正在运行。页面每 3 秒更新，关闭浏览器后作业仍会继续。"));if(data.job?.operation==="report"&&data.job.status==="success")notifications.append(notice("本次报告校验与重建已完成。报告中的分析结论仍需结合其科学边界解读。"));const events=data.events||[];const terminal=new Set(events.filter(e=>["task_finished","task_failed","task_blocked"].includes(e.kind)).map(e=>e.task_id));const resumed=new Set(events.filter(e=>e.kind==="task_resumed").map(e=>e.task_id));const summary=data.summary||{};const summaryCurrent=!data.job||data.job.operation==="report"||data.job.status==="success";$("detail-stats").replaceChildren(stat("运行状态",badge(data.status)),stat("已处理步骤",summaryCurrent?(summary.tasks_executed??terminal.size):terminal.size),stat("本轮复用",resumed.size),stat("已阻断",summaryCurrent?(summary.tasks_blocked??"—"):events.filter(e=>e.kind==="task_blocked").length));const panel=$("events-panel");panel.replaceChildren();const meaningful=events.filter(e=>!["method_resolved","attempt_started","attempt_finished","artifacts_committed"].includes(e.kind));if(!meaningful.length)panel.append(el("div","等待数据导入或后台日志。大型文件导入可能需要一些时间。","empty-state"));for(const event of meaningful.slice(-100).reverse()){const row=el("div",null,"event-row"),body=el("div",null,"event-body");row.append(el("span",new Date(event.timestamp).toLocaleTimeString(),"event-time"));body.append(el("strong",eventLabels[event.kind]||event.kind));const details=event.details||{};body.append(el("small",[event.task_id,capabilities[details.capability]||details.capability,details.method,details.error,details.reason].filter(Boolean).join(" · ")));if(event.kind==="plan_adapted")body.append(el("pre",JSON.stringify(details,null,2)));row.append(body);panel.append(row);}if($("log-output").textContent!==data.log){const atBottom=$("log-output").scrollTop+$("log-output").clientHeight>=$("log-output").scrollHeight-30;$("log-output").textContent=data.log||"暂无后台进程日志。已有 CLI 运行可查看任务事件。";if(atBottom)$("log-output").scrollTop=$("log-output").scrollHeight;}if(!prior||prior.report!==data.report){$("report-output").replaceChildren(data.report?renderMarkdown(data.report):el("p","报告将在分析结束并完成快照校验后生成。"));}const configText=JSON.stringify({summary:data.summary,configuration:data.saved},null,2);if($("config-output").textContent!==configText)$("config-output").textContent=configText;}
function syncDetailActions(data){const active=state.jobs.some(j=>["running","queued","unknown"].includes(j.status)),complete=data.saved?.status==="success"&&data.saved?.import_completed===true;$("reuse-run").disabled=active||!complete;$("download-evidence").disabled=active||!complete;$("download-snapshot").disabled=active||!data.can_report;$("download-methods").disabled=active||!data.can_report;}
async function refreshDetail(){const id=state.selected;if(!id)return;const data=await api("/api/run?id="+encodeURIComponent(id));if(state.selected===id&&state.view==="detail"){renderDetail(data);syncDetailActions(data);}}
async function jobAction(operation){const id=state.selected;if(!id)return;await api("/api/"+operation,{run_id:id});toast(operation==="resume"?"已提交恢复运行；将校验原有配置与产物。":"已提交报告校验与重建。");await refreshRuns();await refreshDetail();}
$("nav-new").onclick=()=>navigate("new");$("history-new").onclick=()=>navigate("new");$("nav-history").onclick=handle(async()=>{navigate("history");await refreshRuns();});$("back-history").onclick=handle(async()=>{navigate("history");await refreshRuns();});$("refresh-runs").onclick=handle(refreshRuns);
$("study-form").addEventListener("input",invalidate);$("study-form").addEventListener("change",invalidate);$("species").addEventListener("change",()=>$("other-species").classList.toggle("hidden",$("species").value!=="other"));$("data").addEventListener("input",()=>clearDataset({preserveMappings:Boolean(state.reuseOriginal)}));
$("study-form").addEventListener("submit",handle(async event=>{event.preventDefault();if(!$("study-form").reportValidity())return;await withButton($("preview"),async()=>{const submitted=JSON.stringify(formValue());$("preview").textContent="正在生成计划…";try{const prepared=await api("/api/preview",JSON.parse(submitted));if(JSON.stringify(formValue())!==submitted){toast("配置已改变，请重新预览。",true);return;}renderPlan(prepared);toast(prepared.can_run?"计划已生成，请查看步骤和提示。":"请检查计划中的问题。",!prepared.can_run);}finally{$("preview").textContent="预览分析计划 →";}});}));
$("start").addEventListener("click",handle(async()=>{if(!state.preview||JSON.stringify(formValue())!==state.previewKey)throw Error("配置已修改，请重新预览");if(!state.designConfirmed)throw Error("请先核对并确认实验设计摘要");state.busy=true;syncStart();try{const response=await api("/api/run",formValue());toast("分析已在后台启动。");await refreshRuns();await openRun(response.run_id);}finally{state.busy=false;syncStart();}}));
$("browse").onclick=handle(async()=>{await browse();$("file-dialog").showModal();});$("parent-folder").onclick=handle(()=>browse(state.folder.parent));$("close-files").onclick=()=>$("file-dialog").close();$("inspect-data").onclick=handle(inspectDataset);
$("export-config").onclick=handle(()=>{downloadBlob(JSON.stringify({webui_schema:1,form:formValue()},null,2),"eacbp-settings.json");toast("设置已导出。数据文件本身未包含在内。");});$("import-config").onclick=()=>$("config-file").click();$("config-file").onchange=handle(async()=>{const file=$("config-file").files[0];if(!file)return;try{if(file.size>1000000)throw Error("设置文件不能超过 1 MB");const value=JSON.parse(await file.text());if(value.webui_schema!==1||!value.form||Array.isArray(value.form)||typeof value.form!=="object")throw Error("请选择通过工作台导出的设置文件");fillForm(value.form);toast("设置已导入，请检查路径并重新预览。");}finally{$("config-file").value="";}});
$("resume").onclick=handle(()=>withButton($("resume"),()=>jobAction("resume")));$("rebuild-report").onclick=handle(()=>withButton($("rebuild-report"),()=>jobAction("report")));$("reuse-run").onclick=handle(reuseRun);$("download-report").onclick=handle(()=>downloadRun("report.md"));$("download-config").onclick=handle(()=>downloadRun("config.json"));$("download-summary").onclick=handle(()=>downloadRun("summary.json"));$("download-evidence").onclick=handle(()=>downloadEndpoint("/api/export?id="+encodeURIComponent(state.selected),state.selected.replaceAll("/","-")+"-evidence.zip"));$("download-snapshot").onclick=handle(()=>downloadEndpoint("/api/snapshot?id="+encodeURIComponent(state.selected),state.selected.replaceAll("/","-")+"-snapshot.json"));$("download-methods").onclick=handle(()=>downloadEndpoint("/api/methods?id="+encodeURIComponent(state.selected),state.selected.replaceAll("/","-")+"-methods.md"));
function selectDetailTab(name){
  for(const button of document.querySelectorAll(".tab")){
    const selected=button.dataset.tab===name;
    button.classList.toggle("active",selected);
    button.setAttribute("aria-selected",String(selected));
    button.tabIndex=selected?0:-1;
    $(button.dataset.tab+"-panel").classList.toggle("hidden",!selected);
  }
  if(name==="results"&&typeof loadResults==="function")loadResults();
}
const detailTabs=[...document.querySelectorAll(".tab")];
for(const [index,button] of detailTabs.entries()){
  const name=button.dataset.tab;
  button.id=name+"-tab";
  button.setAttribute("aria-controls",name+"-panel");
  $(name+"-panel").setAttribute("role","tabpanel");
  $(name+"-panel").setAttribute("aria-labelledby",button.id);
  button.onclick=()=>selectDetailTab(name);
  button.onkeydown=event=>{
    if(!["ArrowLeft","ArrowRight","Home","End"].includes(event.key))return;
    event.preventDefault();
    const next=event.key==="Home"?0:event.key==="End"?detailTabs.length-1:(index+(event.key==="ArrowRight"?1:-1)+detailTabs.length)%detailTabs.length;
    detailTabs[next].focus();selectDetailTab(detailTabs[next].dataset.tab);
  };
}
selectDetailTab("events");
async function bootstrap(){const settings=await api("/api/settings");if(new URLSearchParams(window.location.search).get("directories")==="updated"){$("study-form").reset();$("other-species").classList.add("hidden");clearDataset();window.history.replaceState(null,"",window.location.pathname);toast("目录设置已应用，请重新选择数据并预览分析计划。");}$("environment").textContent="Python 环境已连接";const detail=$("environment-details");detail.replaceChildren();detail.append(el("small","输入工作目录"),el("code",settings.workspace),el("small","结果与历史目录"),el("code",settings.runs_dir));for(const [name,ver] of Object.entries(settings.packages))detail.append(el("small",name+" · "+(ver||"未安装")));await refreshRuns();}
bootstrap().catch(error=>{$("environment").textContent="连接失败";presentError(error);});
let polling=false;setInterval(async()=>{if(polling||document.hidden||state.directorySaving)return;polling=true;try{await refreshRuns();if(state.view==="detail")await refreshDetail();$("connection-error").classList.add("hidden");}catch(error){$("connection-error").textContent="工作台连接异常：\n"+errorText(error);$("connection-error").classList.remove("hidden");}finally{polling=false;}},3000);

$("run-search").addEventListener("input",renderHistory);
$("run-filter").addEventListener("change",renderHistory);
function syncDirectoryControls(){
  const active=state.jobs.some(job=>["queued","running","unknown"].includes(job.status));
  $("save-directories").disabled=Boolean(state.directorySaving)||active;
  $("directory-workspace").disabled=Boolean(state.directorySaving);
  $("directory-runs").disabled=Boolean(state.directorySaving);
  $("directory-job-notice").classList.toggle("hidden",!active);
}
async function openDirectorySettings(){
  const settings=await api("/api/settings");
  await refreshRuns();
  $("directory-workspace").value=settings.workspace;
  $("directory-runs").value=settings.runs_dir;
  $("directory-package").textContent=settings.package_dir||"—";
  $("directory-python").textContent=settings.python_executable||"—";
  $("directory-preferences").textContent=settings.settings_file||"当前会话（未启用持久保存）";
  $("save-directories").textContent=settings.settings_file?"保存并应用":"应用到当前会话";
  $("directory-error").classList.add("hidden");
  syncDirectoryControls();
  $("directory-dialog").showModal();
}
$("nav-directories").onclick=handle(openDirectorySettings);
for(const id of ["close-directories","cancel-directories"])$(id).onclick=()=>{if(!state.directorySaving)$("directory-dialog").close();};
$("directory-dialog").addEventListener("cancel",event=>{if(state.directorySaving)event.preventDefault();});
$("directory-form").addEventListener("submit",async event=>{
  event.preventDefault();
  if(state.directorySaving||!$("directory-form").reportValidity())return;
  state.directorySaving=true;syncDirectoryControls();
  const button=$("save-directories"),prior=button.textContent;
  button.textContent="正在应用…";
  $("directory-error").classList.add("hidden");
  try{
    await api("/api/settings",{workspace:$("directory-workspace").value.trim(),runs_dir:$("directory-runs").value.trim()});
    // The server invalidates all old page tokens. Reload also discards old
    // previews, file selections and run IDs instead of reusing another root.
    window.location.replace("/?directories=updated");
  }catch(error){
    $("directory-error").textContent=errorText(error);
    $("directory-error").classList.remove("hidden");
    state.directorySaving=false;button.textContent=prior;syncDirectoryControls();
  }
});

// Guided research purpose, explicit metadata mappings and the authoritative pre-run checklist.
const purposeInfo={
  overview:{title:"数据与细胞组成",tasks:"数据与实验设计检查 · 质量控制 · 聚类与细胞注释",inputs:"h5ad 文件、物种、组织；已有细胞标签可用于目标细胞筛选。",outputs:"QC 保留情况、聚类与细胞状态结果。"},
  deg:{title:"两组基因表达比较",tasks:"差异表达（DEG）",inputs:"明确的条件列、两个条件值；供体列可支持供体 pseudobulk，PyDESeq2 需要原始整数计数。",outputs:"基因差异表；统计单位和可用重复数由后端审计决定。"},
  abundance:{title:"细胞类型丰度比较",tasks:"差异丰度（differential_abundance）",inputs:"明确的条件列、两个条件值、状态标签、每组至少两个独立供体。",outputs:"按供体汇总的细胞状态比例比较；现有方法不支持配对丰度检验。"},
  trajectory:{title:"细胞状态轨迹探索",tasks:"拟时序分析（trajectory_inference）",inputs:"标准分析模式与一个存在且有生物学依据的根细胞 ID。",outputs:"相对根细胞的拟时序及相关基因结果；方向不等于谱系因果。"},
  batch:{title:"批次影响检查",tasks:"批次校正（integration）与数据审计",inputs:"明确的批次列；当前后端按所选 batch_col 执行现有整合方法。",outputs:"整合后的表示和批次混合指标。"},
  annotation:{title:"本地参考模型注释",tasks:"参考模型注释（cell_annotation / celltypist_local_v1）",inputs:"运行机器上的 CellTypist 模型文件与已安装 celltypist 包。",outputs:"预测标签、冲突标记和注释来源记录。"}
};
state.pendingMappings={};state.pendingConditions={};state.designRequest=0;
function selectedPurpose(){return document.querySelector('input[name="research_purpose"]:checked')?.value||"overview";}
const reuseConfigPaths={
  method_profile:[["method_profile"]],
  advanced_analysis:[["advanced_analysis"],["capability_parameters","deg","counts_layer"],["capability_parameters","deg","allow_x_as_counts"],["capability_parameters","deg","min_donors"],["capability_parameters","deg","alpha"],["capability_parameters","deg","design_formula"]],
  min_genes:[["capability_parameters","qc","min_genes"]],
  max_mito_pct:[["capability_parameters","qc","max_mito_pct"]],
  condition_col:[...["dataset_audit","deg","differential_abundance"].map(name=>["capability_parameters",name,"condition_col"])],
  condition_a:[...["deg","differential_abundance"].map(name=>["capability_parameters",name,"condition_a"])],
  condition_b:[...["deg","differential_abundance"].map(name=>["capability_parameters",name,"condition_b"])],
  donor_col:[...["dataset_audit","deg","differential_abundance"].map(name=>["capability_parameters",name,"donor_col"])],
  batch_col:[...["dataset_audit","deg","differential_abundance","integration"].map(name=>["capability_parameters",name,"batch_col"])],
  paired:[["capability_parameters","deg","paired"]],
  counts_layer:[["capability_parameters","deg","counts_layer"]],
  allow_x_as_counts:[["capability_parameters","deg","allow_x_as_counts"]],
  min_donors:[["capability_parameters","deg","min_donors"]],
  alpha:[["capability_parameters","deg","alpha"]],
  design_formula:[["capability_parameters","deg","design_formula"]],
  root_cell_id:[["capability_parameters","trajectory_inference","root_cell_id"]],
  annotation_use_as_cell_type:[["analysis_extensions","cell_annotation","use_as_cell_type"]],
  celltypist_model_path:[["analysis_extensions","cell_annotation","model_path"]]
};
function removeConfigPath(root,path){let node=root;for(const key of path.slice(0,-1)){if(!node||typeof node!=="object")return;node=node[key];}if(node&&typeof node==="object")delete node[path.at(-1)];}
function formValue(){
  const value={};
  for(const key of ["data","study_id","title","species","tissue","condition_a","condition_b","condition_col","sample_col","donor_col","batch_col","root_cell_id","method_profile","design_formula","counts_layer","celltypist_model_path"]){let current=$(key).value.trim();if(state.reuseOriginal&&!current&&Object.hasOwn(state.pendingMappings||{},key))current=state.pendingMappings[key]||"";if(state.reuseOriginal&&!current&&Object.hasOwn(state.pendingConditions||{},key))current=state.pendingConditions[key]||"";value[key]=current;}
  value.research_purpose=selectedPurpose();
  if(value.species==="other")value.species=$("other-species").value.trim();
  for(const key of ["advanced_analysis","paired","allow_x_as_counts","annotation_use_as_cell_type"])value[key]=$(key).checked;
  for(const key of ["min_genes","max_mito_pct","min_donors","alpha"])value[key]=Number($(key).value);
  value.target_cell_types=$("target_cell_types").value.split(/[,，]/).map(v=>v.trim()).filter(Boolean);
  try{value.config_overrides=JSON.parse($("config_overrides").value||"{}");}catch{$("config_overrides").closest("details").open=true;$("config_overrides").focus();throw Error("高级配置 JSON 格式不正确，请检查逗号、引号和括号。");}
  if(!value.config_overrides||Array.isArray(value.config_overrides)||typeof value.config_overrides!=="object")throw Error("高级配置必须是 JSON 对象");
  if(state.reuseBaseline){for(const [field,paths] of Object.entries(reuseConfigPaths))if(JSON.stringify(value[field])!==JSON.stringify(state.reuseBaseline[field]))for(const path of paths)removeConfigPath(value.config_overrides,path);}
  return value;
}
function fillForm(value){
  $("study-form").reset();
  const defaults={advanced_analysis:false,paired:false,allow_x_as_counts:false,annotation_use_as_cell_type:false,min_donors:2,alpha:0.05,design_formula:"",counts_layer:"counts",celltypist_model_path:"",research_purpose:"overview",config_overrides:{}};
  for(const [key,item] of Object.entries(defaults)){
    if(key==="research_purpose")continue;
    if(["advanced_analysis","paired","allow_x_as_counts","annotation_use_as_cell_type"].includes(key))$(key).checked=item;
    else $(key).value=typeof item==="object"?JSON.stringify(item,null,2):String(item);
  }
  state.pendingMappings={condition_col:"",sample_col:"",donor_col:"",batch_col:"",counts_layer:"counts"};
  state.pendingConditions={condition_a:"",condition_b:""};
  const defaultPurpose=document.querySelector('input[name="research_purpose"][value="overview"]');if(defaultPurpose)defaultPurpose.checked=true;
  $("research-purpose-placeholder")?.remove();
  for(const [key,item] of Object.entries(value||{})){
    if(key==="research_purpose"){
      const radio=[...document.querySelectorAll('input[name="research_purpose"]')].find(node=>node.value===String(item));if(radio)radio.checked=true;continue;
    }
    if(Object.hasOwn(state.pendingMappings,key)){state.pendingMappings[key]=String(item??"");continue;}
    if(Object.hasOwn(state.pendingConditions,key)){state.pendingConditions[key]=String(item??"");continue;}
    if(!$(key))continue;
    if(["advanced_analysis","paired","allow_x_as_counts","annotation_use_as_cell_type"].includes(key))$(key).checked=item===true;
    else if(key==="target_cell_types")$(key).value=Array.isArray(item)?item.join(", "):String(item);
    else if(key==="config_overrides")$(key).value=JSON.stringify(item,null,2);
    else if(key==="species"){const known=["homo_sapiens","mus_musculus"].includes(item);$(key).value=known?item:"other";$("other-species").value=known?"":item;$("other-species").classList.toggle("hidden",known);}
    else if(["condition_a","condition_b"].includes(key))state.pendingConditions[key]=String(item??"");
    else $(key).value=String(item??"");
  }
  clearDataset({preserveMappings:true,skipCapture:true});invalidate();updatePurposeGuide();
}
function jsonLeaves(value,prefix="",out={}){if(value&&typeof value==="object"&&!Array.isArray(value)){for(const key of Object.keys(value).sort())jsonLeaves(value[key],prefix?prefix+"."+key:key,out);}else out[prefix||"$"]=JSON.stringify(value);return out;}
function renderReuseComparison(prepared=null){const root=$("reuse-comparison"),reuse=state.reuseOriginal;if(!reuse){root.replaceChildren();root.classList.add("hidden");return;}root.replaceChildren(el("h2","成功运行配置复用对照"),el("p",reuse.notice),el("p",`原运行：${reuse.original.study_id} · ${reuse.original.run_id}；新研究编号：${$("study_id").value}。原输入：${reuse.original.source_path||"未记录"}；新输入路径保持空白，需重新选择并运行数据检查、计划预览和设计确认。`));if(!reuse.research_purpose_saved)root.append(notice("历史记录未保存研究问题；当前使用“数据与细胞组成”作为页面默认显示，未写回历史配置。","warning"));for(const item of reuse.cleared_external_inputs||[])root.append(notice(`${item.name} 已清除：${item.path}。${item.reason}`,"warning"));if(reuse.unsupported_manifest_sections?.length){const list=el("ul",null,"reuse-diff");for(const item of reuse.unsupported_manifest_sections)list.append(el("li",item));root.append(el("p","新建时不会静默继承以下旧 manifest 字段："),list);}const original=reuse.original.config;let current=prepared?.config;if(!current){try{current=JSON.parse($("config_overrides").value||"{}");}catch{current={};}}const before=jsonLeaves(original),after=jsonLeaves(current),changes=[];for(const key of [...new Set([...Object.keys(before),...Object.keys(after)])].sort())if(before[key]!==after[key])changes.push(`${key}: ${before[key]??"（未设置）"} → ${after[key]??"（未设置）"}`);const heading=el("p",prepared?"本次预检的最终配置与历史最终配置相比：":"当前保留的完整配置与历史最终配置相比：");root.append(heading);if(changes.length){const list=el("ul",null,"reuse-diff");for(const change of changes.slice(0,80))list.append(el("li",change));if(changes.length>80)list.append(el("li",`另有 ${changes.length-80} 项差异。`));root.append(list);}else root.append(el("p","没有配置差异；研究编号、manifest 身份和输入路径仍按新运行重新生成。"));const detail=el("details",null,"reuse-original"),summary=el("summary","查看原运行保存的完整最终配置");detail.append(summary,el("pre",JSON.stringify(original,null,2)));root.append(detail);root.classList.remove("hidden");}
async function reuseRun(){const id=state.selected;if(!id)return;await withButton($("reuse-run"),async()=>{const result=await api("/api/reuse?id="+encodeURIComponent(id));if(state.selected!==id||state.view!=="detail")return;state.reuseOriginal=result;state.reuseBaseline=null;fillForm(result.form);state.reuseBaseline=Object.fromEntries(Object.keys(reuseConfigPaths).map(key=>[key,result.form[key]]));renderReuseComparison();navigate("new");$("reuse-comparison").scrollIntoView({block:"start",behavior:"smooth"});toast("已载入可映射配置；请选择当前数据并重新检查、预览和确认。");});}
function clearDataset(options={}){
  const keep=Boolean(options.preserveMappings);
  if(keep){if(!options.skipCapture){for(const key of ["condition_col","sample_col","donor_col","batch_col","counts_layer"])if($(key)?.value)state.pendingMappings[key]=$(key).value;for(const key of ["condition_a","condition_b"])if($(key)?.value)state.pendingConditions[key]=$(key).value;}}
  else{state.pendingMappings={condition_col:"",sample_col:"",donor_col:"",batch_col:"",counts_layer:"counts"};state.pendingConditions={condition_a:"",condition_b:""};for(const key of ["condition_col","sample_col","donor_col","batch_col"])if($(key))$(key).value="";for(const key of ["condition_a","condition_b"])if($(key))$(key).value="";}
  state.dataset=null;state.designRequest++;$("dataset-info").classList.add("hidden");$("design-overview").replaceChildren();$("design-overview").classList.add("hidden");
  for(const key of ["condition_col","sample_col","donor_col","batch_col"]){if(!$(key))continue;$(key).replaceChildren(new Option(key==="sample_col"?"不显示样本数":"请选择数据列",""));}
  $("condition_a").replaceChildren(new Option("先选择条件列",""));$("condition_b").replaceChildren(new Option("先选择条件列",""));$("condition_a").disabled=true;$("condition_b").disabled=true;
  $("counts_layer").replaceChildren(new Option("counts（后端默认）","counts"));updateGuidance();
}
function installChoices(data){
  state.pendingMappings=state.pendingMappings||{};
  const columns=Object.entries(data.columns||{});
  const specs={condition_col:"请选择条件列",sample_col:"不显示样本数",donor_col:"请选择供体列",batch_col:"不指定批次列"};
  for(const [key,placeholder] of Object.entries(specs)){
    const select=$(key),prior=state.pendingMappings[key]||select.value;select.replaceChildren(new Option(placeholder,""));
    for(const [name,info] of columns){const option=new Option(name,name);if(info.truncated)option.textContent=name+"（取值较多）";select.add(option);}
    if(prior){if(![...select.options].some(option=>option.value===prior)){const stale=new Option("当前文件没有此列："+prior,prior);select.add(stale);}select.value=prior;}
    state.pendingMappings[key]=select.value;
  }
  const layerSelect=$("counts_layer"),priorLayer=state.pendingMappings.counts_layer||layerSelect.value;layerSelect.replaceChildren(new Option("counts（后端默认）","counts"));
  for(const layer of data.layers||[])if(layer!=="counts")layerSelect.add(new Option(layer,layer));
  if(priorLayer&&!([...layerSelect.options].some(option=>option.value===priorLayer)))layerSelect.add(new Option("当前文件没有此 layer："+priorLayer,priorLayer));
  layerSelect.value=priorLayer||"counts";state.pendingMappings.counts_layer=layerSelect.value;
  populateConditionValues();
}
function populateConditionValues(){
  const info=state.dataset?.columns?.[$("condition_col").value];
  for(const key of ["condition_a","condition_b"]){const select=$(key),prior=Object.hasOwn(state.pendingConditions,key)?state.pendingConditions[key]:select.value;select.replaceChildren(new Option(info?"请选择条件值":"先选择条件列",""));select.disabled=!info;
    if(info)for(const value of info.values||[])select.add(new Option(value,value));
    if(prior){if(![...select.options].some(option=>option.value===prior)){select.add(new Option("当前列没有此值："+prior,prior));}select.value=prior;}
    state.pendingConditions[key]=select.value;
  }
}
function updatePurposeGuide(){
  const guide=$("purpose-guide"),info=purposeInfo[selectedPurpose()];if(!guide)return;guide.replaceChildren();
  const title=el("strong",info.title+" · "+info.tasks),line=el("p",null,"purpose-description");line.append(el("span","所需输入："),document.createTextNode(info.inputs),el("br"),el("span","预期输出："),document.createTextNode(info.outputs));guide.append(title,line);
}
function buildDesignOverview(data){
  const root=el("div",null,"design-summary-content"),o=data.overview;
  const grid=el("div",null,"design-stats");
  const metric=(label,value,detail)=>{const item=el("div",null,"design-stat");item.append(el("small",label),el("strong",value),el("small",detail||""));grid.append(item);};
  metric("细胞数 · 观测",Number(o.cells).toLocaleString(),"每个细胞不是一个生物学重复");
  metric("样本 ID",o.samples.count===null?"未映射":Number(o.samples.count).toLocaleString(),o.samples.column||"仅在明确选择样本列后统计");
  metric(o.donors.confirmed?"供体 / 生物学重复":"供体候选 · 未确认",o.donors.count===null?"待确认":Number(o.donors.count).toLocaleString(),o.donors.confirmed?o.donors.column:"后端候选「"+(o.donors.column||"无")+"」；请确认其代表独立重复");
  metric(o.conditions.confirmed?"条件组":"条件候选 · 未确认",o.conditions.column?String(o.conditions.count):"未映射",o.conditions.confirmed?o.conditions.column:"后端候选「"+(o.conditions.column||"无")+"」；需明确选择条件列与对比");
  metric(o.batches.confirmed?"批次":"批次候选 · 未确认",o.batches.count===null?"未映射":String(o.batches.count),o.batches.confirmed?o.batches.column:"后端候选「"+(o.batches.column||"无")+"」；可选");
  root.append(grid);
  if(o.conditions.groups?.length){const table=el("table",null,"design-table"),head=el("thead"),tr=el("tr");["条件","细胞数","样本 ID","供体数"].forEach(v=>tr.append(el("th",v)));head.append(tr);table.append(head);const body=el("tbody");for(const group of o.conditions.groups){const row=el("tr");row.append(el("td",group.condition),el("td",Number(group.cells).toLocaleString()),el("td",group.samples===undefined?"—":String(group.samples)),el("td",group.donors===undefined?"—":String(group.donors)));body.append(row);}table.append(body);root.append(table);}
  root.append(el("small",o.samples.message+" 唯一供体数也只是元数据计数，须结合实验设计确认。","design-disclaimer"));
  if(o.conditions.selected_a||o.conditions.selected_b)root.append(el("small","最终计划对比："+(o.conditions.selected_a||"未设置")+" / "+(o.conditions.selected_b||"未设置"),"design-disclaimer"));
  if(o.conditions.truncated)root.append(notice("条件值超过 100 个，界面只展示前 100 组。","warning"));
  if(o.scope)root.append(notice(o.scope,"warning"));
  return root;
}
function renderDesignOverview(data){
  const root=$("design-overview");root.replaceChildren(buildDesignOverview(data));
  root.classList.remove("hidden");
}
async function refreshDesignOverview(){
  if(!state.dataset||!$("data").value.trim())return;
  const request=++state.designRequest,path=$("data").value.trim();
  try{const result=await api("/api/design",formValue());if(request!==state.designRequest||path!==$("data").value.trim())return;renderDesignOverview(result);}
  catch(error){if(request===state.designRequest){const root=$("design-overview");root.replaceChildren(notice("尚未形成完整映射概览："+error.message,"warning"));root.classList.remove("hidden");}}
}
async function inspectDataset(){
  if(!$("data").value.trim())throw Error("请先选择输入文件");
  await withButton($("inspect-data"),async()=>{const requestedPath=$("data").value.trim(),data=await api("/api/dataset",{data:requestedPath});if(requestedPath!==$("data").value.trim())return;state.dataset=data;
    const root=$("dataset-info");root.replaceChildren(el("strong",Number(data.n_cells).toLocaleString()+" 个细胞 · "+Number(data.n_genes).toLocaleString()+" 个基因"));root.append(el("small","计数层："+(data.layers.join(", ")||"无独立 layer")+" · 空间坐标："+(data.spatial?"有":"无")));
    const details=el("details"),summary=el("summary","查看元数据列");details.append(summary);for(const [name,values] of Object.entries(data.columns))details.append(el("small",name+"："+(values.values||[]).join(", ")+(values.truncated?" …":"")+"（"+values.unique+" 种取值，"+values.missing+" 项缺失）"));root.append(details,el("small",data.notice));root.classList.remove("hidden");installChoices(data);updateGuidance();await refreshDesignOverview();
  });
}
function renderChecklist(prepared,root){
  const list=prepared.design?.checklist;if(!list)return;
  const section=el("section",null,"preflight-checklist"),heading=el("div",null,"preflight-heading");heading.append(el("strong","运行前检查"),el("small",prepared.plan.phase==="before_dataset_audit"?"计划预览阶段：数据审计仍将在运行时执行。":""));section.append(heading);
  const rows=[["阻塞",list.blockers||[],"error"],["警告",list.warnings||[],"warning"],["可选检查",list.optional||[],""]];
  for(const [label,items,kind] of rows){if(!items.length)continue;section.append(el("h3",label));for(const item of items){const row=el("div",null,"preflight-item");row.append(notice(item.message,kind));if(item.action)row.append(el("small","建议："+item.action));const field=$(item.field);if(field){const jump=el("button","定位并修正","text-button");jump.addEventListener("click",()=>{for(let node=field.parentElement;node&&node!==document.body;node=node.parentElement)if(node.tagName==="DETAILS")node.open=true;field.scrollIntoView({behavior:"smooth",block:"center"});field.focus();});row.append(jump);}section.append(row);}}
  root.append(section);
}
function renderPlan(prepared){
  state.preview=prepared;state.previewKey=JSON.stringify(formValue());state.designConfirmed=false;renderReuseComparison(prepared);const root=$("plan-content");root.replaceChildren();$("plan-empty").classList.add("hidden");root.classList.remove("hidden");$("plan-count").textContent=prepared.plan.tasks.length+" 个步骤";
  root.append(notice(prepared.can_run?"请核对设计摘要和计划检查；确认后才可启动。":"检查中有阻塞项，修正后重新预览。",prepared.can_run?"":"warning"));
  if(prepared.design?.overview){const section=el("section",null,"inline-design-summary");section.append(el("h3","实验设计摘要"),buildDesignOverview(prepared.design));root.append(section);}
  renderChecklist(prepared,root);
  for(const error of prepared.plan.errors||[])root.append(notice(error.message||JSON.stringify(error),"error"));if(prepared.missing_packages.length)root.append(notice("当前环境缺少依赖："+prepared.missing_packages.join(", "),"error"));
  for(const missing of prepared.plan.missing_parameters||[])root.append(notice(missing.parameter==="root_cell_id"?"未指定根细胞：运行时将省略标准拟时序及其下游分支。":missing.task_id+" 缺少参数 "+missing.parameter,"warning"));
  const list=el("div",null,"plan-tasks");prepared.plan.tasks.forEach((task,index)=>{const row=el("div",null,"plan-task");row.append(el("strong",String(index+1).padStart(2,"0")+"  "+(capabilities[task.capability]||task.capability)),el("code",task.method));if(task.target_cell_type)row.append(el("small","目标："+task.target_cell_type));list.append(row);});root.append(list);
  const confirmation=el("label",null,"design-confirmation"),check=el("input");check.type="checkbox";check.addEventListener("change",()=>{state.designConfirmed=check.checked;syncStart();updateGuidance();});confirmation.append(check,el("span","我已核对上述细胞、样本、供体、条件和批次映射，并确认此计划符合我的研究设计。"));root.append(confirmation);
  const details=el("details",null,"plan-details");details.append(el("summary","查看最终配置与完整计划"),el("pre",JSON.stringify(prepared,null,2)));root.append(details);syncStart();updateGuidance();$("plan-section").focus({preventScroll:true});$("plan-section").scrollIntoView({block:"start",behavior:"smooth"});
}
function updateGuidance(){
  const purpose=selectedPurpose(),needed=[["data","输入文件"],["study_id","研究编号"],["tissue","组织"]].filter(([id])=>!$(id).value.trim()||!$(id).validity.valid).map(([,label])=>label);
  if($("species").value==="other"&&!$("other-species").value.trim())needed.push("物种名称");$("other-species").required=$("species").value==="other";
  let overrides={};try{overrides=JSON.parse($("config_overrides").value||"{}");}catch{}
  const capability=key=>overrides?.capability_parameters?.[key]||{};
  const effective=(key,field,fallback)=>Object.hasOwn(capability(key),field)?String(capability(key)[field]??"").trim():fallback;
  const a=effective("deg","condition_a",$("condition_a").value.trim()),b=effective("deg","condition_b",$("condition_b").value.trim());
  const condition=effective("deg","condition_col",$("condition_col").value.trim()),donor=effective("deg","donor_col",$("donor_col").value.trim());
  const abundanceCondition=effective("differential_abundance","condition_col",condition),abundanceDonor=effective("differential_abundance","donor_col",donor);
  const batch=effective("integration","batch_col",$("batch_col").value.trim()),root=effective("trajectory_inference","root_cell_id",$("root_cell_id").value.trim());
  let model=$("celltypist_model_path").value.trim();const extension=overrides?.analysis_extensions?.cell_annotation;
  if(extension===false)model="";else if(extension&&typeof extension==="object"&&Object.hasOwn(extension,"model_path"))model=String(extension.model_path??"").trim();
  const annotation=effective("cell_annotation","model_path",model);
  const contrastError=Boolean(a)!==Boolean(b)?"最终配置中的条件 A 和 B 需要同时填写，或同时留空。":a&&a===b?"最终配置中的条件 A 和 B 不能相同。":"";$("condition_a").setCustomValidity(contrastError);$("condition_b").setCustomValidity(contrastError);
  const purposeMissing={deg:[!condition&&"实验条件列",!a&&"条件 A",!b&&"条件 B"],abundance:[!abundanceCondition&&"实验条件列",!a&&"条件 A",!b&&"条件 B",!abundanceDonor&&"供体列"],trajectory:[!root&&"根细胞 ID"],batch:[!batch&&"批次列"],annotation:[!annotation&&"本地 CellTypist 模型"]};
  const required=(purposeMissing[purpose]||[]).filter(Boolean);$("data-step").classList.toggle("complete",needed.length===0);$("data-step-hint").textContent=needed.length?"待填写："+needed.join("、"):state.dataset?state.dataset.n_cells.toLocaleString()+" 细胞 · "+state.dataset.n_genes.toLocaleString()+" 基因":"已填写 · 可检查数据结构";
  $("design-step").classList.toggle("complete",Boolean(a&&b&&!contrastError&&required.length===0));$("design-step-hint").textContent=contrastError?"请完善最终配置中的条件对比":required.length?"此研究问题还需："+required.join("、"):a?a+" / "+b:"请确认实验设计";
  $("method-step-hint").textContent=($("method_profile").value==="standard"?"标准分析":"基础分析")+($("advanced_analysis").checked?" · 高级统计":"");
  $("readiness-title").textContent=needed.length?"还需填写 "+needed.length+" 项":required.length?"研究目的缺少 "+required.length+" 项设计输入":contrastError?"请检查实验设计":state.preview?.can_run?(state.designConfirmed?"实验设计已确认":"请核对并确认实验设计"):"配置已填写，可以预览";
  $("readiness-hint").textContent=needed.length?needed.join(" · "):required.length?required.join(" · "):contrastError||(state.preview?.can_run?(state.designConfirmed?"确认状态仅适用于当前配置；任一输入变化后需重新预览。":"核对运行前清单与实验设计摘要后，勾选确认。") :"预览会检查最终配置、依赖与研究设计前提。");updatePurposeGuide();
}
for(const key of ["condition_col","sample_col","donor_col","batch_col","condition_a","condition_b","counts_layer"]){$(key).addEventListener("change",handle(async()=>{if(key==="condition_col"){state.pendingMappings.condition_col=$(key).value;state.pendingConditions={condition_a:"",condition_b:""};populateConditionValues();}else if(["condition_a","condition_b"].includes(key))state.pendingConditions[key]=$(key).value;else state.pendingMappings[key]=$(key).value;updateGuidance();await refreshDesignOverview();}));}
document.querySelectorAll('input[name="research_purpose"]').forEach(radio=>radio.addEventListener("change",updateGuidance));
$("study-form").addEventListener("change",event=>{if(event.target.id==="celltypist_model_path"||event.target.id==="paired"||event.target.id==="advanced_analysis"||event.target.id==="design_formula"||event.target.id==="min_donors"||event.target.id==="alpha"||event.target.id==="allow_x_as_counts")refreshDesignOverview();});
updatePurposeGuide();updateGuidance();
