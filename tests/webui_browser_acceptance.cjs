"use strict";

// Optional browser regression: node tests/webui_browser_acceptance.cjs
// Use EACBP_PLAYWRIGHT_MODULE for an existing Playwright module path, and
// EACBP_BROWSER_EXECUTABLE for an existing Chromium executable when needed.
// Static assets are real; API fixtures test UI state, not scientific methods.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const http = require("node:http");
const path = require("node:path");
const {chromium} = require(process.env.EACBP_PLAYWRIGHT_MODULE || "playwright");
const {expect} = require((process.env.EACBP_PLAYWRIGHT_MODULE || "playwright")+"/test");

const staticRoot = path.resolve(__dirname, "../eacbp/webui/static");
const outputRoot = path.resolve(__dirname, "../outputs/webui_easy_acceptance");
const assets = {"/":"index.html", "/app.js":"app.js", "/results.js":"results.js", "/style.css":"style.css"};
const types = {".html":"text/html", ".js":"text/javascript", ".css":"text/css"};
const dataset = (input, cells=80) => ({path:input, n_cells:cells, n_genes:40, layers:["counts"], spatial:false,
  columns:{condition:{values:["treated","control"], unique:2, missing:0},
    donor:{values:["D1","D2","D3","D4"], unique:4, missing:0},
    sample_id:{values:["S1","S2","S3","S4"], unique:4, missing:0},
    batch:{values:["run1","run2"], unique:2, missing:0}}, notice:"Synthetic browser fixture; metadata only."});

function design(form){
  const mapped=Boolean(form.condition_col), repeated=Boolean(form.donor_col);
  const blockers=[];
  for(const field of ["condition_col","donor_col","batch_col"]){
    if(!form[field])blockers.push({field,message:"Observed metadata needs explicit confirmation: "+field,
      action:"Select the observed column after confirming its experimental meaning."});
  }
  if(mapped&&(!form.condition_a||!form.condition_b))blockers.push({field:"condition_a",
    message:"Confirm the direction of the two observed conditions.",action:"Select both compared conditions."});
  return {can_run:blockers.length===0, checklist:{blockers,warnings:[],optional:[]}, overview:{cells:80, genes:40,
    samples:{column:form.sample_col||null,count:form.sample_col?4:null,message:"样本 ID 不等于独立重复。"},
    donors:{column:form.donor_col||"donor",count:4,confirmed:repeated,paired_overlap:0,by_condition:{}},
    conditions:{column:form.condition_col||"condition",count:2,confirmed:mapped,missing:0,
      selected_a:form.condition_a||null, selected_b:form.condition_b||null,truncated:false,
      groups:[{condition:"treated",cells:40,donors:2},{condition:"control",cells:40,donors:2}]},
    batches:{column:form.batch_col||null,count:form.batch_col?2:null,confirmed:Boolean(form.batch_col),missing:0},
    scope:"合成浏览器夹具，不代表科学校验。"}};
}

function preview(form){
  const checked=design(form);
  return {form,data:form.data,manifest:{study_id:form.study_id},config:{mode:"real",method_profile:form.method_profile,
    capability_parameters:{deg:{condition_a:form.condition_a,condition_b:form.condition_b}}},
    plan:{valid:true,phase:"before_dataset_audit",tasks:[{task_id:"qc",capability:"qc",method:"basic_qc_v1"},
      {task_id:"deg",capability:"deg",method:"donor_pseudobulk_welch_v1"}],errors:[],missing_parameters:[]},
    preview_token:"synthetic-confirmation-token",missing_packages:[],design:checked,can_run:checked.can_run};
}

async function fixture(browser, viewport={width:1440,height:1000}){
  const page=await browser.newPage({viewport});
  const errors=[], requests={dataset:[],design:[],preview:[],run:[],results:[]};
  page.on("pageerror", error=>errors.push(error.message));
  const controls={dataset:null,design:null,preview:null,runs:[],run:null,results:null};
  await page.route("**/api/**", async route=>{
    const request=route.request(),url=new URL(request.url());
    const form=request.method()==="POST"?request.postDataJSON():{};
    const json=value=>route.fulfill({status:200,contentType:"application/json",body:JSON.stringify(value)});
    if(url.pathname==="/api/settings")return json({workspace:"C:/synthetic",runs_dir:"C:/synthetic/runs",
      package_dir:"C:/EACBP/eacbp",python_executable:"python",packages:{anndata:"test",scanpy:"test"},max_concurrent_jobs:1});
    if(url.pathname==="/api/runs")return json({runs:controls.runs,jobs:[]});
    if(url.pathname==="/api/files")return json({path:"C:/synthetic",parent:null,
      entries:["fixture.h5ad","counts.csv","metadata.tsv"].map(name=>({name,path:name,directory:false,size:1024})),truncated:false});
    if(url.pathname==="/api/dataset"){
      requests.dataset.push(form);if(controls.dataset){const handled=await controls.dataset(route,form,json);if(handled)return;}
      return json(dataset(form.data));
    }
    if(url.pathname==="/api/design"){
      requests.design.push(form);if(controls.design){const handled=await controls.design(route,form,json);if(handled)return;}
      return json(design(form));
    }
    if(url.pathname==="/api/preview"){
      requests.preview.push(form);if(controls.preview){const handled=await controls.preview(route,form,json);if(handled)return;}
      return json(preview(form));
    }
    if(url.pathname==="/api/run"&&request.method()==="POST"){
      requests.run.push(form);return json({job:{status:"queued",run_dir:"C:/synthetic/runs/mock"},run_id:"mock"});
    }
    if(url.pathname==="/api/run")return json(controls.run||{id:"mock",study_id:"mock",title:"Mock analysis",path:"C:/synthetic/runs/mock",
      status:"queued",events:[],log:"",report:"",summary:{},saved:{},can_resume:false,can_report:false,has_report:false});
    if(url.pathname==="/api/results"){
      requests.results.push(url.searchParams.get("id"));return json(controls.results||{tasks:[],warnings:[]});
    }
    throw new Error("Unexpected fixture endpoint: "+request.method()+" "+url.pathname);
  });
  await page.goto(baseUrl);
  await expect(page.locator("#environment")).not.toContainText("正在连接");
  return {page,requests,controls,errors};
}

async function chooseData(page,input="fixture.h5ad"){
  await page.locator("#data").fill(input);await page.locator("#data").press("Tab");
  await expect(page.locator("#dataset-info")).toBeVisible();
}
async function prepareDesign(page){
  await chooseData(page);await page.locator("#tissue").fill("PBMC");
  await page.locator('input[name="research_purpose"][value="deg"]').check();
  await page.locator("#wizard-next").click();
  await expect(page.locator('[data-setup-step="2"]')).toBeVisible();
  await page.locator("#condition_col").selectOption("condition");
  await page.locator("#donor_col").selectOption("donor");
  await page.locator("#condition_a").selectOption("treated");
  await page.locator("#condition_b").selectOption("control");
  await confirmBatch(page);
}
async function confirmBatch(page){
  const details=page.locator("#additional-mappings");
  if(!await details.evaluate(node=>node.open))await details.locator("summary").click();
  await page.locator("#batch_col").selectOption("batch");
}
async function checkNoOverflow(page){
  const dimensions=await page.evaluate(()=>({width:document.documentElement.clientWidth,scroll:document.documentElement.scrollWidth}));
  assert.ok(dimensions.scroll<=dimensions.width+1,JSON.stringify(dimensions));
}
async function checkDialogFits(page,selector){
  await checkNoOverflow(page);
  const dimensions=await page.locator(selector).evaluate(node=>({width:node.clientWidth,scroll:node.scrollWidth,
    left:node.getBoundingClientRect().left,right:node.getBoundingClientRect().right,viewport:innerWidth}));
  assert.ok(dimensions.scroll<=dimensions.width+1&&dimensions.left>=-1&&dimensions.right<=dimensions.viewport+1,JSON.stringify(dimensions));
}
async function openFixtureResults(page,controls,results){
  controls.runs=[{id:"mock",study_id:"mock",title:"Synthetic result guidance",status:"success"}];
  controls.run={id:"mock",study_id:"mock",title:"Synthetic result guidance",path:"C:/synthetic/runs/mock",
    status:"success",events:[],log:"",report:"# Synthetic report",summary:{},
    saved:{status:"success",import_completed:true},can_resume:false,can_report:true,has_report:true};
  controls.results=results;
  await page.locator("#nav-history").click();await page.locator(".run-open").click();
  await page.locator('[data-tab="results"]').click();
  await expect(page.locator("#results-reading-path")).toBeVisible();
}
async function releaseResponse(page,endpoint,release){
  const pending=page.waitForResponse(response=>new URL(response.url()).pathname===endpoint);
  release();await (await pending).finished();
  await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));
}

let baseUrl;
async function main(){
  fs.mkdirSync(outputRoot,{recursive:true});
  const server=http.createServer((request,response)=>{
    const filename=assets[new URL(request.url,"http://localhost").pathname];
    if(!filename){response.writeHead(404);response.end();return;}
    const body=fs.readFileSync(path.join(staticRoot,filename)).toString().replaceAll("__EACBP_TOKEN__","browser-fixture");
    response.setHeader("Content-Type",types[path.extname(filename)]+"; charset=utf-8");response.end(body);
  });
  await new Promise(resolve=>server.listen(0,"127.0.0.1",resolve));
  baseUrl="http://127.0.0.1:"+server.address().port;
  const browser=await chromium.launch({headless:true,...(process.env.EACBP_BROWSER_EXECUTABLE?{executablePath:process.env.EACBP_BROWSER_EXECUTABLE}:{})});
  const passed=[];
  async function run(name,action,viewport){
    const context=await fixture(browser,viewport);
    try{await action(context);assert.deepEqual(context.errors,[]);passed.push(name);console.log("PASS "+name);}
    finally{await context.page.close();}
  }
  try{
    await run("contextual beginner help explains files, design and results without leaving the wizard",async({page})=>{
      await page.locator("#beginner-open").click();await expect(page.locator("#beginner-dialog")).toBeVisible();
      await expect(page.locator('[data-help-panel="files"]')).toBeVisible();
      await expect(page.locator('[data-help-panel="files"]')).toContainText(/h5ad|H5AD/);
      await page.locator('#beginner-dialog [data-help-topic="design"]').click();
      await expect(page.locator('[data-help-panel="design"]')).toBeVisible();
      await expect(page.locator('[data-help-panel="files"]')).toBeHidden();
      await expect(page.locator('[data-help-panel="design"]')).toContainText(/独立|重复/);
      await page.locator('#beginner-dialog [data-help-topic="results"]').click();
      await expect(page.locator('[data-help-panel="results"]')).toBeVisible();
      await expect(page.locator('[data-help-panel="results"]')).toContainText(/结果|报告/);
      await page.locator("#close-beginner").click();await expect(page.locator("#beginner-dialog")).toBeHidden();
      await prepareDesign(page);await page.locator("#beginner-open").click();
      await expect(page.locator('[data-help-panel="design"]')).toBeVisible();
      await page.keyboard.press("Escape");await expect(page.locator("#beginner-dialog")).toBeHidden();
      await expect(page.locator('[data-setup-step="2"]')).toBeVisible();
    });
    await run("inline design checks link blockers to their owning field and explain comparison direction",async({page})=>{
      await chooseData(page);await page.locator("#tissue").fill("PBMC");
      await page.locator('input[name="research_purpose"][value="deg"]').check();await page.locator("#wizard-next").click();
      await expect(page.locator('#design-feedback [data-design-field="condition_col"]')).toBeVisible();
      await page.locator('#design-feedback [data-design-field="condition_col"]').click();
      await expect(page.locator("#condition_col")).toBeFocused();
      if(await page.locator("#additional-mappings").evaluate(node=>node.open))await page.locator("#additional-mappings > summary").click();
      await page.locator('#design-feedback [data-design-field="batch_col"]').click();
      assert.equal(await page.locator("#additional-mappings").evaluate(node=>node.open),true);
      await expect(page.locator("#batch_col")).toBeFocused();
      await page.locator("#condition_col").selectOption("condition");await page.locator("#donor_col").selectOption("donor");
      await page.locator("#condition_a").selectOption("treated");await page.locator("#condition_b").selectOption("control");
      await expect(page.locator("#contrast-explanation")).toContainText(/treated.*control/);
      await expect(page.locator("#contrast-explanation")).toContainText(/正|大于/);
      await expect(page.locator("#contrast-explanation")).toContainText(/负|小于/);
      await page.locator("#condition_a").selectOption("control");await page.locator("#condition_b").selectOption("treated");
      await expect(page.locator("#contrast-explanation")).toContainText(/control.*treated/);
      await confirmBatch(page);await expect(page.locator('#design-feedback [data-design-field]')).toHaveCount(0);
      await page.locator("#wizard-next").click();await expect(page.locator(".design-confirmation input")).toBeVisible();
    });
    await run("expert group overrides defer direction to the checked plan and text edits immediately clear design readiness",async({page,requests})=>{
      const form={data:"expert.h5ad",study_id:"expert_001",species:"homo_sapiens",tissue:"PBMC",research_purpose:"deg",
        condition_col:"condition",donor_col:"donor",batch_col:"batch",condition_a:"treated",condition_b:"control",
        root_cell_id:"cell_001",method_profile:"baseline",min_genes:1,max_mito_pct:20,
        config_overrides:{capability_parameters:{deg:{condition_a:"control",condition_b:"treated"}}}};
      await page.locator("#config-file").setInputFiles({name:"expert-settings.json",mimeType:"application/json",buffer:Buffer.from(JSON.stringify({webui_schema:1,form}))});
      await expect(page.locator("#dataset-info")).toBeVisible();await page.locator("#wizard-next").click();
      await expect(page.locator('[data-setup-step="2"]')).toBeVisible();
      await expect(page.locator("#condition_a")).toHaveValue("treated");await expect(page.locator("#condition_b")).toHaveValue("control");
      await expect(page.locator("#contrast-explanation")).toContainText(/专家配置.*覆盖/);
      await expect(page.locator("#contrast-explanation")).toContainText(/实际比较方向.*实际设计.*结果报告/);
      await expect(page.locator("#contrast-explanation")).not.toContainText(/为正|为负|表达变化/);
      await expect(page.locator("#design-step")).toHaveClass(/(?:^|\s)complete(?:\s|$)/);
      const checkedRequests=requests.design.length;
      await page.locator("#root_cell_id").fill("cell_002");await expect(page.locator("#root_cell_id")).toBeFocused();
      await expect(page.locator("#design-step")).not.toHaveClass(/(?:^|\s)complete(?:\s|$)/);
      assert.equal(requests.design.length,checkedRequests,"Text edits must invalidate the old ready state before blur starts a new check.");
    });
    await run("changing a field immediately removes old design feedback and ignores an older matching form",async({page,requests,controls})=>{
      await prepareDesign(page);await expect(page.locator('#design-feedback [data-design-field]')).toHaveCount(0);
      let release;
      controls.design=async(route,form,json)=>{
        if(form.donor_col!=="")return false;
        await new Promise(resolve=>release=resolve);const stale=design(form);
        stale.checklist.blockers=[{field:"donor_col",message:"STALE_DESIGN_FORM",action:"This older response must be ignored."}];
        await json(stale);return true;
      };
      const before=requests.design.length;
      await page.locator("#donor_col").selectOption("");await expect.poll(()=>requests.design.length).toBeGreaterThan(before);
      await page.locator("#donor_col").selectOption("donor");
      await expect(page.locator('#design-feedback [data-design-field]')).toHaveCount(0);
      // Return to the same donor value after a newer response. The first
      // request remains stale even when its serialized form matches again.
      controls.design=async(route,form,json)=>{
        const current=design(form);current.checklist.blockers=[{field:"donor_col",message:"CURRENT_DESIGN_FORM",action:"Latest checked mapping."}];
        await json(current);return true;
      };
      await page.locator("#donor_col").selectOption("");
      await expect(page.locator("#design-feedback")).toContainText("CURRENT_DESIGN_FORM");
      await releaseResponse(page,"/api/design",release);
      await expect(page.locator("#design-feedback")).toContainText("CURRENT_DESIGN_FORM");
      await expect(page.locator("#design-feedback")).not.toContainText("STALE_DESIGN_FORM");
      let releaseLatest;
      controls.design=async()=>{await new Promise(resolve=>releaseLatest=resolve);return false;};
      await page.locator("#donor_col").selectOption("donor");
      await expect(page.locator("#design-feedback")).not.toContainText("CURRENT_DESIGN_FORM");
      await expect.poll(()=>Boolean(releaseLatest)).toBe(true);await releaseResponse(page,"/api/design",releaseLatest);
      await expect(page.locator('#design-feedback [data-design-field]')).toHaveCount(0);
    });
    await run("late design checks cannot replace a newer inspection of the same data path",async({page,requests,controls})=>{
      let release,first=true;
      controls.design=async(route,form,json)=>{
        if(form.data!=="A.h5ad"||!first)return false;first=false;
        await new Promise(resolve=>release=resolve);const stale=design(form);
        stale.checklist.blockers=[{field:"data",message:"STALE_DESIGN_PATH",action:"This older dataset must be ignored."}];
        await json(stale);return true;
      };
      await chooseData(page,"A.h5ad");await expect.poll(()=>Boolean(release)).toBe(true);
      await chooseData(page,"B.h5ad");await chooseData(page,"A.h5ad");
      await page.locator("#tissue").fill("PBMC");await page.locator("#wizard-next").click();
      await expect(page.locator('#design-feedback [data-design-field="condition_col"]')).toBeVisible();
      assert.equal(requests.dataset.length,3);await releaseResponse(page,"/api/design",release);
      await expect(page.locator("#design-feedback")).not.toContainText("STALE_DESIGN_PATH");
      await expect(page.locator('#design-feedback [data-design-field="condition_col"]')).toBeVisible();
    });
    await run("data request can be copied and selected manually when clipboard is unavailable",async({page})=>{
      await page.locator("#beginner-open").click();
      await page.locator('[data-help-panel="files"] details').filter({has:page.locator("#metadata-request")}).locator("summary").click();
      const request=await page.locator("#metadata-request").inputValue();
      assert.match(request,/h5ad|H5AD|计数/);assert.match(request,/分组|条件/);assert.match(request,/独立|重复/);assert.match(request,/批次/);
      await page.evaluate(()=>Object.defineProperty(navigator,"clipboard",{configurable:true,value:{writeText:async value=>{window.copiedMetadataRequest=value;}}}));
      await page.locator("#metadata-request-copy").click();
      assert.equal(await page.evaluate(()=>window.copiedMetadataRequest),request);
      await expect(page.locator("#metadata-copy-status")).toContainText(/已复制/);
      await page.evaluate(()=>Object.defineProperty(navigator,"clipboard",{configurable:true,value:{writeText:async()=>{throw Error("fixture clipboard denied");}}}));
      await page.locator("#metadata-request-copy").click();await expect(page.locator("#metadata-request")).toBeFocused();
      const selection=await page.locator("#metadata-request").evaluate(node=>({start:node.selectionStart,end:node.selectionEnd,length:node.value.length}));
      assert.deepEqual(selection,{start:0,end:request.length,length:request.length});
      await expect(page.locator("#metadata-copy-status")).toContainText(/手动|选中/);
    });
    await run("result reading path distinguishes verified charts from failed tasks and jumps to the right section",async({page,controls,requests})=>{
      await openFixtureResults(page,controls,{notice:"SYNTHETIC browser fixture; no biological conclusion.",warnings:[],
        tasks:[{task_id:"qc",capability:"qc",status:"success",audited:true,method:"basic_qc_v1"},
          {task_id:"clustering",capability:"clustering",status:"failed",audited:false},
          {task_id:"deg",capability:"deg",status:"success",audited:true,method:"donor_pseudobulk_welch_v1"}],
        qc:{initial_cells:80,retained_cells:70,filtered_cells:10,source_uri:"synthetic://qc"},
        deg:[{task_id:"deg",method:"donor_pseudobulk_welch_v1",statistical_unit:"donor_pseudobulk",
          condition_a:"treated",condition_b:"control",contrast_label:"treated / control",effect_definition:"treated relative to control",
          total_genes:40,shown_genes:2,source_uri:"synthetic://deg",rows:[
            {gene:"FIXTURE_GENE_1",log2_fold_change:1,fdr_q_value:.02,p_value:.005},
            {gene:"FIXTURE_GENE_2",log2_fold_change:-.8,fdr_q_value:.2,p_value:.1}]}]});
      assert.deepEqual(requests.results,["mock"]);
      await expect(page.locator('[data-result-step="quality"] .result-guide-status')).toContainText("已有可读图表");
      await expect(page.locator('[data-result-step="composition"] .result-guide-status')).toContainText("有任务未完成或未通过");
      await expect(page.locator('[data-result-step="deg"] .result-guide-status')).toContainText("已有可读图表");
      await page.locator('[data-result-step="quality"] .result-guide-action').click();
      await expect(page.locator("#result-quality h2")).toBeFocused();
      await page.locator('[data-result-step="composition"] .result-guide-action').click();
      await expect(page.locator("#result-composition h2")).toBeFocused();
      await page.locator('[data-result-step="deg"] .result-guide-action').click();
      await expect(page.locator("#result-deg h2")).toBeFocused();
      assert.equal(await page.locator(".result-gene-terms").evaluate(node=>node.open),false);
      await page.locator(".result-gene-terms > summary").click();
      await expect(page.locator(".result-gene-terms")).toContainText(/FDR/);
      await expect(page.locator(".result-gene-terms")).toContainText(/log2FC/);
      assert.equal(await page.locator(".deg-optional-filters").evaluate(node=>node.open),false);
      await checkNoOverflow(page);await page.screenshot({path:path.join(outputRoot,"results-reading-desktop.png"),fullPage:true});
    });
    await run("missing and unverified results are explained without implying no biological difference",async({page,controls})=>{
      await openFixtureResults(page,controls,{tasks:[{task_id:"deg",capability:"deg",status:"success",audited:false}],warnings:[]});
      await expect(page.locator('[data-result-step="deg"] .result-guide-status')).toContainText("任务已完成，但尚无可展示的已校验结果");
      await expect(page.locator('[data-result-step="quality"] .result-guide-status')).toContainText("当前结果记录未包含此类任务");
      await expect(page.locator('[data-result-step="quality"] .result-guide-status')).toContainText("没有差异");
      await page.locator('[data-result-step="deg"] .result-guide-action').click();
      await expect(page.locator("#result-deg h2")).toBeFocused();
      await page.locator('[data-result-step="quality"] .result-guide-action').click();
      await expect(page.locator('[data-tab="events"]')).toHaveAttribute("aria-selected","true");
      await expect(page.locator("#events-panel")).toBeVisible();
    });
    await run("research completion stays separate from process success and renders gaps as text",async({page,controls})=>{
      await openFixtureResults(page,controls,{tasks:[],warnings:[]});
      await expect(page.locator("#detail-stats")).toContainText("未记录");
      await expect(page.locator("#detail-stats")).not.toContainText("研究要求已完成");
      controls.run.summary={analysis_request_status:"partial",scientific_gaps:[{reason:"<img src=x onerror=alert(1)> missing replication"}]};
      await page.locator("#nav-history").click();await page.locator(".run-open").click();
      await expect(page.locator("#detail-stats")).toContainText("仍有未完成要求");
      await page.locator("#detail-notice details > summary").click();
      await expect(page.locator("#detail-notice")).toContainText("<img src=x onerror=alert(1)>");
      assert.equal(await page.locator("#detail-notice img").count(),0);
      controls.run.summary={analysis_request_status:"complete",scientific_gaps:[]};
      await page.locator("#nav-history").click();await page.locator(".run-open").click();
      await expect(page.locator("#detail-stats")).toContainText("研究要求已完成");
      await expect(page.locator("#detail-notice")).not.toContainText("missing replication");
    });
    await run("mobile beginner help topics and metadata request fit within the viewport",async({page})=>{
      await page.locator("#beginner-open").click();
      for(const topic of ["files","design","results"]){
        await page.locator('#beginner-dialog [data-help-topic="'+topic+'"]').click();
        await expect(page.locator('[data-help-panel="'+topic+'"]')).toBeVisible();
        await checkDialogFits(page,"#beginner-dialog");
        if(topic==="design"){
          await page.locator('[data-help-panel="design"] details > summary').click();
          await checkDialogFits(page,"#beginner-dialog");
        }
      }
      await page.locator('#beginner-dialog [data-help-topic="files"]').click();
      await page.locator('[data-help-panel="files"] details').filter({has:page.locator("#metadata-request")}).locator("summary").click();
      await expect(page.locator("#metadata-request")).toBeVisible();await checkDialogFits(page,"#beginner-dialog");
      await page.screenshot({path:path.join(outputRoot,"beginner-help-mobile.png"),fullPage:true});
      await page.locator("#metadata-request-copy").scrollIntoViewIfNeeded();
      await expect(page.locator("#metadata-request-copy")).toBeInViewport();await checkDialogFits(page,"#beginner-dialog");
      await page.screenshot({path:path.join(outputRoot,"metadata-request-mobile.png"),fullPage:true});
      await page.locator("#close-beginner").click();await checkNoOverflow(page);
    },{width:390,height:844});
    await run("mobile result reading guide keeps task status and navigation usable",async({page,controls})=>{
      await openFixtureResults(page,controls,{warnings:[],tasks:[
        {task_id:"qc",capability:"qc",status:"success",audited:true},
        {task_id:"deg",capability:"deg",status:"pending",audited:false}],
        qc:{initial_cells:80,retained_cells:70,filtered_cells:10}});
      await expect(page.locator('[data-result-step="deg"] .result-guide-status')).toContainText("有任务尚未完成");
      await checkNoOverflow(page);await page.screenshot({path:path.join(outputRoot,"results-reading-mobile.png"),fullPage:true});
      await page.locator('[data-result-step="deg"] .result-guide-action').click();
      await expect(page.locator("#result-deg h2")).toBeFocused();await checkNoOverflow(page);
    },{width:390,height:844});
    await run("guided steps, explicit mapping, automatic preview and confirmation invalidation",async({page,requests})=>{
      await expect(page.locator('[data-setup-step="1"]')).toBeVisible();
      await expect(page.locator('[data-setup-step="2"]')).toBeHidden();
      assert.match(await page.locator("#study_id").inputValue(),/^[A-Za-z0-9][A-Za-z0-9_-]*$/);
      assert.notEqual(await page.locator("#study_id").inputValue(),"study_001");
      await page.locator("#wizard-next").click();await expect(page.locator('[data-setup-step="1"]')).toBeVisible();
      await chooseData(page);await page.locator("#tissue").fill("PBMC");
      await page.locator('[data-wizard-step="3"]').click();
      await expect(page.locator('[data-setup-step="2"]')).toBeVisible();
      assert.equal(requests.preview.length,0);
      await expect(page.locator("#condition_col")).toHaveValue("");await expect(page.locator("#donor_col")).toHaveValue("");
      assert.match(await page.locator('#condition_col option[value="condition"]').innerText(),/treated|control/);
      await page.locator("#condition_col").selectOption("condition");await expect(page.locator("#contrast-fields")).toBeVisible();
      await page.locator("#donor_col").selectOption("donor");await page.locator("#condition_a").selectOption("treated");
      await page.locator("#condition_b").selectOption("control");await confirmBatch(page);await page.locator("#wizard-next").click();
      await expect(page.locator(".design-confirmation input")).toBeVisible();assert.equal(requests.preview.length,1);
      assert.equal(await page.locator("#analysis-settings").evaluate(node=>node.open),false);
      await expect(page.locator("#start")).toBeDisabled();await page.locator(".design-confirmation input").check();
      await expect(page.locator("#start")).toBeEnabled();await checkNoOverflow(page);
      await page.locator("#analysis-settings > summary").click();await page.locator("#min_genes").fill("110");
      await expect(page.locator("#start")).toBeDisabled();await page.locator("#preview").click();
      await expect(page.locator(".design-confirmation input")).toBeVisible();await expect(page.locator(".design-confirmation input")).not.toBeChecked();
      await page.locator(".design-confirmation input").check();await page.locator("#start").click();
      await expect.poll(()=>requests.run.length).toBe(1);assert.equal(requests.run[0].min_genes,110);assert.equal(requests.run[0].preview_token,"synthetic-confirmation-token");
    });
    await run("file picker inspects automatically",async({page,requests})=>{
      await page.locator("#browse").click();await page.locator(".file-entry").filter({hasText:"fixture.h5ad"}).click();
      await expect(page.locator("#dataset-info")).toBeVisible();assert.equal(requests.dataset.length,1);
      await expect(page.locator("#data")).toHaveValue("fixture.h5ad");
    });
    await run("CSV selectors stay in the import dialog without inspecting tables as H5AD",async({page,requests})=>{
      await page.locator(".data-help > summary").click();await page.locator("#csv-import-open").click();
      await page.locator("#browse-import-counts").click();await page.locator(".file-entry").filter({hasText:"counts.csv"}).click();
      await expect(page.locator("#import-counts")).toHaveValue("counts.csv");await expect(page.locator("#import-dialog")).toBeVisible();
      await page.locator("#browse-import-metadata").click();await page.locator(".file-entry").filter({hasText:"metadata.tsv"}).click();
      await expect(page.locator("#import-metadata")).toHaveValue("metadata.tsv");await expect(page.locator("#confirm-import")).toBeDisabled();
      await expect(page.locator("#data")).toHaveValue("");assert.equal(requests.dataset.length,0);
      await page.locator("#cancel-import").click();await page.locator("#browse").click();
      await page.locator(".file-entry").filter({hasText:"counts.csv"}).click();await expect(page.locator("#import-dialog")).toBeVisible();
      assert.equal(requests.dataset.length,0);
    });
    await run("imported settings retain explicit mappings and comparison direction",async({page,requests})=>{
      const form={data:"imported.h5ad",study_id:"imported_001",species:"homo_sapiens",tissue:"PBMC",
        research_purpose:"deg",condition_col:"condition",donor_col:"donor",batch_col:"batch",condition_a:"treated",condition_b:"control",
        method_profile:"baseline",min_genes:1,max_mito_pct:20};
      await page.locator("#config-file").setInputFiles({name:"settings.json",mimeType:"application/json",buffer:Buffer.from(JSON.stringify({webui_schema:1,form}))});
      await expect(page.locator("#dataset-info")).toBeVisible();await expect(page.locator("#study_id")).toHaveValue("imported_001");
      await page.locator("#wizard-next").click();await expect(page.locator("#condition_col")).toHaveValue("condition");
      await expect(page.locator("#donor_col")).toHaveValue("donor");await expect(page.locator("#condition_a")).toHaveValue("treated");
      await expect(page.locator("#batch_col")).toHaveValue("batch");
      await expect(page.locator("#condition_b")).toHaveValue("control");await page.locator("#wizard-next").click();
      await expect(page.locator(".design-confirmation input")).toBeVisible();assert.equal(requests.preview.at(-1).method_profile,"baseline");
      assert.equal(requests.preview.at(-1).condition_a,"treated");assert.equal(requests.preview.at(-1).condition_b,"control");
    });
    await run("late inspection cannot replace a newer selection of the same path",async({page,requests,controls})=>{
      let release,occurrences=0;
      controls.dataset=async(route,form,json)=>{
        if(form.data!=="A.h5ad")return false;
        occurrences++;if(occurrences===1){await new Promise(resolve=>release=resolve);await json(dataset(form.data,12));}
        else await json(dataset(form.data,56));return true;
      };
      await page.locator("#data").fill("A.h5ad");await page.locator("#data").press("Tab");
      await expect.poll(()=>requests.dataset.length).toBe(1);
      await chooseData(page,"B.h5ad");await chooseData(page,"A.h5ad");
      await expect(page.locator("#dataset-info")).toContainText("56");await releaseResponse(page,"/api/dataset",release);
      await expect.poll(()=>page.locator("#inspect-data").isDisabled()).toBe(false);
      await expect(page.locator("#dataset-info")).toContainText("56");await expect(page.locator("#dataset-info")).not.toContainText("12 个细胞");
    });
    await run("late preview cannot revoke the newer confirmed plan",async({page,requests,controls})=>{
      let release,first=true;
      controls.preview=async(route,form,json)=>{if(!first)return false;first=false;await new Promise(resolve=>release=resolve);
        const stale=preview(form);stale.can_run=false;stale.design.checklist.blockers=[{field:"data",message:"STALE_PREVIEW"}];await json(stale);return true;};
      await prepareDesign(page);await page.locator("#wizard-next").click();await expect.poll(()=>requests.preview.length).toBe(1);
      await page.locator("#analysis-settings > summary").click();await page.locator("#min_genes").fill("101");
      await page.locator("#min_genes").fill("100");
      // Native form submission can overlap the pending HTTP request; disabled
      // refresh buttons must not make an older response authoritative again.
      await page.locator("#study-form").evaluate(node=>node.requestSubmit());
      await expect(page.locator(".design-confirmation input")).toBeVisible();await page.locator(".design-confirmation input").check();
      await expect(page.locator("#start")).toBeEnabled();await releaseResponse(page,"/api/preview",release);
      await expect.poll(()=>page.locator("#preview").isDisabled()).toBe(false);
      await expect(page.locator("#start")).toBeEnabled();await expect(page.locator("#plan-content")).not.toContainText("STALE_PREVIEW");
    });
    await run("invalid numeric controls reveal their step and collapsed settings",async({page,requests})=>{
      await prepareDesign(page);await page.locator("#statistics-settings > summary").click();
      await page.locator("#min_donors").fill("1");await page.locator("#statistics-settings > summary").click();
      await page.locator("#wizard-next").click();await expect(page.locator('[data-setup-step="2"]')).toBeVisible();
      assert.equal(await page.locator("#statistics-settings").evaluate(node=>node.open),true);
      await expect(page.locator("#min_donors")).toBeFocused();assert.equal(requests.preview.length,0);
      await page.locator("#min_donors").fill("2");await page.locator("#wizard-next").click();
      await expect(page.locator(".design-confirmation input")).toBeVisible();
      await page.locator("#analysis-settings > summary").click();await page.locator("#min_genes").fill("-1");
      await page.locator("#analysis-settings > summary").click();await page.locator("#preview").click();
      assert.equal(await page.locator("#analysis-settings").evaluate(node=>node.open),true);
      await expect(page.locator("#min_genes")).toBeFocused();assert.equal(requests.preview.length,1);
    });
    await run("unconfirmed batch blocks the plan and correction opens the owning step and details",async({page})=>{
      await chooseData(page);await page.locator("#tissue").fill("PBMC");
      await page.locator('input[name="research_purpose"][value="deg"]').check();await page.locator("#wizard-next").click();
      await page.locator("#condition_col").selectOption("condition");await page.locator("#donor_col").selectOption("donor");
      await page.locator("#condition_a").selectOption("treated");await page.locator("#condition_b").selectOption("control");
      if(await page.locator("#additional-mappings").evaluate(node=>node.open))await page.locator("#additional-mappings > summary").click();
      await page.locator("#wizard-next").click();await expect(page.locator("#start")).toBeDisabled();
      await expect(page.locator(".preflight-item .text-button")).toBeVisible();await page.locator(".preflight-item .text-button").click();
      await expect(page.locator('[data-setup-step="2"]')).toBeVisible();
      assert.equal(await page.locator("#additional-mappings").evaluate(node=>node.open),true);
      await expect(page.locator("#batch_col")).toBeFocused();await expect(page.locator("#start")).toBeDisabled();
      await page.locator("#batch_col").selectOption("batch");await page.locator("#wizard-next").click();
      await expect(page.locator(".design-confirmation input")).toBeVisible();await page.locator(".design-confirmation input").check();
      await expect(page.locator("#start")).toBeEnabled();
    });
    await run("reinspecting the same input invalidates a previously confirmed plan",async({page,requests})=>{
      await prepareDesign(page);await page.locator("#wizard-next").click();
      await expect(page.locator(".design-confirmation input")).toBeVisible();await page.locator(".design-confirmation input").check();
      await expect(page.locator("#start")).toBeEnabled();await page.locator('[data-wizard-step="1"]').click();
      await page.locator("#inspect-data").click();await expect.poll(()=>requests.dataset.length).toBe(2);
      await expect(page.locator("#inspect-data")).toBeEnabled();await expect(page.locator("#start")).toBeDisabled();
    });
    await run("desktop wizard layout and all three step screenshots",async({page})=>{
      await checkNoOverflow(page);await page.screenshot({path:path.join(outputRoot,"step1-desktop.png"),fullPage:true});
      await prepareDesign(page);await checkNoOverflow(page);await page.screenshot({path:path.join(outputRoot,"step2-desktop.png"),fullPage:true});
      await page.locator("#wizard-next").click();await expect(page.locator(".design-confirmation input")).toBeVisible();
      await checkNoOverflow(page);await page.screenshot({path:path.join(outputRoot,"step3-desktop.png"),fullPage:true});
    });
    await run("mobile wizard remains usable without horizontal overflow",async({page})=>{
      await checkNoOverflow(page);
      await page.screenshot({path:path.join(outputRoot,"step1-mobile.png"),fullPage:true});
      await prepareDesign(page);await checkNoOverflow(page);await page.screenshot({path:path.join(outputRoot,"step2-mobile.png"),fullPage:true});
      await page.locator("#wizard-next").click();
      await expect(page.locator(".design-confirmation input")).toBeVisible();await checkNoOverflow(page);
      await page.screenshot({path:path.join(outputRoot,"step3-mobile.png"),fullPage:true});
    },{width:390,height:844});
    const result={notice:"SYNTHETIC ONLY; browser UI state fixtures, no biological interpretation",passed:passed.length,tests:passed};
    fs.writeFileSync(path.join(outputRoot,"browser-acceptance.json"),JSON.stringify(result,null,2));
    console.log(JSON.stringify(result,null,2));
  }finally{await browser.close();await new Promise(resolve=>server.close(resolve));}
}
main().catch(error=>{console.error(error);process.exitCode=1;});
