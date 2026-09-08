suppressPackageStartupMessages({
  library(testthat)
  pkgload::load_all('packages/omicsCore', quiet=TRUE)
  pkgload::load_all('packages/omicsApp', quiet=TRUE)
})
Sys.setenv(OMICSAPP_DATA_DIR=tempfile('audit-projects-'))
probe <- function(label, code) {
  cat('\nPROBE:',label,'\n')
  tryCatch(force(code),error=function(e)cat('PROBE_ERROR:',conditionMessage(e),'\n'))
}
set.seed(802)
e <- matrix(rnorm(30*8,10,1),nrow=30,dimnames=list(paste0('g',1:30),paste0('s',1:8)))
e[1:12,5:8] <- e[1:12,5:8]+4
md <- data.frame(group=rep(c('A','B'),each=4),row.names=colnames(e))
fd <- data.frame(feature_id=rownames(e),feature_symbol=c('TP53','JUN','FOS','STAT1','STAT3','IL6','IFNG','CDK1','CCNB1','BCL2',paste0('SYM',11:30)),row.names=rownames(e))
x <- omics_input(e,md,fd,omics_type='proteomics',assay_type='normalized_intensity')
x <- record_step(x,'read_omics',list(path='raw/audit.xlsx',type='excel',omics_type='proteomics',assay_type='normalized_intensity'))
p <- omics_project('Project A',list(prot=x))
b <- run_diff(x,method='ttest',group_col='group',control_group='A',case_group='B')
b <- omicsApp:::stamp_from_launch(b,omicsApp:::capture_launch(p,'diff',tags='prot',params=b$params,activation=1L))
p$bundles <- list(diff=b)
p <- register_analysis(p,b)
probe('worker project still carries registry, preprocessing and history',{
  big <- p
  legacy <- b
  legacy$results$diff_object <- raw(4*1024^2)
  big$analyses[[b$provenance_binding$analysis_id]] <- legacy
  big$bundles$diff <- legacy
  big <- set_preprocessing(big,'prot',x,list(impute_method='none'))
  big$retired <- list(list(results=list(diff_result_df=b$results$diff_result_df)))
  stripped <- omicsApp:::project_for_worker(big,assays=FALSE)
  cat('EXPERIMENT_MATRIX_DIM=',paste(dim(stripped$experiments$prot$expr_mat),collapse='x'),'\n',sep='')
  cat('REGISTRY_MODEL_BYTES=',length(stripped$analyses[[1]]$results$diff_object),'\n',sep='')
  cat('PREPROCESSING_MATRIX_DIM=',paste(dim(stripped$preprocessing$prot$cleaned$expr_mat),collapse='x'),'\n',sep='')
  cat('RETIRED_COUNT=',length(stripped$retired),'\n',sep='')
  cat('WORKER_PROJECT_MIB=',future:::objectSize(stripped)/1024^2,'\n',sep='')
})
probe('report grid default arguments and stale callbacks',{
  jobs <- list()
  pr <- shiny::reactiveVal(p)
  with_mocked_bindings(
    shiny::testServer(omicsApp:::report_view_server,args=list(current_project=pr),{
      session$setInputs(run_grid=1)
      cat('GRID_JOBS=',length(jobs),'\n',sep='')
      cat('HAS_ADMISSION_GATE=',is.function(jobs[[1]]$admissible),'\n',sep='')
      worker_calls <- list()
      with_mocked_bindings(jobs[[1]]$func(),
        run_enrichment_grid=function(diff_bundle,on_progress=NULL,...) {
          worker_calls[[1]] <<- list(...)
          new_analysis_bundle('run_enrichment_grid',params=list(types='ora'),results=list(grid_log=data.frame(status='empty'),enrich_result_df=omicsCore:::new_enrich_result_template()))
        },.package='omicsCore')
      cat('GRID_EXTRA_ARGS=',paste(names(worker_calls[[1]]),collapse=','),'\n',sep='')
      pb <- omics_project('Project B',list(prot=x))
      pr(pb)
      session$flushReact()
      fake <- new_analysis_bundle('run_enrichment_grid',params=list(types='ora'),results=list(grid_log=data.frame(status='empty'),enrich_result_df=omicsCore:::new_enrich_result_template()))
      jobs[[1]]$on_success(fake)
      session$flushReact()
      cat('LATE_A_RESULT_IN_B_VIEW=',identical(grid_bundle()$provenance_binding$project_id,p$metadata$project_id),'\n',sep='')
      cat('RESULT_VALID_IN_B=',result_inputs_current(grid_bundle()$provenance_binding,pr()),'\n',sep='')
      pr(NULL)
      session$flushReact()
      cat('GRID_SURVIVES_NEW_PROJECT=',!is.null(grid_bundle()),'\n',sep='')
      jobs[[1]]$on_error('late error from Project A')
      cat('LATE_ERROR_DISPLAY=',grid_error(),'\n',sep='')
    }),
    run_async=function(func,on_success,on_error,...){
      jobs[[length(jobs)+1L]] <<- c(list(func=func,on_success=on_success,on_error=on_error),list(...))
      invisible(NULL)
    },.package='omicsApp')
})
probe('grid discards parameter provenance and exported script omits grid',{
  resolved <- list()
  grid <- with_mocked_bindings(run_enrichment_grid(b,types='ora',databases='hallmark',ora_directions='both',p_cutoff=0.001,effect_cutoff=2,output_p_cutoff=1),
    run_enrichment=function(diff_bundle,type,database,direction,...) {
      resolved <<- list(...)
      new_analysis_bundle('run_enrichment',params=resolved,results=list(enrich_result_df=omicsCore:::new_enrich_result_template()))
    },.package='omicsCore')
  cat('ACTUAL_CUTOFF=',resolved$p_cutoff,'\n',sep='')
  cat('RECORDED_CUTOFF=',if(is.null(grid$params$p_cutoff))'<missing>' else grid$params$p_cutoff,'\n',sep='')
  cat('RECORDED_EFFECT_CUTOFF=',if(is.null(grid$params$effect_cutoff))'<missing>' else grid$params$effect_cutoff,'\n',sep='')
  pg <- p
  pg$bundles$enrich_all <- grid
  script <- export_script(pg,include_plots=FALSE)
  cat('SCRIPT_HAS_GRID_CALL=',any(grepl('run_enrichment_grid',script,fixed=TRUE)),'\n',sep='')
  cat('SCRIPT_HAS_OMISSION_NOTE=',any(grepl('grid|enrich_all',script)),'\n',sep='')
  writeLines(script,'audit-results/grid-replay-generated.R')
})
probe('backend failure mislabeled empty by grid',{
  errb <- b
  errb$results$diff_result_df$adj_p_value <- 0.001
  errb$results$diff_result_df$p_value <- 0.0001
  g <- with_mocked_bindings(
    run_enrichment_grid(errb,types='ora',databases='hallmark',ora_directions='both'),
    enricher=function(...)stop('audit backend failure'),.package='clusterProfiler')
  print(g$results$grid_log)
  cat('GRID_WARNINGS=',paste(g$warnings,collapse='|'),'\n',sep='')
})
probe('workbook omits per-layer QC and integrated secondary differential',{
  pp <- p
  pp <- add_experiment(pp,'prot_second',x)
  q1 <- run_qc(x,impute_method='none',outlier_method='none')
  q2 <- run_qc(x,impute_method='min',outlier_method='none')
  pp <- set_preprocessing(pp,'prot',q1$results$cleaned_input,q1$params)
  pp <- set_preprocessing(pp,'prot_second',q2$results$cleaned_input,q2$params)
  pp$bundles$qc <- q2
  pp$analyses$secondary_diff <- b
  ss <- workbook_sheets(pp)
  cat('PREPROCESSED_LAYERS=',length(pp$preprocessing),'\n',sep='')
  cat('QC_SAMPLE_SHEETS=',sum(grepl('QC samples',names(ss),fixed=TRUE)),'\n',sep='')
  cat('DIFFERENTIAL_SHEETS=',sum(grepl('Differential',names(ss),fixed=TRUE)),'\n',sep='')
  cat('PARAMS_SECOND_LAYER_QC=',any(ss[['Analysis parameters']]$analysis=='qc: prot'),'\n',sep='')
})
cat('\nPROBES_FINISHED\n')
