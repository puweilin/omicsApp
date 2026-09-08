source('audit-results/probe_fixtures.R')
probe('integration restarts after writing its own result',{
  rn <- omics_input(e,md,fd,omics_type='rnaseq',assay_type='logcpm')
  pp <- add_experiment(p,'rnaseq',rn)
  jobs <- list()
  with_mocked_bindings(shiny::testServer(omicsApp:::app_server,{
    current_project(pp)
    activation(1L)
    session$setInputs('diff-layer'='prot','diff-method'='ttest','diff-group_col'='group','diff-control'='A','diff-case'='B')
    session$setInputs('diff-rerun'=1)
    diffs <- which(vapply(jobs,function(j)identical(j$message,'Running differential analysis...'),logical(1)))
    cat('JOB_LABELS=',paste(vapply(jobs,function(j)j$message,''),collapse='|'),'\n',sep='')
    if(length(diffs)) {
      j <- jobs[[tail(diffs,1)]]
      Sys.sleep(1.1)
      out <- j$func()
      if(is.null(j$admissible)||isTRUE(j$admissible('success')$ok))j$on_success(out)
      session$flushReact()
    }
    for(k in 1:3){
      ints <- which(vapply(jobs,function(j)identical(j$message,'Running multi-omics integration...'),logical(1)))
      cat('INTEGRATION_JOBS_BEFORE_RELEASE_',k,'=',length(ints),'\n',sep='')
      if(!length(ints))break
      j <- jobs[[tail(ints,1)]]
      Sys.sleep(1.1)
      out <- j$func()
      if(is.null(j$admissible)||isTRUE(j$admissible('success')$ok))j$on_success(out)
      session$flushReact()
      cat('INTEGRATION_PRESENT_',k,'=',!is.null(current_project()$bundles$integration),'\n',sep='')
      cat('REGISTRY_COUNT_',k,'=',length(current_project()$analyses),'\n',sep='')
    }
    ints <- which(vapply(jobs,function(j)identical(j$message,'Running multi-omics integration...'),logical(1)))
    cat('INTEGRATION_JOBS_FINAL=',length(ints),'\n',sep='')
  }),run_async=function(func,on_success,on_error,message='Running...',...){
    jobs[[length(jobs)+1L]] <<- c(list(func=func,on_success=on_success,on_error=on_error,message=message),list(...))
    invisible(NULL)
  },.package='omicsApp')
})
probe('integration adjusted-model metadata changes binding after stamping',{
  options(shiny.allowoutputreads=TRUE)
  xx <- x
  xx$meta_df$age <- c(20,30,25,35,22,32,27,37)
  rn <- omics_input(round(2^xx$expr_mat),xx$meta_df,fd,omics_type='rnaseq',assay_type='raw_count')
  pp <- omics_project('Adjusted',list(prot=xx,rnaseq=rn))
  bb <- run_diff(xx,method='limma',group_col='group',control_group='A',case_group='B',covariates='age')
  bb <- omicsApp:::stamp_from_launch(bb,omicsApp:::capture_launch(pp,'diff',tags='prot',params=bb$params,activation=1L))
  pp$bundles <- list(diff=bb)
  pp <- register_analysis(pp,bb)
  shiny::testServer(omicsApp:::integration_view_server,args=list(current_project=shiny::reactiveVal(pp),diff_bundle=shiny::reactiveVal(bb),activation=shiny::reactiveVal(1L)),{
    session$setInputs(rerun=1)
    z <- integration_bundle()
    cat('INTEGRATION_ERROR=',integration_error(),'\n',sep='')
    cat('INTEGRATION_COMPUTED=',!is.null(z),'\n',sep='')
    cat('SIGNATURE_EQUALS_ACTUAL_PARAMS=',identical(z$provenance_binding$params_signature,params_signature(z$params)),'\n',sep='')
    q <- register_analysis(pp,secondary_bundle())
    q <- register_analysis(q,z)
    q$bundles$integration <- z
    print(validate_analysis_registry(q))
    q <- retire_invalid(q)
    cat('ADJUSTED_INTEGRATION_RETAINED=',!is.null(q$bundles$integration),'\n',sep='')
  })
})
cat('\nFOURTH_PROBES_FINISHED\n')
