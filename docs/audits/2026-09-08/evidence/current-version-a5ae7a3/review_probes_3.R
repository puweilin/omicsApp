source('audit-results/probe_fixtures.R')
probe('new project cannot reuse the staged identical import',{
  tmp <- tempfile(fileext='.xlsx')
  openxlsx::write.xlsx(list(expression=data.frame(feature_id=rownames(e),e,check.names=FALSE),
    metadata=data.frame(sample_id=rownames(md),md,row.names=NULL),feature_annot=fd),file=tmp)
  shiny::testServer(omicsApp:::app_server,{
    session$setInputs('import-omics_type'='proteomics',
      'import-file'=list(name='old.xlsx',datapath=tmp,size=file.info(tmp)$size,type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'),
      'import-normalize'=FALSE)
    session$setInputs('import-confirm'=1)
    cat('PROJECT_EXISTS_AFTER_IMPORT=',!is.null(current_project()),'\n',sep='')
    session$setInputs('project-confirm_new_project'=1)
    cat('PROJECT_EMPTY_AFTER_NEW=',is.null(current_project()),'\n',sep='')
    session$setInputs('import-confirm'=2)
    cat('PROJECT_STILL_EMPTY_AFTER_RECONFIRM=',is.null(current_project()),'\n',sep='')
    cat('IMPORT_MODULE_STILL_CONFIRMED=',!is.null(imported_input()),'\n',sep='')
  })
})
probe('QC preview reruns on unrelated project result writes',{
  pr <- shiny::reactiveVal(p)
  calls <- 0L
  real_qc <- omicsCore::run_qc
  with_mocked_bindings(
    shiny::testServer(omicsApp:::qc_view_server,args=list(current_project=pr),{
      session$flushReact()
      before <- calls
      q <- pr()
      q$metadata$audit_unrelated_change <- 'new differential completed'
      pr(q)
      session$flushReact()
      cat('QC_RUNS_FOR_UNRELATED_UPDATE=',calls-before,'\n',sep='')
    }),
    run_qc=function(...){calls <<- calls+1L;real_qc(...)},.package='omicsCore')
})
cat('\nTHIRD_PROBES_FINISHED\n')
