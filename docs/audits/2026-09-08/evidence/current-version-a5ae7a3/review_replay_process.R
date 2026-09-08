source('audit-results/probe_fixtures.R')
probe('exported scale-override script executed in a fresh R process', {
  store <- tempfile('replay-store-')
  dir.create(store)
  Sys.setenv(OMICSAPP_DATA_DIR=store)
  tmp <- tempfile(fileext='.xlsx')
  mat <- round(e)
  openxlsx::write.xlsx(list(
    expression=data.frame(feature_id=rownames(mat),mat,check.names=FALSE),
    metadata=data.frame(sample_id=rownames(md),md,row.names=NULL),
    feature_annot=fd),file=tmp)
  shiny::testServer(omicsApp:::import_view_server, {
    session$setInputs(omics_type='rnaseq', file=list(
      name='rounded-logcpm.xlsx',datapath=tmp,size=file.info(tmp)$size,
      type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'))
    session$setInputs(assay_type='logcpm',confirm=1)
    got <- confirmed_input()
    pp <- omics_project('Scale override',list(rnaseq=got))
    pp$bundles <- list(diff=run_diff(got,method='limma',
      group_col='group',control_group='A',case_group='B'))
    script <- file.path(store,'replay.R')
    writeLines(export_script(pp,include_plots=FALSE),script)
    lib <- normalizePath('audit-results/library')
    res <- callr::r(function(script,store,lib) {
      .libPaths(c(lib,.libPaths()))
      setwd(store)
      tryCatch({source(script);'success'},error=function(e)conditionMessage(e))
    },args=list(script=script,store=store,lib=lib))
    cat('ORIGINAL_DIFF_SUCCEEDED=',!is.null(pp$bundles$diff),'\n',sep='')
    cat('FRESH_PROCESS_REPLAY=',res,'\n',sep='')
  })
})
