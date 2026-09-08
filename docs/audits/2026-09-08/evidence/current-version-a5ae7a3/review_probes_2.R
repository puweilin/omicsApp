source('audit-results/probe_fixtures.R')
probe('import scale override is missing from replay provenance',{
  tmp <- tempfile(fileext='.xlsx')
  mat <- round(e)
  openxlsx::write.xlsx(list(expression=data.frame(feature_id=rownames(mat),mat,check.names=FALSE),
    metadata=data.frame(sample_id=rownames(md),md,row.names=NULL),
    feature_annot=fd),file=tmp)
  shiny::testServer(omicsApp:::import_view_server,{
    session$setInputs(omics_type='rnaseq',file=list(name='rounded-logcpm.xlsx',datapath=tmp,size=file.info(tmp)$size,type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'))
    cat('INFERRED_SCALE=',parsed()$input$assay_type,'\n',sep='')
    session$setInputs(assay_type='logcpm',confirm=1)
    got <- confirmed_input()
    cat('CONFIRMED_SCALE=',got$assay_type,'\n',sep='')
    read_step <- Filter(function(z)identical(z$op,'read_omics'),provenance_steps(got))[[1]]
    cat('READ_PROVENANCE_SCALE=',read_step$args$assay_type,'\n',sep='')
    pp <- omics_project('Scale override',list(rnaseq=got))
    bb <- run_diff(got,method='limma',group_col='group',control_group='A',case_group='B')
    pp$bundles <- list(diff=bb)
    lines <- export_script(pp,include_plots=FALSE)
    writeLines(lines,'audit-results/scale-override-replay.R')
    replayed <- do.call(read_omics,c(list(path=tmp),read_step$args[setdiff(names(read_step$args),'path')]))$input
    cat('REPLAY_INPUT_SCALE=',replayed$assay_type,'\n',sep='')
    result <- tryCatch({run_diff(replayed,method='limma',group_col='group',control_group='A',case_group='B');'success'},error=function(e)conditionMessage(e))
    cat('REPLAY_DIFF=',result,'\n',sep='')
  })
})
probe('same-modality layers resolve the wrong input in script',{
  y <- x
  y$expr_mat[1:10,5:8] <- y$expr_mat[1:10,5:8]-8
  pp <- omics_project('Two proteomics layers',list(skin=x,follicle=y))
  bb <- run_diff(y,method='limma',group_col='group',control_group='A',case_group='B')
  bb <- omicsApp:::stamp_from_launch(bb,omicsApp:::capture_launch(pp,'diff',tags='follicle',params=bb$params,activation=1L))
  pp$bundles <- list(diff=bb)
  pp <- register_analysis(pp,bb)
  cat('TRUE_LAYER=',bb$input_info$layer_tag,'\n',sep='')
  cat('SCRIPT_RESOLVES_LAYER=',omicsCore:::resolve_tag(pp,bb),'\n',sep='')
  lines <- export_script(pp,include_plots=FALSE)
  writeLines(lines,'audit-results/wrong-layer-replay.R')
  i <- grep('diff <- run_diff',lines,fixed=TRUE)
  cat('REPLAY_CALL=',paste(lines[i+0:2],collapse=' '),'\n',sep='')
  wrong <- run_diff(x,method='limma',group_col='group',control_group='A',case_group='B')
  cat('MAX_EFFECT_DIFFERENCE=',max(abs(bb$results$diff_result_df$effect-wrong$results$diff_result_df$effect)),'\n',sep='')
})
probe('clearing a project leaves the previous import staged',{
  pr <- shiny::reactiveVal(p)
  tmp <- tempfile(fileext='.xlsx')
  openxlsx::write.xlsx(list(expression=data.frame(feature_id=rownames(e),e,check.names=FALSE),
    metadata=data.frame(sample_id=rownames(md),md,row.names=NULL),feature_annot=fd),file=tmp)
  shiny::testServer(omicsApp:::import_view_server,args=list(current_project=pr),{
    session$setInputs(omics_type='proteomics',file=list(name='old.xlsx',datapath=tmp,size=file.info(tmp)$size,type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'),normalize=FALSE)
    session$setInputs(confirm=1)
    cat('CONFIRMED_BEFORE_RESET=',!is.null(confirmed_input()),'\n',sep='')
    pr(NULL)
    session$flushReact()
    cat('CONFIRMED_AFTER_RESET=',!is.null(confirmed_input()),'\n',sep='')
    cat('PARSED_AFTER_RESET=',!is.null(parsed()$input),'\n',sep='')
  })
})
cat('\nSECOND_PROBES_FINISHED\n')
