source('audit-results/probe_fixtures.R')
options(shiny.allowoutputreads=TRUE)
probe('same-modality integration uses second layer twice',{
  y <- x
  y$expr_mat[1:10,5:8] <- y$expr_mat[1:10,5:8]-8
  pp <- omics_project('Two proteomics layers',list(skin=x,follicle=y))
  bb <- run_diff(y,method='limma',group_col='group',control_group='A',case_group='B')
  bb <- omicsApp:::stamp_from_launch(bb,omicsApp:::capture_launch(pp,'diff',tags='follicle',params=bb$params,activation=1L))
  pp$bundles <- list(diff=bb)
  pp <- register_analysis(pp,bb)
  shiny::testServer(omicsApp:::integration_view_server,args=list(current_project=shiny::reactiveVal(pp),diff_bundle=shiny::reactiveVal(bb),activation=shiny::reactiveVal(1L)),{
    session$setInputs(rerun=1)
    cat('TRUE_PRIMARY_LAYER=',bb$input_info$layer_tag,'\n',sep='')
    cat('CHOSEN_PRIMARY_LAYER=',can_run()$primary_tag,'\n',sep='')
    cat('SECONDARY_LAYER=',can_run()$secondary_tag,'\n',sep='')
    cat('INTEGRATION_ERROR=',integration_error(),'\n',sep='')
    z <- integration_bundle()$results$integration_df
    cat('INTEGRATION_ROWS=',nrow(z),'\n',sep='')
    cat('MAX_REPORTED_EFFECT_DIFFERENCE=',max(abs(z$effect_a-z$effect_b)),'\n',sep='')
  })
})
cat('\nFIFTH_PROBES_FINISHED\n')
