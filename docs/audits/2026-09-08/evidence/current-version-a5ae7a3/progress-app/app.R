pkgload::load_all('../../packages/omicsCore', quiet=TRUE)
pkgload::load_all('../../packages/omicsApp', quiet=TRUE)
future::plan(future::multisession, workers=2)
shiny::shinyApp(ui=shiny::fluidPage(shiny::actionButton('go','Go'),shiny::textOutput('status')),
server=function(input,output,session){
  cat('ALLOWOUTPUTREADS=',getOption('shiny.allowoutputreads'), '\n'); st <- shiny::reactiveVal('idle')
  output$status <- shiny::renderText(st())
  shiny::observeEvent(input$go,{
    path <- tempfile()
    read_progress <- omicsApp:::progress_reader(path)
    omicsApp:::run_async(
      omicsApp:::detached_call(function(){
        for(i in 1:6){cat(sprintf('%d\t6\tstep%d\n',i-1,i),file=path,append=TRUE);Sys.sleep(1)}
        unlink(path)
        'done'
      },path=path),
      on_success=function(res)st(res),on_error=function(msg)st(msg),
      progress_poll=function(){
        cat(sprintf('POLL_AT %s\n',format(Sys.time(),'%H:%M:%OS3')))
        read_progress()
      },message='audit progress')
  })
})
