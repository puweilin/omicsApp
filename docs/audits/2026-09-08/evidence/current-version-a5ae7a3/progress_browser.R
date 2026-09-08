Sys.setenv(NOT_CRAN='true')
suppressPackageStartupMessages(library(shinytest2))
appdir <- normalizePath('audit-results/progress-app')
a <- AppDriver$new(appdir,load_timeout=60000,name='audit-progress')
on.exit(a$stop(),add=TRUE)
a$run_js("document.querySelector('#go').click()")
for(i in 1:9){
 Sys.sleep(1)
 cat('TICK',i,'PROGRESS=',a$get_js("document.body.innerText + ' BARS=' + Array.from(document.querySelectorAll('.progress-bar')).map(x=>x.style.width).join('|')"),'\n')
}
print(a$get_logs())
a$stop()
