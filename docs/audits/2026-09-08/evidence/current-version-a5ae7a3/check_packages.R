Sys.setenv(NOT_CRAN='true', '_R_CHECK_FORCE_SUGGESTS_'='false')
audit_lib <- normalizePath('audit-results',mustWork=TRUE)
audit_lib <- file.path(audit_lib,'library')
dir.create(audit_lib,showWarnings=FALSE)
.libPaths(c(audit_lib,.libPaths()))
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
install.packages('packages/omicsCore',repos=NULL,type='source',lib=audit_lib,quiet=TRUE)
for(p in c('omicsCore','omicsApp')) {
  cat('\nCHECK_PACKAGE:',p,'\n')
  r <- rcmdcheck::rcmdcheck(file.path('packages',p),
    args=c('--no-manual','--no-tests','--no-vignettes'),error_on='never',
    check_dir=file.path('audit-results',paste0(p,'-check')))
  saveRDS(r,file.path('audit-results',paste0(p,'-check.rds')))
  print(r)
  cat('CHECK_COUNTS:',p,'errors=',length(r$errors),'warnings=',length(r$warnings),'notes=',length(r$notes),'\n')
}
