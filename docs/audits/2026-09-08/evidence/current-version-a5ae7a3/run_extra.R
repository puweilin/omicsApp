suppressPackageStartupMessages(library(testthat))
Sys.setenv(NOT_CRAN='true', OMICSAPP_DATA_DIR=tempfile('audit-extra-store-'),
  OMICSAPP_FOLLICLE_ROOT='/Users/puweilin/SciProject/CHISSS',
  OMICSAPP_PROTEOMICS_XLSX='/Users/puweilin/SciProject/CHISSS/data/SkinProteomics/Proteomics_Data.xlsx',
  OMICSCORE_PERF_TESTS='1')
for(p in c('omicsCore','omicsApp')) {
 filter <- if(p=='omicsCore') 'golden-rnaseq-follicle|perf-budget' else 'golden-raw-to-result'
 cat('\nEXTRA_PACKAGE:',p,'\n')
 r <- test_local(file.path('packages',p),filter=filter,reporter='summary',stop_on_failure=FALSE)
 saveRDS(r,file.path('audit-results',paste0(p,'-extra.rds')))
 d <- as.data.frame(r)
 write.csv(d[,!vapply(d,is.list,logical(1)),drop=FALSE],file.path('audit-results',paste0(p,'-extra.csv')),row.names=FALSE)
 print(colSums(d[,intersect(c('nb','failed','skipped','error','warning','passed'),names(d)),drop=FALSE]))
}
