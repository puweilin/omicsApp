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
