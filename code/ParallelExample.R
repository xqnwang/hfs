library(magrittr)
library(doParallel)
library(foreach)
library(parallel)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
MonARCH <- FALSE
foreach(i = 1:5) %dopar% {
  if (MonARCH){
    path <- "~/.local/share/r-miniconda/envs/r-reticulate/bin/python3.8"
    setwd(Sys.glob(file.path("~/wm15/", "*", "hfs")))
  } else{
    path <- "~/Library/r-miniconda-arm64/bin/python3.10"
  }
  reticulate::use_python(path, required = T)
  reticulate::source_python("Python/example.py")
  test(i)
}
stopCluster(cl)  
