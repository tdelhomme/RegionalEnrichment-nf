args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}

argsL <- as.list(do.call("cbind", parseArgs(args))[c(F,T)])
names(argsL) <- as.list(do.call("cbind", parseArgs(args))[c(T,F)])
args <- argsL;rm(argsL)

if(! is.null(args$help)) {
  cat("
      Mandatory arguments:
      --mut_type                - Mutation type (if not SNV, we do not have the nt at risk in the regression)
      Optional arguments:
      --help \n\n")
  q(save="no")
}
if(is.null(args$mut_type)) {stop("Option --mut_type should be provided")} else{mut_type=args$mut_type}

print(paste(date(), " INFO: mutation type to consider is: ", mut_type, sep=""))

library(MASS)

input_tables = list.files(path = ".", pattern = "*bin*")

i=1
for(t in input_tables){
  dattmp = read.table(t, h=T, sep="\t", stringsAsFactors = F)
  if(!is.na(dattmp[1,1])){
    if(i==1) { dat = dattmp } else { dat = rbind(dat, dattmp) } 
    i=i+1
  }
}
dat$lnNt = log(dat$ntAtRisk)
dat$bin = as.numeric(factor(dat$bin))

sink(paste(gsub("_reformat_bin1.txt", "", input_tables[1]), "_out.txt", sep=""), append = TRUE, type = "output")

if(mut_type == "SNV") { 
  res = MASS::glm.nb(formula = "count ~ bin + offset(lnNt)", data = dat) 
} else {
  print("INFO: we do not consider the nt at risk in the regression and apply simple regression")
  bins = names(table(dat$bin))
  dat2 = data.frame(count = unlist(lapply(bins, function(b) sum(dat[which(dat$bin == b), "count"]))), bin = bins)
  res = cor.test(dat$count, dat$bin,method="pearson")
}
closeAllConnections()

save(res, file = paste(gsub("SVs_", "", gsub("SNVs_", "", gsub("_reformat_bin1.txt","", input_tables[1]))), ".Rdata", sep=""))

