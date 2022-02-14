library(MASS)

input_tables = list.files(path = ".", pattern = "*txt")

for(t in input_tables){
  i=1
  dattmp = read.table(t, h=T, sep="\t", stringsAsFactors = F)
  if(!is.na(dattmp[1,1])){
    if(i==1) { dat = dattmp } else { dat = rbind(dat, dattmp) } 
    i=i+1
  }
}
dat$lnNt = log(dat$ntAtRisk)
dat$bin = order(dat$bin)

sink(paste(gsub("_reformat_bin1.txt", "", input_tables[1]), "_out.txt", sep=""), append = TRUE, type = "output")

res = MASS::glm.nb(formula = "count ~ bin + offset(lnNt)", data = dat)
closeAllConnections()

save(res, file = paste(gsub("SVs_", "", gsub("SNVs_", "", gsub("_reformat_bin1.txt","", input_tables[1]))), ".Rdata", sep=""))