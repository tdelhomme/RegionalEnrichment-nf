library(MASS)

input_tables = list.files(path = ".", pattern = "*txt")

for(t in input_tables){
  if(t == input_tables[[1]]) { dat = read.table(t, h=T, sep="\t", stringsAsFactors = F) } else {
    dat = rbind(dat,  read.table(t, h=T, sep="\t", stringsAsFactors = F))
  }
}
dat$lnNt = log(dat$ntAtRisk)
dat$bin = as.factor(dat$bin)

sink(paste(gsub("_reformat_bin1.txt", "", input_tables[1]), "_out.txt", sep=""), append = TRUE, type = "output")

MASS::glm.nb(formula = "count ~ bin + offset(lnNt)", data = dat)
closeAllConnections()
