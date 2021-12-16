library(MASS)

input_tables = list.files(path = ".", pattern = "*txt")

for(t in input_tables){
  if(t == input_tables[[1]]) { dat = read.table(t, h=T, sep="\t", stringsAsFactors = F) } else {
    dat = rbind(dat,  read.table(t, h=T, sep="\t", stringsAsFactors = F))
  }
}
dat$lnNt = log(dat$ntAtRisk)
dat$bin = as.factor(dat$bin)

MASS::glm.nb(formula = "count ~ bin + lnNt", data = dat)
