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

library(ggplot2)
library(ggpubr)

alldat = list.files(path = ".", pattern = ".Rdata")
celllines = unique(unlist(lapply(alldat, function(x) unlist(strsplit(x, "_"))[1])))

plot_list = list()
for(cl in celllines){
  dat_cl = alldat[which(grepl(cl , alldat))]
  allsm = unlist(lapply(dat_cl, function(x) { tmp=unlist(strsplit(gsub(".Rdata", "", x), paste(cl, "_", sep=""))); tmp[length(tmp)] } ))
  for(sm in allsm){
    load(alldat[grepl(paste(cl, "_", sm, ".Rdata", sep=""), alldat)])
    if(mut_type == "SNV"){
      val = res$coefficients["bin"]
      ddtmp = data.frame(sm = sm , coef=val, cimin = confint(res)["bin",1], cimax = confint(res)["bin",2])
    } else { # correlation test when not SNVs
      val = as.numeric(t$estimate)
      ddtmp = data.frame(sm = sm , coef=val, cimin = as.numeric(t$conf.int[1]), cimax = as.numeric(t$conf.int[2]))
    }
    if(match(sm, allsm)==1){ dd = ddtmp } else { dd = rbind(dd, ddtmp) }
  }
  assign(paste(cl, "_dat", sep=""), dd)
  dd$type = rep("CT", nrow(dd))
  dd[which(grepl("^H", dd$sm)), "type"] = "Helium"
  dd[which(grepl("^P", dd$sm)), "type"] = "Proton"

  pp <- ggplot(dd, aes(x=sm, y=coef, ymin=cimin, ymax=cimax, color=type)) +
  geom_hline(aes(yintercept=0, color=type)) +
  geom_errorbar(width=0.1) + 
  geom_point() +
  coord_flip() +
  theme_classic() + ggtitle(cl) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  plot_list[[match(cl, celllines)]] = pp
}

pdf("plot.pdf", width = 12, height = 4)
ggarrange(plotlist=plot_list, ncol=3)
dev.off()
