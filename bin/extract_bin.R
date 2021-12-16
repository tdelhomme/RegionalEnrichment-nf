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
      --input_table               - Input tab-del file containing mutations, with columns sample, chr, start, end, ref, alt
      --bin_bed                   - Bed file for that bin
      --output_table              - Name of the output table
      Optional arguments:
      --help \n\n")
  q(save="no")
}
if(is.null(args$input_table)) {stop("Option --input_table should be provided")} else{input_table=args$input_table}
if(is.null(args$bin_bed)) {stop("Option --bin_bed should be provided")} else{bin_bed=args$bin_bed}
if(is.null(args$output_table)) {stop("Option --output_table should be provided")} else{output_table=args$output_table}
print(paste("INFO: input table is: ", input_table))
print(paste("INFO: bin bed file is: ", bin_bed))

library(GenomicRanges)
library(MutationalPatterns)
library(liftOver)

# chain file for the liftover
path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain") # if not present, download here: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
ch = import.chain(path)

df_bed = read.table(bin_bed) ; names(df_bed) = c("chr", "start", "end")
gr19 = makeGRangesFromDataFrame(df = df_bed)
# liftover
seqlevelsStyle(gr19) = "UCSC"  # necessary
grbins = unlist(liftOver(gr19, ch)) # liftover returns a list of Granges objects

df = read.table(input_table, h=T, sep="\t")
grmuts = makeGRangesFromDataFrame(df)

mutsbin = GenomicRanges::intersect(grmuts, grbins)

res = data.frame(mutsbin)
res$ref = df[match(start(mutsbin), df$start), "ref"]
res$alt = df[match(start(mutsbin), df$start), "alt"]
colnames(res)[which(colnames(res) == "seqnames")] = "chr"

write.table(res, file = output_table, quote = FALSE, row.names = F, sep = "\t")
