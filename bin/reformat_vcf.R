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
      --input_muts                - Input tab-del file containing mutations, with columns sample, chr, start, end, ref, alt
      --bin_bed                   - Bed file for that bin
      --output_table              - Name of the output table
      Optional arguments:
      --help \n\n")
  q(save="no")
}
if(is.null(args$input_muts)) {stop("Option --input_muts should be provided")} else{input_muts=args$input_muts}
if(is.null(args$bin_bed)) {stop("Option --bin_bed should be provided")} else{bin_bed=args$bin_bed}
if(is.null(args$output_table)) {stop("Option --output_table should be provided")} else{output_table=args$output_table}

library(helperMut)
library(MutationalPatterns)
library(GenomicRanges)
library(Biostrings)
library(liftOver)
library("BSgenome.Hsapiens.UCSC.hg38")

# chain file for the liftover
path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain") # if not present, download here: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
ch = import.chain(path)

# mutations
muts = read.table(input_muts, sep="\t", h=T, stringsAsFactors = F)

# nucleotides at risk
genome = BSgenome.Hsapiens.UCSC.hg38
df_bed = read.table(bin_bed) ; names(df_bed) = c("chr", "start", "end")
gr19 = makeGRangesFromDataFrame(df = df_bed)
k = 3 # length of the nucleotide context
gr19 = helperMut::extend(x = gr19, upstream = k, downstream = k)

# liftover and get bins counts
seqlevelsStyle(gr19) = "UCSC"  # necessary
gr = unlist(liftOver(gr19, ch)) # liftover returns a list of Granges objects
genome_sql = seqlengths(genome)
seqlengths(gr) = suppressWarnings(genome_sql[seqlevels(gr)])
gr = GenomicRanges::trim(gr)
reg = Biostrings::getSeq(genome,gr)
oligo_counts = Biostrings::oligonucleotideFrequency(x = reg, width = k)
oligo_counts = apply(oligo_counts,2,sum)

# observed mutation types
df = read.table(input_muts, h=T, sep="\t")
grmuts = makeGRangesFromDataFrame(df)
nut3_context = MutationalPatterns::type_context(grmuts, "BSgenome.Hsapiens.UCSC.hg38") # get nut3 context for mutations PASS in the sample, SNVs or indels
all_mut96 = paste(substr(nut3_context$context,1,1), df$ref, substr(nut3_context$context,3,3), ">", df$alt, sep = "") # 96 nut

df_muts = data.frame(MS = names(table(all_mut96)), count = as.numeric(table(all_mut96)), 
                     ntAtRisk = oligo_counts[unlist(lapply(names(table(all_mut96)), function(x) unlist(strsplit(x, ">"))[1]))],
                     bin = basename(bin_bed))

write.table(df_muts, file = output_table, quote = FALSE, row.names = F, sep = "\t")
