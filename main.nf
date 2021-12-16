#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2021 IRB Barcelona

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  RegionalEnrichment-nf: nextflow pipeline to perfom regional enrichment analysis "
log.info "                                                                         "
log.info "-------------------------------------------------------------------------"
log.info "Copyright (C) 2021 IRB Barcelona"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-------------------------------------------------------------------------"
log.info ""

params.help = null

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '  USAGE              '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run main.nf --VCF_folder --regions_info'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --mut_folder                   FOLDER         Input folder containing txt files (one per sample, resuming all mutations).'
    log.info '    --bin_files                    FILE           File containing one line per bin path (bed file) to consider'
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --output_folder                FOLDER         Output folder containing the results (default=REA_output).'
    log.info 'Flags:'
    log.info '    --help                                        Display this message'
    log.info ''
    exit 0
}

params.mut_folder = null
params.bin_files = null
params.output_folder = "REA_output"

if(params.mut_folder == null | params.bin_files == null ){
  exit 1, "Please specify each of the following parameters: --mut_folder and --bin_files"
}

mut_files = Channel.fromPath( params.mut_folder+'/*.txt' )
bin_files=file(params.bin_files)

process reformat_vcf {

  tag {tmptag}

  input:
  file muts from mut_files
  file bin_files

  output:
  file '*_reformat*.txt' into tables

  shell:
  tmptag = muts.baseName.replace(".txt", "")
  '''
  file=!{muts}
  i=1
  while IFS= read -r bin
  do
      Rscript !{baseDir}/bin/extract_bin.R --input_table=!{muts} --bin_bed=$bin --output_table=tmp.txt
      Rscript !{baseDir}/bin/reformat_vcf.R --input_muts=tmp.txt --output_table=${file/.txt/_reformat_bin${i}.txt} --bin_bed=$bin
      let "i++"
  done < !{bin_files}
  '''

}

process REA {

  tag {table_tag}

  input:
  file table from tables

  output:
  file '_reformat.txt' into res

  shell:
  '''
  !{baseDir}/bin/REA.R
  '''

}
