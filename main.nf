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
    log.info '    --VCF_folder                   FOLDER         Input folder containing VCF files to process (extension .vcf.gz).'
    log.info '    --bin_files                    FILE           File containing one line per bin path (bed file) to consider'
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --output_folder                FOLDER         Output folder containing the results (default=REA_output).'
    log.info 'Flags:'
    log.info '    --help                                        Display this message'
    log.info ''
    exit 0
}

params.VCF_folder = null
params.bin_files = null
params.output_folder = "REA_output"

if(params.VCF_folder == null | params.bin_files == null ){
  exit 1, "Please specify each of the following parameters: --VCF_folder and --bin_files"
}

vcfs = Channel.fromPath( params.VCF_folder+'/*.vcf.gz' )
bin_files=file(params.bin_files)

process splitVCF {

  tag {vcf_tag}

  input:
  file vcf from vcfs

  output:
  file '*_split.vcf.gz' into sm_vcf

  shell:
  vcf_tag = vcf.baseName.replace(".vcf","")
  '''
  for sample in `bcftools query -l !{vcf}`; do
   bcftools view -c1 -Oz -s $sample -o ${sample}_split.vcf.gz !{vcf}
  done
  '''

}

process reformat_vcf {

  tag {vcf_tag}

  input:
  file vcf from sm_vcf
  file bin_files

  output:
  file '_reformat.txt' into tables

  shell:
  vcf_tag = vcf.replace(".vcf.gz","")
  '''
  file=!{vcf}
  while IFS= read -r bin
  do
      bcftools -R $bin $file | bgzip -c > tmp.vcf.gz
      !{baseDir}/bin/reformat_vcf.R --input_vcf=tmp.vcf.gz --output_table=${file/.vcf*/_reformat.txt} --bin_bed=$bin
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
  table_tag = table.baseName
  '''
  !{baseDir}/bin/REA.R --input_table=!{table}
  '''

}
