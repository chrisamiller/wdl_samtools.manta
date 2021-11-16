version 1.0

#######################
## SAMTOOLS WDL STEP ##
#######################
task samtools {
  input {
    File tumor_bam
    File tumor_bam_bai
    String output_filename_tumor
    File fusion_sites
  }
  Int cores = 4
  Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
  Float regions_size = size([fusion_sites], "GB")
  Int size_needed_gb = 10 + round(bam_size + regions_size)
  runtime {
    memory: "16GB"
    cpu: cores
    docker: "chrisamiller/docker-genomic-analysis:latest"
    disks: "local-disk ~{size_needed_gb} SSD"
  }
  command <<<
    set -o pipefail
    set -o errexit
    #Working through the Tumor BAM
    samtools view -H ~{tumor_bam} > tmp.sam
    cat ~{fusion_sites} | while read chr start stop ; do samtools view ~{tumor_bam} $chr:$start-$stop >> input.sam ; done
    sort -S 8G input.sam | uniq >> tmp.sam
    samtools sort -O bam -o ~{output_filename_tumor} tmp.sam
    samtools index ~{output_filename_tumor}
    rm tmp.sam
    rm input.sam
  >>>
  output {
    File sorted_bam_tumor = "tumor.filtered.sorted.bam"
    File sorted_bam_tumor_bai = "tumor.filtered.sorted.bam.bai"
  }
}




####################
## MANTA WDL STEP ##
####################
task manta {
  input {
    File tumor_bam
    File tumor_bam_bai
    File manta_config
    File reference
    File reference_fai
    File reference_dict
    Boolean non_wgs = false
    Boolean output_contigs = false
  }
  Int cores = 1
  Float ref_size = size([reference, reference_fai, reference_dict], "GB")
  Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
  Int size_needed_gb = 10 + round(ref_size + bam_size)
  runtime {
    docker: "mgibio/manta_somatic-cwl:1.6.0"
    cpu: cores
    ram: "16GB"
    bootDiskSizeGb: 10
    disks: "local-disk ~{size_needed_gb} SSD"
  }
  command <<<
    /usr/bin/python /usr/bin/manta/bin/configManta.py \
    ~{if non_wgs then "--exome" else ""} \
    ~{if output_contigs then "--outputContig" else ""} \
    --config ~{manta_config} \
    --referenceFasta ~{reference} \
    --tumorBam ~{tumor_bam} \
    --runDir $PWD \
    && /usr/bin/python runWorkflow.py -m local -g 8 \
    -j ~{cores}
  >>>
  output {
    File? diploid_variants = "results/variants/diploidSV.vcf.gz"
    File? diploid_variants_tbi = "results/variants/diploidSV.vcf.gz.tbi"
    File? somatic_variants = "results/variants/somaticSV.vcf.gz"
    File? somatic_variants_tbi = "results/variants/somaticSV.vcf.gz.tbi"
    File all_candidates = "results/variants/candidateSV.vcf.gz"
    File all_candidates_tbi = "results/variants/candidateSV.vcf.gz.tbi"
    File small_candidates = "results/variants/candidateSmallIndels.vcf.gz"
    File small_candidates_tbi = "results/variants/candidateSmallIndels.vcf.gz.tbi"
    File? tumor_only_variants = "results/variants/tumorSV.vcf.gz"
    File? tumor_only_variants_tbi = "results/variants/tumorSV.vcf.gz.tbi"
  }
}




##################
## WDL WORKFLOW ##
##################
workflow wf {
  input {
    File tumor_bam
    File tumor_bam_bai
    String output_filename_tumor
    File fusion_sites
    File manta_config
    File reference
    File reference_fai
    File reference_dict
    Boolean manta_non_wgs = true
    Boolean? manta_output_contigs
  }
  call samtools {
    input:
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    output_filename_tumor=output_filename_tumor,
    fusion_sites=fusion_sites
  }
  call manta {
    input:
    tumor_bam=samtools.sorted_bam_tumor,
    tumor_bam_bai=samtools.sorted_bam_tumor_bai,
    manta_config=manta_config,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    non_wgs=manta_non_wgs,
    output_contigs=manta_output_contigs
  }
}
