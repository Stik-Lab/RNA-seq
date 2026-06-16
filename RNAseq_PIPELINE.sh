#!/bin/bash

#SBATCH --job-name=RNA-seq
#SBATCH --mem=50gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=RNA-seq_%A-%a.log
#SBATCH --array=1-N

# ========== VARIABLES ==========
names=( X ) # Replace X with actual sample names
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

path_fq='path_to_fastq_files'
path_bam="path_to_bam_files"
path_bw="path_to_bigwig_files"
genomeIndex='path_to_STARindexgenome'
GTFfile='path_to_annotation_file'
effective_genome_size=2913022398

BAM=${path_bam}/STAR_${describer}Aligned.sortedByCoord.out.bam
BW=${path_bw}/${describer}.bw

for dir in "${path_bam}" "${path_bw}" ; do
  if [ ! -d "${dir}" ]; then
    mkdir -p "${dir}"
  fi
done


# ========== MODULES ==========
module load fastqc-0.11.9-gcc-11.2.0-dd2vd2m
module load STAR/2.7.6a-GCC-11.2.0
module load samtools-1.12-gcc-11.2.0-n7fo7p2
module load deepTools/3.5.1-foss-2021b

# ========== STEP 1: FASTQC ==========

echo "................................................................ 1. START_FASTQC ${describer} ................................................................"
FQC_DONE=$(ls ${path_fq}/${describer}*_fastqc.zip 2>/dev/null | head -n 1)
if [ -n "${FQC_DONE}" ]; then
    echo "SKIP FastQC: reports already exist for ${describer}"
else
  fastqc ${path_fq}/${describer}*.f*q.gz -o ${path_fq}
fi

echo "................................................................ 1. END_FASTQC ${describer} ................................................................"

# ========== STEP 2: ALIGNMENT ==========

echo "................................................................ 3. START_ALIGNMENT ${describer} ................................................................"

if [ -s "${BAM}" ] && [ -s "${BAM}.bai" ]; then
    echo "SKIP STAR: BAM already exists for ${describer}"
else
  STAR --genomeDir ${genomeIndex} \
      --sjdbGTFfile ${GTFfile} \
      --readFilesCommand zcat  \
      --readFilesIn ${path_fq}/${describer}_*1.f*q.gz ${path_fq}/${describer}_*2.f*q.gz \
      --runThreadN 4 \
      --outFileNamePrefix ${path_bam}/STAR_${describer} \
      --outSAMtype BAM SortedByCoordinate --outTmpDir tmp_${describer} \
      --outWigStrand Unstranded \
      --quantMode GeneCounts

  samtools index ${path_bam}/STAR_${describer}Aligned.sortedByCoord.out.bam
fi

echo "................................................................ 2. END_ALIGNMENT ${describer} ................................................................"


# ========== STEP 3: BAM TO BIGWIG =========

echo "................................................................ 9. START_BAMCOVERAGE ${describer} ................................................................"

if [ -s "${BW}" ]; then
    echo "SKIP bamCoverage: bigwig already exists for ${describer}"
else

  bamCoverage --bam ${path_bam}/STAR_${describer}Aligned.sortedByCoord.out.bam \
        --outFileName ${path_bw}/${describer}.bw \
        --effectiveGenomeSize 2913022398 \
        --outFileFormat bigwig --binSize 1 \
        --normalizeUsing RPGC
fi

echo "................................................................ 9. END_BAMCOVERAGE ${describer} ................................................................"


