#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import subprocess as sp
import argparse
import os
import sys


def e(command):
    if os.environ.get('verbose', False):
        sys.stderr.write(command+'\n')
    sp.run(command, shell=True) #ejecuta el comado en la consola.

DOCKER='docker run --rm -w /out -v $PWD:/out -u $(id -u ${USER}):$(id -g ${USER})'


def process_sample(fastq1, fastq2, data_dir, out_dir, sample, species, reference_filename='h37rv.fna',
                   cpus=4, crop=300, headcrop=18, slidingwindow=10, windowqa=5):
    trimming_reads(fastq1, fastq2, out_dir, cpus, crop, headcrop, slidingwindow, windowqa)
    mapping(sample, species, reference_filename, data_dir, out_dir, cpus)
    haplotype(out_dir, data_dir, reference_filename)
            
def trimming_reads(fastq1, fastq2, out_dir, cpus=4, crop=300, 
                   headcrop=18, slidingwindow=10, windowqa=5):
    docker=f'{DOCKER} -v {out_dir}:{out_dir} '
    cmd=f"{docker} -v {fastq1}:{fastq1} -v {fastq2}:{fastq2} staphb/trimmomatic:0.38 \
        java -jar /Trimmomatic-0.38/trimmomatic-0.38.jar PE \
        -threads {cpus} {fastq1} {fastq2} {out_dir}/1.fastq.gz  \
        {out_dir}/1U.fastq.gz {out_dir}/2.fastq.gz {out_dir}/2U.fastq.gz \
	ILLUMINACLIP:/Trimmomatic-0.38/adapters/TruSeq2-PE.fa:2:30:10 \
        ILLUMINACLIP:/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:10 \
	ILLUMINACLIP:/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 \
        ILLUMINACLIP:/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10  \
	CROP:{crop}  HEADCROP:{headcrop}  TRAILING:3 SLIDINGWINDOW:{slidingwindow}:{windowqa} MINLEN:36"
#    print(cmd)
    e(cmd)
    
    
def mapping(sample, species, reference_filename, data_dir, out_dir, cpus=4):
    docker=f'{DOCKER} -v {out_dir}:{out_dir} -v {data_dir}:{data_dir} '
    cmd=f"{docker} staphb/bwa:0.7.17 bwa mem -o {out_dir}/mapping.sam -t {cpus} -M \
        -R '@RG\\tID:group1\\tSM:{sample}\\tPL:illumina\\tLB:{species}' \
        {data_dir}/{reference_filename} {out_dir}/1.fastq.gz {out_dir}/2.fastq.gz "
    e(cmd)
    e(f'{docker} staphb/samtools:1.13 bash -c "samtools view -b -F 4 \
      {out_dir}/mapping.sam | samtools sort -  > {out_dir}/mapped_reads_raw.bam"')
    e(f'{docker} broadinstitute/gatk:4.2.2.0 gatk MarkDuplicates \
      -INPUT {out_dir}/mapped_reads_raw.bam -OUTPUT {out_dir}/dedup.bam -METRICS_FILE {out_dir}/metrics.txt')
    e(f'{docker} staphb/samtools:1.13 samtools sort {out_dir}/dedup.bam > {out_dir}/final.bam')
    e(f'{docker} staphb/samtools:1.13 samtools index {out_dir}/final.bam')    
    e(f'{docker} staphb/samtools:1.13 samtools flagstat {out_dir}/final.bam > {out_dir}/flagstat.txt')
    
def haplotype(out_dir, data_dir, reference_filename):
    docker=f'{DOCKER} -v {out_dir}:{out_dir} -v {data_dir}:{data_dir}'
    cmd=f"{docker} broadinstitute/gatk:4.2.2.0 gatk HaplotypeCaller -ERC GVCF  \
        -R {data_dir}/{reference_filename} -ploidy 2 -I {out_dir}/final.bam \
        --output-mode EMIT_ALL_ACTIVE_SITES -O {out_dir}/raw.vcf.gz"
    e(cmd)
    

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-d','--data', help= 'Directorio donde se guarda la \
                        referencia y los archivos aux.', default='./data')
    parser.add_argument('-o','--out', help= 'Directorio donde se guarda los \
                        archivos de salida', default='./out_data')
    parser.add_argument('-cpus','--cpus', help= '', default=4)
    parser.add_argument('-fastq1', '--fastq1', help= '', required=True)
    parser.add_argument('-fastq2', '--fastq2', help= '', required=True)
    parser.add_argument('-sa', '--sample', help= '', default=None)
    parser.add_argument('-sp', '--species', help= '', default='H37Rv')
    parser.add_argument('-rf', '--reference_filename', help= '', default='h37rv.fna')
    parser.add_argument('--verbose', help= '', action='store_true')
    args=parser.parse_args()
    os.environ['verbose']=str(args.verbose)
    if not os.path.exists(args.data):
        os.makedirs(args.data)
        
    
    assert os.path.exists(args.data), f'No se pudo crear {args.data}'
    assert os.path.exists(args.fastq1), f'No existe el archivo {args.fastq1}'
    assert os.path.exists(args.fastq2), f'No existe el archivo {args.fastq2}'
    if not args.sample:
        args.sample=args.fastq1.split('_')[0]
    process_sample(os.path.abspath(args.fastq1), os.path.abspath(args.fastq2), 
                   os.path.abspath(args.data), os.path.abspath(args.out), args.sample, 
                   args.species, args.reference_filename, args.cpus)
    


if __name__=='__main__':
    main()