#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 17:57:46 2021

@author: fleer
"""


import subprocess as sp
import argparse
import os
import sys
from glob import glob


def e(command):
    if os.environ.get('verbose', False):
        sys.stderr.write(command+'\n')
    sp.run(command, shell=True) #ejecuta el comado en la consola.

DOCKER='docker run --rm -w /out -v $PWD:/out -u $(id -u ${USER}):$(id -g ${USER})'

#{out_dir}/raw.vcf.gz
def procces_samples(data_dir, output_dir, vcf_dir, reference_filename):
    vcf_files=[]
    
    for x in glob(vcf_dir + "/*vcf*"):
        if x.endswith(".vcf") or x.endswith(".vcf.gz") or x.endswith(".gvcf") or x.endswith(".gvcf.gz"):
            vcf_files.append(x)
    if not vcf_files:
        raise FileNotFoundError(f'no .vcf or .vcf.gz files where found at {vcf_dir}')
    
    vcfs = " ".join([f"--variant {x}" for x in vcf_files])
    cmdx = f"""docker run -u $(id -u ${{USER}}):$(id -g ${{USER}})  \
        -v $PWD:/out  -v {data_dir}:{data_dir} -v {vcf_dir}:{vcf_dir} -v {output_dir}:{output_dir}\
        -w /out broadinstitute/gatk:4.2.2.0 \
        gatk CombineGVCFs -R {data_dir}/{reference_filename} {vcfs} -O {output_dir}/combined.raw.vcf.gz"""
    e(cmdx)
    cmdx = f"""docker run -u $(id -u ${{USER}}):$(id -g ${{USER}})  \
        -v $PWD:/out -v {data_dir}:{data_dir}\
        -v {output_dir}:{output_dir}\
        -w /out broadinstitute/gatk:4.2.2.0 gatk GenotypeGVCFs \
                    -R "{data_dir}/{reference_filename}" -ploidy 2 \
                    -V "{output_dir}/combined.raw.vcf.gz" \
                    -O "{output_dir}/variant.vcf" """
    e(cmdx)
    cmdx=f'java -jar {data_dir}/snpEff/snpEff.jar -stats {output_dir}/snpEff.html \
        h37rv {output_dir}/variant.vcf > {output_dir}/variant.ann.vcf 2>> {output_dir}/snpEff.log'
    e(cmdx)
 


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-d','--data', help= 'Directorio donde se guarda la \
                         referencia y los archivos aux.', default='./data')
    parser.add_argument('-o','--out', help= 'Directorio donde se guarda los \
                         archivos de salida', default='./out_data')
    parser.add_argument('-i','--vcfs_dir', help= 'Directorio donde se encuentran \
                        los vcf crudos', required=True)
    parser.add_argument('-rf', '--reference_filename', help= '', default='h37rv.fna')
    parser.add_argument('--verbose', help= '', action='store_true')
    args=parser.parse_args()
    # os.environ['verbose']=str(args.verbose)
    if not os.path.exists(args.out):
        os.makedirs(args.out)
       
    
    assert os.path.exists(args.out), f'No se pudo crear {args.out}'
    assert os.path.exists(args.data), f'No existe el directorio {args.data}'
    output_dir=os.path.abspath(args.out)
    data_dir=os.path.abspath(args.data)
    vcf_dir=os.path.abspath(args.vcfs_dir)
    print(f'Running vcf_dir:{vcf_dir}\n data_dir:{data_dir}\n output_dir:{output_dir}\n')
    process=procces_samples(data_dir, output_dir, vcf_dir, args.reference_filename)


if __name__=='__main__':
    main()