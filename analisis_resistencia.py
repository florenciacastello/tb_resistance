#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 15:06:59 2022

@author: fleer
"""
import argparse
import pandas as pd
from collections import defaultdict
import sys
import math
from Bio.SeqUtils import seq1
import json
import subprocess as sp
import argparse
import os



def e(command):
    if os.environ.get('verbose', False):
        sys.stderr.write(command+'\n')
    sp.run(command, shell=True) #ejecuta el comado en la consola.

DOCKER='docker run --rm -w /out -v $PWD:/out -u $(id -u ${USER}):$(id -g ${USER})'


''' la parte de andy en proces reference y el resto aca'''

def resistance_notcovered(data_dir, out_dir, sample, resistance_db):
    docker=f'{DOCKER} -v {out_dir}:{out_dir} -v {data_dir}:{data_dir}'
    cmd=f'{docker} staphb/bedtools:2.30.0 bedtools intersect -a {data_dir}/{resistance_db} -wa \
    -b {out_dir}/{sample}_final_merged.bed > {out_dir}/{sample}_resistance_notcovered.bed' 
    e(cmd)


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-d','--data', help= 'Directorio donde se guarda la \
                        referencia y los archivos aux.', default='./data')
    parser.add_argument('-o','--out', help= 'Directorio donde se guarda los \
                        archivos de salida', default='./out_data3')
    parser.add_argument('-sa', '--sample', help= '', default=None)
    parser.add_argument('-rdb','--res_db', help= 'Base de datos donde se guardan\
                        las resistencias', default='andytb.bed')
    parser.add_argument('--verbose', help= '', action='store_true')
    args=parser.parse_args()
    os.environ['verbose']=str(args.verbose)
    if not os.path.exists(args.data):
        os.makedirs(args.data)
    assert os.path.exists(args.data), f'No se pudo crear {args.data}'
    if not args.sample:
        print(f'{args.sample} does not exist')
    resistance_notcovered(os.path.abspath(args.data), os.path.abspath(args.out), args.sample, args.res_db)
    
    


if __name__=='__main__':
    main()