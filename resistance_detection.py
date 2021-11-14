#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
primero leer el csv (pandas)
procesar el csv. cuales son los genes que afectados dan resistencia
leer posicion por posicion y ver 1)si la posicion esta en ungen de interes
        2) si esta cerca de un gen de interes  
        3)si la variante esta tal cual reportada
        4)si el cambio es distinto pero la posicion es la misma
        5) si rompe la proteina (stop, start o frameship)

reportar cuando esta en baja proporcion

va anecesitar andytb(pandas), y un vcf.ann y un directorio de out
argparse
buscar como leer vcf 
'''
import argparse
import pandas as pd
import defaultdict
import sys

'''
import vcf
def read(variant):
    reader = vcf.Reader(open(variant))
    df = pd.DataFrame([vars(r) for r in reader])
    out = df.merge(pd.DataFrame(df.INFO.tolist()), left_index=True, right_index=True)
    return out
o vcfpy
'''

def resistance_detector(resistance_tb_path):
    """ resistance tb columns: 
    GeneID/GeneName/Drug/NucleotidePosH37/NucleotidePosGene/REF/ALT/AApos/AAref
    AAalt/Effect/MultipleMutation/PhylogeneticMarker/Source"""
    genes_truncados=defaultdict(list)
    genes_sinaa=defaultdict(list)
    resistance_tb=pd.read_csv(resistance_tb_path)
    for i, r in resistance_tb.iterrows():
        if r.AAalt in ['STOP', 'fs']:
            entry={'row':r, 'Drug':r.Drug}
            genes_truncados[r.GeneID].append(entry)
        elif r.AAalt=='-':
            entry={'row':r, 'Drug':r.Drug}
            genes_sinaa[r.NucleotidePosH37].append(entry)
    return resistance_tb, genes_truncados, genes_sinaa

def process_vcf(vcf, resistance_tb_path):
     resistance_tb, genes_truncados, genes_sinaa=resistance_detector(resistance_tb_path)
     with open(vcf) as h:
         for line in h:
             if not line.startswith("#"):             
                vec = line.split()
                pos = int(vec[1])
                ref = vec[3]
                alts = {i: x for i, x in enumerate(vec[4].split(","))}
                
                gt_options = {k + 1: v for k, v in alts.items()}
                gt_options[0] = ref
                gt_options[99] = "N"                
                format_f = vec[8]  # GT:AD:DP:GQ:P
                gt_index = [i for i, x in enumerate(format_f.split(":")) if x == "GT"]
                if not gt_index:
                    sys.stderr.write("not all pos/samples have a GT field")
                    sys.exit(2)
                gt_index = gt_index[0]
                ad_index = [i for i, x in enumerate(format_f.split(":")) if x == "AD"]
                ad_index = ad_index[0]
                dp_index = [i for i, x in enumerate(format_f.split(":")) if x == "DP"]
                dp_index = dp_index[0]
                low_freq = False
                sample=vec[9].split(':')
                gt=sample[gt_index].replace('|', '/').split('/')
                ad=sample[ad_index].split(',')
                dp=sample[dp_index]
                if gt[0]!='.':
                    '''asumimos que el segundo es el ALT pero hay que actualzarlo 
                    para procesar ambos'''
                    alt=gt_options[int(gt[1])]
                    variante_ann=procesar_variante(pos, posgen, gen, posaa, aa, alt, \
                                      resistance_tb, genes_truncados, genes_sinaa)
                

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-d','--data', help= 'Directorio donde se guarda la \
                        referencia y los archivos aux.', default='./data')
    parser.add_argument('-o','--out', help= 'Directorio donde se guarda los \
                        archivos de salida', default='./out_resistance')
    parser.add_argument('-sa', '--sample', help= '', default=None)
    parser.add_argument('-va', '--vcf_ann', help= '', default='variant.ann.vcf')
    parser.add_argument('-r', '--reference', help= '', default='andytb.csv')
    args=parser.parse_args()
    


if __name__=='__main__':
    main()