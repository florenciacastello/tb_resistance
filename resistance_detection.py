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
from collections import defaultdict
import sys
from Bio.SeqUtils import seq1
import json
from tqdm import tqdm


def resistance_detector(resistance_tb_path):
    """ resistance tb columns: 
    GeneID/GeneName/Drug/NucleotidePosH37/NucleotidePosGene/REF/ALT/AApos/AAref
    AAalt/Effect/MultipleMutation/PhylogeneticMarker/Source"""
    genes_truncados=defaultdict(list)
    genes_sinaa=defaultdict(list)
    resistance_tb=pd.read_csv(resistance_tb_path) 
    
    resistance_tb=resistance_tb[resistance_tb['Drug']!='-']
    
    def splitpositions(x) :
        if "-" not in str(x):
            return int(str(x).split('/')[0])   
        else:
            return "-"
    resistance_tb['NucleotidePosH37']=[  splitpositions(x)
                                       for x in resistance_tb.NucleotidePosH37
                                       ]
    resistance_tb['NucleotidePosGene']=[  splitpositions(x)
                                       for x in resistance_tb.NucleotidePosGene
                                       ]
    #resistance_tb['AApos']=resistance_tb['AApos'].replace('-',None).astype(int)
    #resistance_tb['NucleotidePosH37']=resistance_tb['NucleotidePosH37'].replace('-',None)
    #resistance_tb['NucleotidePosGene']=resistance_tb['NucleotidePosGene'].replace('-',None).astype(int)
    for i, r in resistance_tb.iterrows():
        if r.AAalt in ['STOP', 'fs']:
            entry={'row':r, 'Drug':r.Drug}
            genes_truncados[r.GeneID].append(entry)
        elif r.AAalt=='-':
            entry={'row':r, 'Drug':r.Drug}
            genes_sinaa[r.NucleotidePosH37].append(entry)
    resistance_tb['json']=[x.to_dict() for _, x in resistance_tb.iterrows()]
    return resistance_tb, genes_truncados, genes_sinaa

def process_line(resistance_tb, genes_truncados, genes_sinaa, samples, line):
    vec = line.split()
    pos = int(vec[1])
    ref = vec[3]
    alts = {i: x for i, x in enumerate(vec[4].split(","))}
    
    gt_options = {k + 1: v for k, v in alts.items()}
    gt_options[0] = ref
    #gt_options[99] = "N"                
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
    variant_final=[]
    for index, sample_name in enumerate(samples, 9):
        sample=vec[index].split(':')
        gt_str=sample[gt_index].replace('|', '/')
        gt=gt_str.split('/')
        ad=sample[ad_index].split(',')
        dp=sample[dp_index]
        if (gt[0]!='.' and gt_str!='0/0'):
            try:
                
                '''asumimos que el segundo es el ALT pero hay que actualzarlo 
                para procesar ambos'''
                alt=gt_options[int(gt[1])]
                info=vec[7]
                ann=info.split('ANN=')[1].split(';')[0].split(',')[0].split('|')
                variant_type=ann[1]
                impact=ann[2]
                gene_name=ann[3]
                lt=ann[4]
                gene_variant=ann[9][2:]
                if ann[11]:
                    posgen=ann[11].split('/')[0] 
                else:
                    posgen=0
                if ann[10]:
                    protein_variant=ann[10][2:]
                    posaa=ann[13].split('/')[0]
                    #posaa=str(math.ceil(float(posgen)/3))
                    #print('----', posaa)
                    variant_reference=['missense_variant', 'stop_gained', 
                                       'frameshift_variant', 'start_lost', 
                                       'synonymous_variant']
                    if set(variant_type.split('&')) & set(variant_reference):
                        aa_ref=protein_variant.split(posaa)[0]
                        if posaa in protein_variant:
                            aa_alt=protein_variant.split(posaa)[1]
                        elif protein_variant.endswith('fs'):
                            aa_alt='fs'
                        elif protein_variant.endswith('*'):
                            aa_alt='*'
                        else:
                            posaa=None
                        
                    else:
                        posaa=None
                        aa_ref=None
                        aa_alt=None
                else:
                    protein_variant=None
                    posaa=None
                    aa_ref=None
                    aa_alt=None
                    
                
                variante_ann=procesar_variante(pos, ref, posgen, variant_type, impact,  lt, posaa, aa_ref, aa_alt, alt,
                                  resistance_tb, genes_truncados, genes_sinaa, gene_name)
                if variante_ann: 
                    variante_ann['AD']={allele:int(ad[index]) for index, allele in gt_options.items()}
                    variante_ann['ID']=sample_name
                    variante_ann['variant_type']=variant_type
                    variante_ann['gene_name']=gene_name
                    variante_ann['Pos']=pos
                    variante_ann['Pos_gen']=posgen
                    variante_ann['Ref']=ref
                    variante_ann['Alt']=alt
                    variante_ann['aa_pos']=posaa
                    variante_ann['aa_ref']=aa_ref
                    variante_ann['aa_alt']=aa_alt
                    '''1-4, 8, 9  ANN=T|missense_variant|MODERATE|dnaA|Rv0001|transcript|Rv0001|protein_coding|1/1|c.71C>T|p.Pro24Leu|71/1524|71/1524|24/507||'''
                    variant_final.append(variante_ann)
                    #print(variante_ann)
            except:
                print('Error')
                print(line)
                print(gt_options, ad)
                print(ann)
                raise
    return variant_final

def process_vcf(vcf, resistance_tb_path, out):
     resistance_tb, genes_truncados, genes_sinaa=resistance_detector(resistance_tb_path)
     with open(vcf) as h, open (out, "w") as hw:
         first=True
         hw.write('[')
         for line in tqdm(h):
             line=line.strip()
             if not line.startswith("##"):
                if line.startswith("#"):
                    samples=line.split()[9:]
                else: 
                    try:
                        record=process_line(resistance_tb, genes_truncados, genes_sinaa, samples, line)
                        if record:
                            for r in record:
                                if not first:
                                    hw.write(',')
                                    #pass
                                else:
                                    first=False
                                json.dump(r, hw, indent=4)
                    except:
                        print(line)
                        raise
         hw.write(']')
         #print(genes_truncados)
                        
                    
def procesar_variante(pos, ref, posgen, variant_type, impact, lt, posaa, aa_ref, 
                      aa_alt, alt, resistance_tb, genes_truncados, genes_sinaa, gene_name):
    ''' 1) variante tal cual en genoma (conf alta)
        2) variante tal cual en prot (conf alta)
        3) variante tal cual en el gen (conf alta)
        4) frameshift/STOP de alto impacto en gen truncado (conf alta)
        5) variante en pos reportada (conf media)
        6) si esta dentro de un gen de interes y es de impacto moderado (del/ins/NOsinonim)/alto (rumor)
        7) si esta dentro de un gen con impacto bajo (sinonima/mod codon no expr) (rumor)
        8) si estan cerca de un gen con bajo/moderado impacto (rumor)'''
    
    ''' 1) variante tal cual en genoma (conf alta)'''
    df_resistance_tb=resistance_tb[(resistance_tb['NucleotidePosH37']==int(pos))&
                                    (resistance_tb['REF']==ref)&
                                    (resistance_tb['ALT']==alt)]
    if len(df_resistance_tb):
        return {'confidence':'high', 'Drugs':list(set(df_resistance_tb['Drug'])), 'db_entries':
                list(df_resistance_tb['json'])}
   
    """2) variante tal cual en prot (conf alta)
    pensar mejor npara el caso de delec e inserciones de aa"""
    if aa_alt:        
        aa_ref=seq1(aa_ref)
        aa_alt=seq1(aa_alt)
        df_resistance_tb=resistance_tb[(resistance_tb['GeneName']==gene_name) &
                                        (resistance_tb['AApos']==int(posaa)) &
                                        (resistance_tb['AAref']==aa_ref)&
                                        (resistance_tb['AAalt']==aa_alt)]
        if len(df_resistance_tb):
            return {'confidence':'high', 'Drugs':list(set(df_resistance_tb['Drug'])), 'db_entries':
                    list(df_resistance_tb['json'])}
    
    
    '''3) variante tal cual en el gen (conf alta)'''
    df_resistance_tb=resistance_tb[(resistance_tb['GeneName']==gene_name) &
                                    (resistance_tb['NucleotidePosGene']==int(posgen))&
                                    (resistance_tb['REF']==ref)&
                                    (resistance_tb['ALT']==alt)]
    if len(df_resistance_tb):
        return {'confidence':'high', 'Drugs':list(set(df_resistance_tb['Drug'])), 'db_entries':
                list(df_resistance_tb['json'])}
    
    '''4) frameshift/STOP de alto impacto en gen truncado (conf alta)'''
    trunc=['*', 'fs']
    if aa_alt in trunc: # agregar star_lost 
        if lt in genes_truncados:
            return {'confidence':'high', 'Drugs':list(set([x['Drug'] for x in genes_truncados[lt]])), 'db_entries':
                    [r for r in genes_truncados[lt]]}    
                
    '''5) variante en pos reportada (conf media)'''
    df_resistance_tb=resistance_tb[(resistance_tb['NucleotidePosH37']==int(pos))&
                                    (resistance_tb['GeneID']==lt)]
    if len(df_resistance_tb):
        if 'synonymous' not in variant_type:
            return {'confidence':'medium', 'Drugs':list(set(df_resistance_tb['Drug'])), 'db_entries':
                    list(df_resistance_tb['json'])}    
    df_resistance_tb=resistance_tb[(resistance_tb['NucleotidePosGene']==int(posgen))&
                                    (resistance_tb['GeneID']==lt)]
    if len(df_resistance_tb):
        if 'synonymous' not in variant_type:
            return {'confidence':'medium', 'Drugs':list(set(df_resistance_tb['Drug'])), 'db_entries':
                    list(df_resistance_tb['json'])}    
    
    '''6) si esta dentro de un gen de interes y es de impacto moderado (del/ins/NOsinonim)/alto (rumor)'''
    rumor_impact=['MODERATE', 'HIGH']
    if impact in rumor_impact:
        if lt in set(resistance_tb['GeneID']):
            df_resistance_tb=resistance_tb[resistance_tb['GeneID']==lt]
            return {'confidence':'low', 'Drugs':list(set(df_resistance_tb['Drug'])), 'db_entries':
                    list(df_resistance_tb['json'])}  
            
    '''7) si esta dentro de un gen con impacto bajo (sinonima/mod codon no expr) '''
    if 'LOW' in impact:
        if lt in set(resistance_tb['GeneID']):
            df_resistance_tb=resistance_tb[resistance_tb['GeneID']==lt]
            return {'confidence':'low', 'Drugs':list(set(df_resistance_tb['Drug'])), 'db_entries':
                    list(df_resistance_tb['json'])}  

    '''8) si estan cerca de un gen con bajo/moderado impacto (rumor)'''
    rumor_impact=['MODIFIER', 'LOW']
    if 'upstream_gene_variant' in variant_type:
        if impact in rumor_impact:
            if lt in set(resistance_tb['GeneID']):
                df_resistance_tb=resistance_tb[resistance_tb['GeneID']==lt]
                return {'confidence':'low', 'Drugs':list(set(df_resistance_tb['Drug'])), 'db_entries':
                            list(df_resistance_tb['json'])} 

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-d','--data', help= 'Directorio donde se guarda la \
                        referencia y los archivos aux.', default='./data')
    parser.add_argument('-o','--out', help= 'Directorio donde se guarda los \
                        archivos de salida', default='./out_resistance.json')
    parser.add_argument('-sa', '--sample', help= '', default=None)
    parser.add_argument('-va', '--vcf_ann', help= '', default='./out_data/variant.ann.vcf')
    parser.add_argument('-r', '--reference', help= '', default='andytb.csv')
    
    args=parser.parse_args()
    vcf_procesado=process_vcf(args.vcf_ann, args.reference, args.out)
    
    


if __name__=='__main__':
    main()