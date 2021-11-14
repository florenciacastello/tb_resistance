import subprocess as sp
import argparse
import os


def e(command):
    sp.run(command, shell=True) #ejecuta el comado en la consola.

docker='docker run --rm -w /out -v $PWD:/out -u $(id -u ${USER}):$(id -g ${USER})'


def prepare_ref(data_dir):
    cmd=f'wget -O {data_dir}/h37rv.gbk.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gbff.gz"'
    print('Descargando anotacion')
    if not os.path.exists(f'{data_dir}/h37rv.gbk'):
        e(cmd)
        print('Descargando referencia')
    cmd=f'wget -O {data_dir}/h37rv.fna.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz"'
    if not os.path.exists(f'{data_dir}/h37rv.fna'):
        e(cmd)
        print('Referencia descargada')
    print('Indexando referencia')
    e(f'gunzip {data_dir}/h37rv.fna.gz')
    e(f'gunzip {data_dir}/h37rv.gbk.gz')
    cmd=f"sed 's|LOCUS       NC_000962  |LOCUS       NC_000962.3|' {data_dir}/h37rv.gbk|sed 's|ACCESSION   NC_000962|ACCESSION   NC_000962.3|'> {data_dir}/genes.gbk"
    e(cmd)
    cmd=f'wget -O {data_dir}/snpEff.zip "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip"'
    if not os.path.exists(f'{data_dir}/snpEff'):
        e(cmd)
        e(f'unzip {data_dir}/snpEff.zip -d {data_dir}') #buscar como pasar el output a data_dir
    e(f'{docker} staphb/bwa:0.7.17 bwa index {data_dir}/h37rv.fna')
    e(f'{docker} staphb/samtools:1.13 samtools faidx {data_dir}/h37rv.fna')
    e(f'{docker} staphb/samtools:1.13 samtools dict -o {data_dir}/h37rv.dict {data_dir}/h37rv.fna')
    print('Done')
    cmd=f'mkdir -p {data_dir}/snpEff/data/h37rv'
    e(cmd)
    cmd=f'cp {data_dir}/genes.gbk {data_dir}/snpEff/data/h37rv/'
    if not os.path.exists(f'{data_dir}/snpEff/data/h37rv/genes.gbk'):
        e(cmd)
        cmd=f'echo "h37rv.genome : h37rv" >> {data_dir}/snpEff/snpEff.config'
        e(cmd)
        cmd=f'echo "  h37rv.chromosomes : NC_000962.3" >> {data_dir}/snpEff/snpEff.config'
        e(cmd)
        cmd= f'cd {data_dir}/snpEff && java -jar snpEff.jar build -genbank -v h37rv'
        e(cmd)

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-d','--data', help= 'Directorio donde se guarda la referencia y los archivos aux.', default='./data')
 #   parser.add_argument('-cpus','--cpus', help= '', default=4)
  #  parser.add_argument('-fastq1', '--fastq1', help= '', required=True)
  #  parser.add_argument('-fastq2', '--fastq2', help= '', requiere=True)
  #  parser.add_argument('-rd', '--reference_dir', help= '', default=m)
  #  parser.add_argument('-rf', '--reference_filename', help= '', required=True)
  #  parser.add_argument('-wd', '--', help= '', required=True)
    args=parser.parse_args()
    if not os.path.exists(args.data):
        os.makedirs(args.data)
    
    assert os.path.exists(args.data), f'No se pudo crear {args.data}'
    prepare_ref(args.data)


if __name__=='__main__':
    main()