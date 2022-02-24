# tb_resistance



## Pasos: 

Debe estar Git, python 3 y Docker instalado en la PC. (https://docs.docker.com/engine/install/ubuntu/)
y deben estar hechos los pasos post instalacion:
(https://docs.docker.com/engine/install/linux-postinstall/)

 	$ sudo usermod -aG docker $USER  ##posiblemente haya que cerrar todas las terminales luego de correr este comando y volver a abrirlas. Para chequear que este todo bien correr:
	
	$ docker run hello-world       



por ultimo correr:
	
	$ pip install --user biopython 

0) Descarga de scripts
	
	$ git clone https://github.com/florenciacastello/tb_resistance.git
  	$ cd tb_resistance/ #esta sera la carpeta de trabajo
   
### El proceso consta de 5 scripts:
1) Prepare_reference.py (descarga y prepara la referencia de H37Rv del NCBI)

	**$ python3 prepare_reference.py**

	Esto por defecto genera una carpeta ./data/ en la carpeta de trabajo.
	
2) Process_sample.py (Prepara las muestras tomando fastq1 de R1 y R2. Para más de 1 puede meterse en bash. Para que este punto funcione los archivos deben estar en la misma carpeta tb_resistance/ o algun subdirectorio.)

	**$ python3 process_sample.py -fastq1 muestra1_R1.fastq -fastq2 muestra1_R2.fastq**
 
 	 Repetir este paso con CADA muestra. (la carpeta de salida por defecto es ./vcfs en ella deposita al muestra1.vcf para cada muestra.)
	 
3) Process_group.py (toma la carpeta de ./vcfs y los une en un archivo “variant.ann.vcf”. Este ultimo es el archivo ya procesado que contiene a todas las muestras y ya se encuentra anotadas las variantes.)

	**$ python3 process_group.py -i vcfs/**

4) Resistance_detection.py (toma el “variant.ann.vcf” y lo compara contra las referencias y bases originales buscando resistencias. Tiene como salida un .json ‘out_resistance’. actualmente solo trabaja con la primer muestra, en proceso para multiples)

	**$ python3 resistance_detection.py**



5) Analisis_resistencia.py (en proceso…

------------------------------------------------------------------------------------------------------------------------------------------
**Borrador para los scripts**
input: .fastq pair end, lista de resistencias fenotipicas. 
output: .vcf (anotado), resistencia, comparativa

--->vcf multimuestra usando como ref H37Rv y anotado con SNPEff. Da la diferencia de la ref contra las variantes.
--->resistencia: mutaciones asociadas a una resistencia. 
--->comparativa: diferencias entre MUESTRAS de mi proy. 2 tablas: mutaciones asociadas a la resistencia, falsos pos+neg.


etapa1: procesar cada una de las muestras. 
-->2 scipts:
	*procesa una cepa de punta a punta. (input: fastq1/fastq2. output: .vcf de haplotipos, .bam mapeo)
	*Toma un directorio de muestras (.fastq) y le aplica el script anterior a todos.

etapa2: analizar el grupo. comparativa. 
--> 
	*realiza el vcf multimuestra (junta todos los archivos en uno).
	*obtener la resistencia.
	*obtener reportes.
	
