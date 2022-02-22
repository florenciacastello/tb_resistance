# tb_resistance



## Pasos: 

Debe estar Docker instalado en la CPU. (https://docs.docker.com/engine/install/ubuntu/)
y deben estar hechos los pasos post instalacion:
(https://docs.docker.com/engine/install/linux-postinstall/)
 	$ sudo usermod -aG docker $USER  (posiblemente haya que cerrar todas las terminales luego de correr este comando y volver a abrirlas. Para chequear qye este todo bien correr:
	$ docker run hello-world        )                          

0) git clone 

### El proceso consta de 5 scripts:
1) Prepare_reference.py (descarga y prepara la referencia de H37Rv del NCBI)

	**$ python3 prepare_reference.py**

2) Process_sample.py (Prepara las muestras tomando fastq1 de R1 y R2. Para más de 1 puede meterse en bash)

	**$ python3 process_sample.py -fastq1 muestra1_R1.fastq -fastq2 muestra1_R2.fastq**

3) Process_group.py (toma la carpeta de VCFs crudos y los une en un archivo “variant.ann.vcf”)

	**$ python3 process_group.py -i vcfs/**

4) Resistance_detection.py (toma el “variant.ann.vcf” y lo compara contra las referencias y bases originales buscando resistencias. Tiene como salida un .json ‘out_resistance’. actualmente solo trabaja con la primer muestra, en proceso para multiples)

	**$ python3 resistance_detection.py**

5) Analisis_resistencia.py (en proceso…

------------------------------------------------------------------------------------------------------------------------------------------

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
	
