# tb_resistance

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
