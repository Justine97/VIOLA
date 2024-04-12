Run EXOMISER
================

We run **Exomiser** on IRCAN server. <br> To run Exomiser, we need
*positive_wes-analysis-batch-commands.txt* file which contains paths of
all config files of all the patients of the cohort.

``` bash
cd /scratch/bin/exomiser/

java -Xms2g -Xmx4g -jar exomiser-cli-13.1.0.jar --batch /home/jlabory/data/positive_wes_vcf/yml_file/maf_0.01/positive_wes-analysis-batch-commands.txt --spring.config.location=/scratch/bin/exomiser/application.properties
```
