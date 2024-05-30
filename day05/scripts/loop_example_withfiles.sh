#!/bin/bash

#########################
## EXAMPLE WITH FOLDER ##

FOLDER=/path/to/your/data_files/fastq
for index in ${ls $FOLDER}
do
echo ${index}

done



