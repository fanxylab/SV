#!/usr/bin/bash
#*********************************************************
# @File    : bwaAlign.sh			         *				
# @Author  : WANG Jun                                    *
# @Date    : 2021/03                                     *
# @E-mail  : oucwj@outlook.com                           *               
#*********************************************************

OUT=$1
REF=$2
IN=$3
echo $OUT
echo $REF
echo $IN

cd $OUT
bwa index $REF
cd $IN
for cell in $(ls -d */ | cut -f1 -d'/') ;
do
fqPath=$IN/$cell/SS2/CutAdapt
file1=$fqPath/*trimmed_1.fq.gz
file2=$fqPath/*trimmed_2.fq.gz
echo $file1
samfile=$OUT/$cell'.sam'
tmpfile=$OUT/$cell'_tmp.bam'
bamfile=$OUT/$cell'.bam'
sortedbam=$OUT/$cell'.sorted.bam'
bwa mem ${REF} ${file1} ${file2} > $samfile
samtools view -bS $samfile > $tmpfile
samtools view -b -F 4 $tmpfile > $bamfile
rm $tmpfile
rm $samfile
samtools sort -@ 20 -o $sortedbam $bamfile
samtools index $sortedbam
done

