#!/usr/bin/env bash

mkdir -p merged

fmsuff1=_L001_R1_001.fastq.gz
fmsuff2=$(sed 's/L001/L002/' <<< ${fmsuff1})
fmsuff3=$(sed 's/L001/L003/' <<< ${fmsuff1})
fmsuff4=$(sed 's/L001/L004/' <<< ${fmsuff1})

suff1=_1.fastq.gz
suff2=$(sed 's/1/2/' <<< ${suff1})

rvmsuff1=$(sed 's/1/2/' <<< ${fmsuff1})
rvmsuff2=$(sed 's/L001/L002/' <<< ${rvmsuff1})
rvmsuff3=$(sed 's/L001/L003/' <<< ${rvmsuff1})
rvmsuff4=$(sed 's/L001/L004/' <<< ${rvmsuff1})

count=`ls -1 *${fmsuff1} | wc -l`

num=0

if [ $count != 0 ];
then
for rd1 in *${fmsuff1}
do
	num=$(( $num + 1 ))

	base=$(basename ${rd1} ${fmsuff1})
	cat ${base}${fmsuff1} ${base}${fmsuff2} ${base}${fmsuff3} ${base}${fmsuff4} > merged/${base}${suff1}
	cat ${base}${rvmsuff1} ${base}${rvmsuff2} ${base}${rvmsuff3} ${base}${rvmsuff4} > merged/${base}${suff2}
done

echo "Merged ${num} files, see the 'merged' directory"

else 
	echo "There is no file with ${fmsuff1} suffix in the directory"
fi
