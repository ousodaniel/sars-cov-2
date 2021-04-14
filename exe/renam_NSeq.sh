#!/usr/bin/env bash

mkdir -p merged

fmsuff1=_L001_R1_001.fastq.gz

count=`ls -1 *${fmsuff1} | wc -l`

num=0

suff1=_1.fastq.gz
suff2=$(sed 's/1/2/' <<< ${suff1})

if [ $count != 0 ];
then
for rd1 in *${fmsuff1}
do
  	rrd1=$(sed 's/R1/R2/' <<< $rd1)
        num=$(( $num + 1 ))
        rd2=$(sed 's/L001/L002/' <<< $rd1)
        rd3=$(sed 's/L001/L003/' <<< $rd1)
        rd4=$(sed 's/L001/L004/' <<< $rd1)

        rrd2=$(sed 's/L001/L002/' <<< $rrd1)
        rrd3=$(sed 's/L001/L003/' <<< $rrd1)
        rrd4=$(sed 's/L001/L004/' <<< $rrd1)

        base=$(basename ${rd1} ${fmsuff1})
        cat $rd1 $rd2 $rd3 $rd4 > merged/${base}_1.fastq.gz
        cat $rrd1 $rrd2 $rrd3 $rrd4 > merged/${base}_2.fastq.gz
done

echo "${num} files to be merged were founs, see output($num) in the 'merged' directory"

else
    	echo "There is no file with ${fmsuff1} suffix in the directory"
fi
