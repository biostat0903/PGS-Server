pop=EUR
external=/home/yasheng/comprsWeb/scripts/external.sh
extsumm=/home/yasheng/comprsWeb/example_data/ext/extsumm.ldsc.gz
esteff=/home/yasheng/comprsWeb/example_data/output/lassosum_esteff.txt,1,2,3
outpath=/home/yasheng/comprsWeb/example_data/ext/

sh ${external} -s ${esteff} -e ${extsumm} -p ${pop} -o ${outpath}