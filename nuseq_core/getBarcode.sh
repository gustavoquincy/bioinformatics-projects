#!/bin/zsh

for i in 6 7 8 9 10 11 12 13
do 
awk -v var="$i" -F '|' 'FNR==var { print $3,$5 }' bcfile.csv >> barcodeFile.txt
done 

