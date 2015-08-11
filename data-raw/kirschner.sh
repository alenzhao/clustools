DATA="inst/extdata"

# download data
wget -O $DATA/GSE65525.tar 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE65525&format=file'
mkdir $DATA/kirschner
tar xvC $DATA/kirschner -f $DATA/GSE65525.tar
bzip2 -d $DATA/kirschner/*

head -1 $DATA/kirschner/GSM1599494_ES_d0_main.csv | sed 's/[^,]//g' | wc -c > $DATA/kirschner/d0.txt
head -1 $DATA/kirschner/GSM1599497_ES_d2_LIFminus.csv | sed 's/[^,]//g' | wc -c > $DATA/kirschner/d2.txt
head -1 $DATA/kirschner/GSM1599498_ES_d4_LIFminus.csv | sed 's/[^,]//g' | wc -c > $DATA/kirschner/d4.txt
head -1 $DATA/kirschner/GSM1599499_ES_d7_LIFminus.csv | sed 's/[^,]//g' | wc -c > $DATA/kirschner/d7.txt

sed -i '' 's/[^,]*,//' $DATA/kirschner/GSM1599497_ES_d2_LIFminus.csv
sed -i '' 's/[^,]*,//' $DATA/kirschner/GSM1599498_ES_d4_LIFminus.csv
sed -i '' 's/[^,]*,//' $DATA/kirschner/GSM1599499_ES_d7_LIFminus.csv

paste -d"," $DATA/kirschner/GSM1599494_ES_d0_main.csv $DATA/kirschner/GSM1599497_ES_d2_LIFminus.csv $DATA/kirschner/GSM1599498_ES_d4_LIFminus.csv $DATA/kirschner/GSM1599499_ES_d7_LIFminus.csv > $DATA/kirschner/reads.txt
