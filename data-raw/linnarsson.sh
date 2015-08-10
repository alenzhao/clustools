DATA="inst/extdata"

# download data
wget -O $DATA/linnarsson.txt 'http://linnarssonlab.org/blobs/cortex/expression_mRNA_17-Aug-2014.txt'

# remove some headers
sed -i '' '1d' $DATA/linnarsson.txt
sed -i '' '2,10d' $DATA/linnarsson.txt

# assign labels
sed -i '' '1s/ #//' $DATA/linnarsson.txt
