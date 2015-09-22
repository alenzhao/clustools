DATA="inst/extdata"

# download data
wget -O $DATA/GSE51372_readCounts.txt.gz 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE51372&format=file&file=GSE51372%5FreadCounts%2Etxt%2Egz'
gunzip $DATA/GSE51372_readCounts.txt.gz
# sed -i '' 's/://' $DATA/GSE51372_readCounts.txt
# sed -i '' 's/"//' $DATA/GSE51372_readCounts.txt
# sed -i '' 's/-//' $DATA/GSE51372_readCounts.txt
# sed -i '' 's/\.//' $DATA/GSE51372_readCounts.txt
# sed -i '' 's/(//' $DATA/GSE51372_readCounts.txt
# sed -i '' 's/)//' $DATA/GSE51372_readCounts.txt
# sed -i '' 's/\///' $DATA/GSE51372_readCounts.txt
# sed -i '' 's/_//' $DATA/GSE51372_readCounts.txt
# sed -i '' 's/,//' $DATA/GSE51372_readCounts.txt
