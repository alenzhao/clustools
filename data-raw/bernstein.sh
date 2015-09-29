DATA="inst/extdata"

# download data
wget -O $DATA/bernstein.txt.gz 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE57872&format=file&file=GSE57872%5FGBM%5Fdata%5Fmatrix%2Etxt%2Egz'
gunzip $DATA/bernstein.txt.gz
