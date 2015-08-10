DATA="inst/extdata"

# download data
wget -P $DATA 'http://www.nature.com/nature/journal/v509/n7500/extref/nature13173-s4.txt'

# transpose matrix
awk  '{
for (f = 1; f <= NF; f++)
a[NR, f] = $f
}
NF > nf { nf = NF }
END {
for (f = 1; f <= nf; f++)
for (r = 1; r <= NR; r++)
printf a[r, f] (r==NR ? RS : FS)
}' $DATA/nature13173-s4.txt > $DATA/quake.txt
rm $DATA/nature13173-s4.txt

# remove some headers
sed -i '' '1,3d' $DATA/quake.txt

# assign labels
sed -i '' '1s/Clara/1/g' $DATA/quake.txt
sed -i '' '1s/AT1/2/g' $DATA/quake.txt
sed -i '' '1s/AT2/3/g' $DATA/quake.txt
sed -i '' '1s/BP/4/g' $DATA/quake.txt
sed -i '' '1s/ciliated/5/g' $DATA/quake.txt
sed -i '' '1s/"putative_cell_type" //' $DATA/quake.txt
