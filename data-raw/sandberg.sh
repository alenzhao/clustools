DATA="inst/extdata"

# download data
wget -O $DATA/GSE45719_RAW.tar 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE45719&format=file'
mkdir $DATA/sandberg
tar xvC $DATA/sandberg -f $DATA/GSE45719_RAW.tar
gunzip $DATA/sandberg/*

# cells:
#  1   1-4 - zygote
#  2   5-12 - early2cell
#  3   13-24 - mid2cell
#  4   25-34 - late2cell
#  5   35-48 - 4cell
#  6   49-85 - 8cell
#  7   86-135 - 16cell
#  8   136-178 - earlyblast
#  9   179-238 - midblast
#  10   239-268 - lateblast

# Raw READS
paste $DATA/sandberg/GSM111*_zy* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/1_reads-zy.txt
sed -i '' '1s/reads/1/g' $DATA/sandberg/1_reads-zy.txt
paste $DATA/sandberg/GSM111*_early2cell_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/2_reads-early2cell.txt
sed -i '' '1s/reads/2/g' $DATA/sandberg/2_reads-early2cell.txt
paste $DATA/sandberg/GSM111*_mid2cell_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/3_reads-mid2cell.txt
sed -i '' '1s/reads/3/g' $DATA/sandberg/3_reads-mid2cell.txt
paste $DATA/sandberg/GSM111*_late2cell_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/4_reads-late2cell.txt
sed -i '' '1s/reads/4/g' $DATA/sandberg/4_reads-late2cell.txt
paste $DATA/sandberg/GSM111*_4cell_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/5_reads-4cell.txt
sed -i '' '1s/reads/5/g' $DATA/sandberg/5_reads-4cell.txt
paste $DATA/sandberg/*_8cell_*-* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/6_reads-8cell.txt
sed -i '' '1s/reads/6/g' $DATA/sandberg/6_reads-8cell.txt
paste $DATA/sandberg/GSM111*_16cell_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/7_reads-16cell.txt
sed -i '' '1s/reads/7/g' $DATA/sandberg/7_reads-16cell.txt
paste $DATA/sandberg/GSM111*_earlyblast_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/8_reads-earlyblast.txt
sed -i '' '1s/reads/8/g' $DATA/sandberg/8_reads-earlyblast.txt
paste $DATA/sandberg/GSM111*_midblast_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/91_reads-midblast.txt
sed -i '' '1s/reads/9/g' $DATA/sandberg/91_reads-midblast.txt
paste $DATA/sandberg/GSM111*_lateblast_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/92_reads-lateblast.txt
sed -i '' '1s/reads/10/g' $DATA/sandberg/92_reads-lateblast.txt

awk -F"\t" '{if ($1) print $1}' $DATA/sandberg/GSM1112767_zy2_expression.txt > $DATA/sandberg/gene_names.txt

paste $DATA/sandberg/*_reads-* > $DATA/sandberg/all-reads.txt
paste $DATA/sandberg/gene_names.txt $DATA/sandberg/all-reads.txt > $DATA/sandberg/sandberg-all-reads.txt
sed -i '' '1s/^#//' $DATA/sandberg/sandberg-all-reads.txt

# RPKMs
paste $DATA/sandberg/GSM111*_zy* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/1_reads-zy.txt
sed -i '' '1s/RPKM/1/g' $DATA/sandberg/1_reads-zy.txt
paste $DATA/sandberg/GSM111*_early2cell_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/2_reads-early2cell.txt
sed -i '' '1s/RPKM/2/g' $DATA/sandberg/2_reads-early2cell.txt
paste $DATA/sandberg/GSM111*_mid2cell_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/3_reads-mid2cell.txt
sed -i '' '1s/RPKM/3/g' $DATA/sandberg/3_reads-mid2cell.txt
paste $DATA/sandberg/GSM111*_late2cell_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/4_reads-late2cell.txt
sed -i '' '1s/RPKM/4/g' $DATA/sandberg/4_reads-late2cell.txt
paste $DATA/sandberg/GSM111*_4cell_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/5_reads-4cell.txt
sed -i '' '1s/RPKM/5/g' $DATA/sandberg/5_reads-4cell.txt
paste $DATA/sandberg/*_8cell_*-* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/6_reads-8cell.txt
sed -i '' '1s/RPKM/6/g' $DATA/sandberg/6_reads-8cell.txt
paste $DATA/sandberg/GSM111*_16cell_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/7_reads-16cell.txt
sed -i '' '1s/RPKM/7/g' $DATA/sandberg/7_reads-16cell.txt
paste $DATA/sandberg/GSM111*_earlyblast_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/8_reads-earlyblast.txt
sed -i '' '1s/RPKM/8/g' $DATA/sandberg/8_reads-earlyblast.txt
paste $DATA/sandberg/GSM111*_midblast_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/91_reads-midblast.txt
sed -i '' '1s/RPKM/9/g' $DATA/sandberg/91_reads-midblast.txt
paste $DATA/sandberg/GSM111*_lateblast_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > $DATA/sandberg/92_reads-lateblast.txt
sed -i '' '1s/RPKM/10/g' $DATA/sandberg/92_reads-lateblast.txt

paste $DATA/sandberg/*_reads-* > $DATA/sandberg/all-rpkms.txt
paste $DATA/sandberg/gene_names.txt $DATA/sandberg/all-rpkms.txt > $DATA/sandberg/sandberg-all-rpkms.txt
sed -i '' '1s/^#//' $DATA/sandberg/sandberg-all-rpkms.txt
