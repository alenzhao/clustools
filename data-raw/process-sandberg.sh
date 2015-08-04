# Raw READS

# cluster 1
paste inst/extdata/GSE45719_RAW/GSM111*_zy* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/1_reads-zy.txt
sed -i '' '1s/reads/1/g' inst/extdata/GSE45719_RAW/1_reads-zy.txt
# cluster 2
paste inst/extdata/GSE45719_RAW/GSM111*_early2cell_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/2_reads-early2cell.txt
sed -i '' '1s/reads/2/g' inst/extdata/GSE45719_RAW/2_reads-early2cell.txt
# cluster 3
paste inst/extdata/GSE45719_RAW/GSM111*_mid2cell_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/3_reads-mid2cell.txt
sed -i '' '1s/reads/3/g' inst/extdata/GSE45719_RAW/3_reads-mid2cell.txt
# cluster 4
paste inst/extdata/GSE45719_RAW/GSM111*_late2cell_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/4_reads-late2cell.txt
sed -i '' '1s/reads/4/g' inst/extdata/GSE45719_RAW/4_reads-late2cell.txt
# cluster 5
paste inst/extdata/GSE45719_RAW/GSM111*_4cell_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/5_reads-4cell.txt
sed -i '' '1s/reads/5/g' inst/extdata/GSE45719_RAW/5_reads-4cell.txt
# cluster 6
paste inst/extdata/GSE45719_RAW/GSM111*_8cell_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/6_reads-8cell.txt
sed -i '' '1s/reads/6/g' inst/extdata/GSE45719_RAW/6_reads-8cell.txt
# cluster 7
paste inst/extdata/GSE45719_RAW/GSM111*_16cell_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/7_reads-16cell.txt
sed -i '' '1s/reads/7/g' inst/extdata/GSE45719_RAW/7_reads-16cell.txt
# cluster 8
paste inst/extdata/GSE45719_RAW/GSM111*_earlyblast_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/8_reads-earlyblast.txt
sed -i '' '1s/reads/8/g' inst/extdata/GSE45719_RAW/8_reads-earlyblast.txt
# cluster 9
paste inst/extdata/GSE45719_RAW/GSM111*_midblast_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/91_reads-midblast.txt
sed -i '' '1s/reads/9/g' inst/extdata/GSE45719_RAW/91_reads-midblast.txt
# cluster 10
paste inst/extdata/GSE45719_RAW/GSM111*_lateblast_* | awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/92_reads-lateblast.txt
sed -i '' '1s/reads/10/g' inst/extdata/GSE45719_RAW/92_reads-lateblast.txt

paste inst/extdata/GSE45719_RAW/*_reads-* > inst/extdata/GSE45719_RAW/all-reads.txt
rm inst/extdata/GSE45719_RAW/*_reads-*

# RPKMs

# cluster 1
paste inst/extdata/GSE45719_RAW/GSM111*_zy* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/1_reads-zy.txt
sed -i '' '1s/RPKM/1/g' inst/extdata/GSE45719_RAW/1_reads-zy.txt
# cluster 2
paste inst/extdata/GSE45719_RAW/GSM111*_early2cell_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/2_reads-early2cell.txt
sed -i '' '1s/RPKM/2/g' inst/extdata/GSE45719_RAW/2_reads-early2cell.txt
# cluster 3
paste inst/extdata/GSE45719_RAW/GSM111*_mid2cell_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/3_reads-mid2cell.txt
sed -i '' '1s/RPKM/3/g' inst/extdata/GSE45719_RAW/3_reads-mid2cell.txt
# cluster 4
paste inst/extdata/GSE45719_RAW/GSM111*_late2cell_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/4_reads-late2cell.txt
sed -i '' '1s/RPKM/4/g' inst/extdata/GSE45719_RAW/4_reads-late2cell.txt
# cluster 5
paste inst/extdata/GSE45719_RAW/GSM111*_4cell_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/5_reads-4cell.txt
sed -i '' '1s/RPKM/5/g' inst/extdata/GSE45719_RAW/5_reads-4cell.txt
# cluster 6
paste inst/extdata/GSE45719_RAW/GSM111*_8cell_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/6_reads-8cell.txt
sed -i '' '1s/RPKM/6/g' inst/extdata/GSE45719_RAW/6_reads-8cell.txt
# cluster 7
paste inst/extdata/GSE45719_RAW/GSM111*_16cell_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/7_reads-16cell.txt
sed -i '' '1s/RPKM/7/g' inst/extdata/GSE45719_RAW/7_reads-16cell.txt
# cluster 8
paste inst/extdata/GSE45719_RAW/GSM111*_earlyblast_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/8_reads-earlyblast.txt
sed -i '' '1s/RPKM/8/g' inst/extdata/GSE45719_RAW/8_reads-earlyblast.txt
# cluster 9
paste inst/extdata/GSE45719_RAW/GSM111*_midblast_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/91_reads-midblast.txt
sed -i '' '1s/RPKM/9/g' inst/extdata/GSE45719_RAW/91_reads-midblast.txt
# cluster 10
paste inst/extdata/GSE45719_RAW/GSM111*_lateblast_* | awk '{for (i = 3; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > inst/extdata/GSE45719_RAW/92_reads-lateblast.txt
sed -i '' '1s/RPKM/10/g' inst/extdata/GSE45719_RAW/92_reads-lateblast.txt

paste inst/extdata/GSE45719_RAW/*_reads-* > inst/extdata/GSE45719_RAW/all-rpkms.txt
rm inst/extdata/GSE45719_RAW/*_reads-*
