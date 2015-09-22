DATA="inst/extdata"

# download data
wget -O $DATA/GSE55291_RAW.tar 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55291&format=file'
mkdir $DATA/kim
tar xvC $DATA/kim -f $DATA/GSE55291_RAW.tar
gunzip $DATA/kim/*

# labels
# GSM1333674	iPS 1
# GSM1333675	iPS 2
# GSM1333676	iPS 3
# GSM1333677	iPS 4
# GSM1333678	iPS 5
# GSM1333679	iPS 6
# GSM1333680	iPS 7
# GSM1333681	iPS 8
# GSM1333682	iPS 9
# GSM1333683	iPS 10
# GSM1333684	fES 1
# GSM1333685	fES 2
# GSM1333686	fES 3
# GSM1333687	fES 4
# GSM1333688	fES 5
# GSM1333689	fES 6
# GSM1333690	fES 7
# GSM1333691	fES 8
# GSM1333692	TTF 1
# GSM1333693	TTF 2
# GSM1333694	TTF 3
# GSM1333695	TTF 4
# GSM1333696	Wk2 1
# GSM1333697	Wk2 2
# GSM1333698	Wk2 3
# GSM1333699	Wk2 4
# GSM1333700	mES 1
# GSM1333701	mES 2
# GSM1333702	mES 3
# GSM1333703	mES 4
# GSM1333704	mES 5
# GSM1333705	mES 6
# GSM1333706	iPS 2i 1
# GSM1333707	iPS 2i 2
# GSM1333708	iPS 2i 3
# GSM1333709	iPS 2i 4
# GSM1333710	iPS 2i 5
# GSM1333711	iPS 2i 6
# GSM1333712	iPS 2i 7
# GSM1333713	iPS 2i 8
# GSM1333714	ES
# GSM1333715	TTF
# GSM1333716	iPS
# GSM1333717	iPS KD1
# GSM1333718	iPS KD2
# GSM1333719	iPS sRNA
# GSM1570590	Wk2 5
# GSM1570591	Wk2 6
# GSM1570592	Wk2 7
# GSM1570593	Wk2 8
# GSM1570594	Wk2 9
# GSM1570595	Wk2 10
# GSM1570596	Wk2 11
# GSM1570597	Wk2 12
# GSM1570598	Wk2 13
# GSM1570599	Wk2 14
# GSM1570600	iPS 11
# GSM1570601	iPS 12
# GSM1570602	iPS 13
# GSM1570603	iPS 14
# GSM1570604	iPS 15
# GSM1570605	iPS 16
# GSM1570606	iPS 17
# GSM1570607	iPS 18
# GSM1570608	iPS 19
# GSM1570609	iPS 20
# GSM1570610	TTF 5
# GSM1570611	TTF 6
# GSM1570612	TTF 7
# GSM1570613	TTF 8
# GSM1570614	TTF 9
# GSM1570615	TTF 10
# GSM1570616	TTF 11
# GSM1570617	TTF 12
# GSM1570618	TTF 13
# GSM1570619	TTF 14
# GSM1570620	ES 15
# GSM1570621	ES 16
# GSM1570622	ES 17
# GSM1570623	ES 18
# GSM1570624	ES 19
# GSM1570625	ES 20
# GSM1570626	ES 21
# GSM1570627	Wk2 15
# GSM1570628	iPS 21
# GSM1570629	iPS 22
# GSM1570630	iPS 23
# GSM1570631	iPS 24
# GSM1570632	iPS 25
# GSM1570633	iPS 26
# GSM1570634	iPS 27
# GSM1570635	iPS 28
# GSM1570636	iPS 29
# GSM1570637	iPS 30
# GSM1570638	iPS 31

# rename files
for i in $(seq 74 83)
do
    find $DATA/kim -iname "GSM13336$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/iPS_$i.txt
done

for i in $(seq 84 91)
do
    find $DATA/kim -iname "GSM13336$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/fES_$i.txt
done

for i in $(seq 92 95)
do
    find $DATA/kim -iname "GSM13336$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/TTF_$i.txt
done

for i in $(seq 96 99)
do
    find $DATA/kim -iname "GSM13336$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/Wk2_$i.txt
done

for i in $(seq 0 5)
do
    find $DATA/kim -iname "GSM133370$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/mES_$i.txt
done

for i in $(seq 0 5)
do
    find $DATA/kim -iname "GSM133370$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/mES_$i.txt
done

for i in $(seq 706 713)
do
    find $DATA/kim -iname "GSM1333$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/iPS2i_$i.txt
done

# find $DATA/kim -iname "GSM1333714*" -print0 | xargs -0 -I {} mv {} $DATA/kim/ES_$i.txt
# find $DATA/kim -iname "GSM1333715*" -print0 | xargs -0 -I {} mv {} $DATA/kim/TTF_$i.txt
# find $DATA/kim -iname "GSM1333716*" -print0 | xargs -0 -I {} mv {} $DATA/kim/iPS_$i.txt
# find $DATA/kim -iname "GSM1333717*" -print0 | xargs -0 -I {} mv {} $DATA/kim/iPSKD1_$i.txt
# find $DATA/kim -iname "GSM1333718*" -print0 | xargs -0 -I {} mv {} $DATA/kim/iPSKD2_$i.txt
# find $DATA/kim -iname "GSM1333719*" -print0 | xargs -0 -I {} mv {} $DATA/kim/iPSsRNA_$i.txt

for i in $(seq 590 599)
do
    find $DATA/kim -iname "GSM1570$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/Wk2_$i.txt
done

for i in $(seq 600 609)
do
    find $DATA/kim -iname "GSM1570$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/iPS_$i.txt
done

for i in $(seq 610 619)
do
    find $DATA/kim -iname "GSM1570$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/TTF_$i.txt
done

for i in $(seq 620 626)
do
    find $DATA/kim -iname "GSM1570$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/ES_$i.txt
done

find $DATA/kim -iname "GSM1570627*" -print0 | xargs -0 -I {} mv {} $DATA/kim/Wk2_$i.txt

for i in $(seq 628 638)
do
    find $DATA/kim -iname "GSM1570$i*" -print0 | xargs -0 -I {} mv {} $DATA/kim/iPS_$i.txt
done

# combine files

# 1 - TTF cells (Tail-tip fibroblast) - 14
# 2 - Wk2 cells (Transitional cells) - 15
# 3 - iPS cells - 31
# 4 - ES cells - 21

paste $DATA/kim/TTF* | awk '{for (i = 10; i <= NF; i += 13) printf ("%s%c", $i, i + 13 <= NF ? "\t" : "\n");}' > $DATA/kim/1-TTF.txt
sed -i '' '1s/FPKM/1/g' $DATA/kim/1-TTF.txt

paste $DATA/kim/Wk2* | awk '{for (i = 10; i <= NF; i += 13) printf ("%s%c", $i, i + 13 <= NF ? "\t" : "\n");}' > $DATA/kim/2-Wk2.txt
sed -i '' '1s/FPKM/2/g' $DATA/kim/2-Wk2.txt

paste $DATA/kim/iPS_* | awk '{for (i = 10; i <= NF; i += 13) printf ("%s%c", $i, i + 13 <= NF ? "\t" : "\n");}' > $DATA/kim/3-iPS.txt
sed -i '' '1s/FPKM/3/g' $DATA/kim/3-iPS.txt

paste $DATA/kim/*ES* | awk '{for (i = 10; i <= NF; i += 13) printf ("%s%c", $i, i + 13 <= NF ? "\t" : "\n");}' > $DATA/kim/4-ES.txt
sed -i '' '1s/FPKM/4/g' $DATA/kim/4-ES.txt

awk -F"\t" '{if ($5) print $5}' $DATA/kim/iPS_638.txt > $DATA/kim/gene_names.txt
paste $DATA/kim/1-TTF.txt $DATA/kim/2-Wk2.txt $DATA/kim/3-iPS.txt $DATA/kim/4-ES.txt > $DATA/kim/all.txt
paste $DATA/kim/gene_names.txt $DATA/kim/all.txt > $DATA/kim/kim-all-fpkm.txt
