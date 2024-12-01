# Evlutionary and Functional Analysis of EF-hand Motif Genes

## *The evolutionary, structural, and functional diversity of EF-hand motifs, NP_001364925.1*

## Lab 3: Finding homologs with BLAST KEY

### Purpose of the Lab 3
*Use BLAST to identify homologous sequences (orthologs and paralogs) of a query protein from multiple species using NP_001364925.1*

*Call Lab 3*
```
git clone https://github.com/Bio312/lab03-$MYGIT

cd lab03-$MYGIT
```

*BLAST a globin protein against the database to identify homologs using NP_001364925.1*

*Create a new working directory*
```
mkdir ~/lab03-$MYGIT/globins

cd ~/lab03-$MYGIT/globins

pwd
```
*Download, NP_001364925.1, my protein data from NCBI, and perform BLAST*
```
ncbi-acc-download -F fasta -m protein "NP_001364925.1"

blastp -db ../allprotein.fas -query NP_001364925.1.fa -outfmt 0 -max_hsps 1 -out globins.blastp.typical.out
```
*Perform a BLAST search, and request tabular output*
```
blastp -db ../allprotein.fas -query NP_000549.1.fa  -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out globins.blastp.detail.out

grep -c H.sapiens globins.blastp.detail.out
```

*Filter our output file to satisfy this requirement*
```
awk '{if ($6< 1e-30)print $1 }' globins.blastp.detail.out > globins.blastp.detail.filtered.out
```

*Count the total number of hits in the BLAST results after the filter*

```
wc -l globins.blastp.detail.filtered.out
```
*To know how may paralogs are found in each specie*
```
grep -o -E "^[A-Z]\.[a-z]+" globins.blastp.detail.filtered.out  | sort | uniq -c
```

### Save and Push (End of the Lab 3)
```
cd ~/lab03-$MYGIT

history > ~/lab03-$MYGIT/lab3.commandhistory.txt

find . -size +5M | sed 's|^\./||g' | cat >> .gitignore; awk '!NF || !seen[$0]++' .gitignore

git add .

git commit -a -m "Adding all new data files I generated in AWS to the repository."

git pull --rebase=false --no-edit

git push 
```


## Lab 4: Gene family sequence alignment
### Purpose of the Lab 4
*Aligning gene family sequences, aalyzing quality, and exploring evoloutionary relationships*

*Call Lab 4*
```
git clone https://github.com/Bio312/lab04-$MYGIT

cd lab04-$MYGIT
```

*Obtain  the sequences for the globin gene mem bers that were BLAST hits to my gene with an e-value less than 1e-30*
```
mkdir ~/lab04-$MYGIT/globins

cd ~/lab04-$MYGIT/globins

pwd
```

*Here is the squkit command to obtain the sequences that are in the BLAST output file*
```
seqkit grep --pattern-file ~/lab03-$MYGIT/globins/globins.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/globins/globins.homologs.fas
```

*Perform a global multiple sequence alignment in muscle and save it as a PDF*
```
muscle -align ~/lab04-$MYGIT/globins/globins.homologs.fas -output ~/lab04-$MYGIT/globins/globins.homologs.al.fas

alv -kli  ~/lab04-$MYGIT/globins/globins.homologs.al.fas | less -RS
```

*(for majority option)*
```
alv -kli --majority ~/lab04-$MYGIT/globins/globins.homologs.al.fas | less -RS
```

*using R script for printing the alignment*
```
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R  ~/lab04-$MYGIT/globins/globins.homologs.al.fas
```

### Get more information about the alignment
*Calculate the width(length) of the alignment*
```
alignbuddy  -al  ~/lab04-$MYGIT/globins/globins.homologs.al.fas

```

*To calculate the length of the alignment after removing any column with gaps*
```
alignbuddy -trm all  ~/lab04-$MYGIT/globins/globins.homologs.al.fas | alignbuddy  -al
```

*To calculate the length of the alignment after removing invariant (completely conserved) position*
```
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/globins/globins.homologs.al.fas | alignbuddy  -al
```

### Calculate the average percent identity
*Use t_coffee to calculate the average percent identity*
```
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/globins/globins.homologs.al.fas -output sim
```

*repeat calculating the average percent identity alignbuddy*
```
alignbuddy -pi ~/lab04-$MYGIT/globins/globins.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```

### Save and Push (end of the Lab 4)
```
cd ~/lab04-$MYGIT

find . -size +5M | sed 's|^\./||g' | cat >> .gitignore; awk '!NF || !seen[$0]++' .gitignore

git add .

git commit -a -m "Adding all new data files I generated in AWS to the repository."

git pull --rebase=false --no-edit

git push 
```

## Lab 5: Gene Family Phylogeny using IQ-TREE
### Purpose of the Lab 5

## Lab 6. Reconciling a Gene and Species Tree
### Purpose of the Lab 6

## Lab 8: Protein Domain Prediction
### Purpose of the Lab 8
