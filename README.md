# main title here

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

## Lab 5: Gene Family Phylogeny using IQ-TREE
### Purpose of the Lab 5

## Lab 6. Reconciling a Gene and Species Tree
### Purpose of the Lab 6

## Lab 8: Protein Domain Prediction
### Purpose of the Lab 8
