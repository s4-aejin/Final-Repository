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
*Making phylogenetic trees for gene familes using IQ-TREE, rooting methods, a d visualization tools*

*Call Lab 5*
```
git clone https://github.com/Bio312/lab05-$MYGIT

cd lab05-$MYGIT
```

*Making a directory in lab 5 for the globin tree*
```
mkdir ~/lab05-$MYGIT/globins

cd ~/lab05-$MYGIT/globins  
```

*Make sure that sequences with duolicate labels are removed and create a clean alignment file for use in Lab 5* 
```
sed 's/ /_/g'  ~/lab04-$MYGIT/globins/globins.homologs.al.fas | seqkit grep -v -r -p "dupelabel" >  ~/lab05-$MYGIT/globins/globins.homologsf.al.fas
```

*Use IQ-Tree to calculate the optim amino acid substitution model an d amin o acid frequen cies and estimate branch lengths as they go*
```
iqtree -s ~/lab05-$MYGIT/globins/globins.homologsf.al.fas -bb 1000 -nt 2 
```

*Substitution model - read the tree in a Newick format*
```
nw_display ~/lab05-$MYGIT/globins/globins.homologsf.al.fas.treefile
```

*Look at the unrooted tree - use R script to look at the unrooted tree with a graphical display*
```
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R  ~/lab05-$MYGIT/globins/globins.homologsf.al.fas.treefile ~/lab05-$MYGIT/globins/globins.homologsf.al.fas.treefile.pdf 0.4 15
```

*Root in the optimal phylogeny*

*1. Midpoint rooting*
```
gotree reroot midpoint -i ~/lab05-$MYGIT/globins/globins.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/globins/globins.homologsf.al.mid.treefile

nw_order -c n ~/lab05-$MYGIT/globins/globins.homologsf.al.mid.treefile  | nw_display -
```

*make graphic image*
```
nw_order -c n ~/lab05-$MYGIT/globins/globins.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s  >  ~/lab05-$MYGIT/globins/globins.homologsf.al.mid.treefile.svg -
```

*make PDF*
```
convert  ~/lab05-$MYGIT/globins/globins.homologsf.al.mid.treefile.svg  ~/lab05-$MYGIT/globins/globins.homologsf.al.mid.treefile.pdf
```

*2. Branch lengths*
*switching the view to a cladogram*
```
nw_order -c n ~/lab05-$MYGIT/globins/globins.homologsf.al.mid.treefile | nw_topology - | nw_display -s  -w 1000 > ~/lab05-$MYGIT/globins/globins.homologsf.al.midCl.treefile.svg -

convert ~/lab05-$MYGIT/globins/globins.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/globins/globins.homologsf.al.midCl.treefile.pdf
```

*3. Outgroup rooting*
*specifying an outgroup, a group that we know based on outside evidence is more distantly relat4ed than any of lineages in the ingroup*
```
nw_reroot ~/lab05-$MYGIT/globins/globins.homologsf.al.fas.treefile H.sapiens_HBG1_hemoglobin_subunit_gamma1 H.sapiens_HBG2_hemoglobin_subunit_gamma2 H.sapiens_HBB_hemoglobin_subunit_beta H.sapiens_HBD_hemoglobin_subunit_delta >~/lab05-$MYGIT/globins/globins.homologsf.outgroupbeta.treefile

nw_order -c n ~/lab05-$MYGIT/globins/globins.homologsf.outgroupbeta.treefile | nw_topology - | nw_display -s  -w 1000 > ~/lab05-$MYGIT/globins/globins.homologsf.outgroupbeta.treefile.svg -

convert ~/lab05-$MYGIT/globins/globins.homologsf.outgroupbeta.treefile.svg ~/lab05-$MYGIT/globins/globins.homologsf.outgroupbeta.treefile.pdf
```

### Save and Push (end of the Lab 5)
```
history > lab5.commandhistory.txt

cd ~/lab05-$MYGIT

find . -size +5M | sed 's|^\./||g' | cat >> .gitignore; awk '!NF || !seen[$0]++' .gitignore

git add .

git commit -a -m "Adding all new data files I generated in AWS to the repository."

git pull --no-edit

git push 
```

## Lab 6. Reconciling a Gene and Species Tree
### Purpose of the Lab 6
*reconciling gene and species trees to identify duplication and loss events and understand gene family evolution*

*Call Lab 6 and a directory called globins*
```
git clone https://github.com/Bio312/lab06-$MYGIT

cd lab06-$MYGIT

cd ~/lab06-$MYGIT/globins

ls ~/lab06-$MYGIT/globins/globins.homologsf.outgroupbeta.treefile 

```

*Call tree from globins directory and drop the beta globin (outgroup) species because we are onlyh intterested in alpha globins (ingroup)*
```
gotree prune -i ~/lab06-$MYGIT/globins/globins.homologsf.outgroupbeta.treefile -o ~/lab06-$MYGIT/globins/globins.homologsf.pruned.treefile H.sapiens_HBG1_hemoglobin_subunit_gamma1 H.sapiens_HBG2_hemoglobin_subunit_gamma2 H.sapiens_HBB_hemoglobin_subunit_beta H.sapiens_HBD_hemoglobin_subunit_delta
```

*Create a new folder*
```
mkdir ~/lab06-$MYGIT/mygenefamily

cp ~/lab05-$MYGIT/mygenefamily/mygenefamily.homologs.al.mid.treefile ~/lab06-$MYGIT/mygenefamily/mygenefamily.homologs.al.mid.treefile
```

*Reconcile the gene and species tree using Notnug*
```
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/globins/globins.homologsf.pruned.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/globins/

nw_display ~/lab05-$MYGIT/species.tre

grep NOTUNG-SPECIES-TREE ~/lab06-$MYGIT/globins/globins.homologsf.pruned.treefile.rec.ntg | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -
```

*Generate a RecPhyloXML object and view the gene within species tree via thridkind
```
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/globins/globins.homologsf.pruned.treefile.rec.ntg --include.species

thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/globins/globins.homologsf.pruned.treefile.rec.ntg.xml -o  ~/lab06-$MYGIT/globins/globins.homologsf.pruned.treefile.rec.svg

convert  -density 150 ~/lab06-$MYGIT/globins/globins.homologsf.pruned.treefile.rec.svg ~/lab06-$MYGIT/globins/globins.homologsf.pruned.treefile.rec.pdf
```

### Save and Push (end of the Lab 6)
```
history > lab6.commandhistory.txt

cd ~/lab06-$MYGIT

find . -size +8M | sed 's|^\./||g' | cat >> .gitignore; awk '!NF || !seen[$0]++' .gitignore

git add.

git status

git commit -a -m "Adding all new data files I generated in AWS to the repository."

git pull --no-edit

git push 
```

## Lab 8: Protein Domain Prediction
### Purpose of the Lab 8
*Identyfying protein domains using RPS-BLAST, mapping them onto phylogenetic trees, and analyzing domain evolution in gene familes*

*Call Lab 8*
```
git clone https://github.com/Bio312/lab08-$MYGIT

cd lab08-$MYGIT
```

*Make a directory for the globin sequences and make a copy of raw unaligned sequen ces, removing the asterisk (stop codon)*
```
mkdir ~/lab08-$MYGIT/globins && cd ~/lab08-$MYGIT/globins

sed 's/*//' ~/lab04-$MYGIT/globins/globins.homologs.fas > ~/lab08-$MYGIT/globins/globins.homologs.fas
```

*We have Pfam database and will run RPS-BLAST*
```
rpsblast -query ~/lab08-$MYGIT/globins/globins.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/globins/globins.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```

*Using gene tree from lab 5, we can plot the predicted Pfam domains on the phylogeny.*
*The script that Dr.Rest wrote, plotTreeAndDomains.r, was used.*
```
cp ~/lab05-$MYGIT/globins/globins.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/globins

Rscript  --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/globins/globins.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/globins/globins.rps-blast.out ~/lab08-$MYGIT/globins/globins.tree.rps.pdf
```

*The domain-on-tree graphic that we just produced using Dr.Rest's R script and loot at the predited domains in more detail.
Can find which Pfam domain annotation is most commonly found, the shortest annotated protein domain or the longest one, and the best e-vlaue*

```
mlr --inidx --ifs "\t" --opprint  cat ~/lab08-$MYGIT/globins/globins.rps-blast.out | tail -n +2 | less -S

cut -f 1 ~/lab08-$MYGIT/globins/globins.rps-blast.out | sort | uniq -c

cut -f 6 ~/lab08-$MYGIT/globins/globins.rps-blast.out | sort | uniq -c

awk '{a=$4-$3;print $1,'\t',a;}' ~/lab08-$MYGIT/globins/globins.rps-blast.out |  sort  -k2nr

cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/globins/globins.rps-blast.out 
```


### Save and Push (end of the Lab 8)

```
history > lab8.commandhistory.txt

cd ~/lab08-$MYGIT

find . -size +5M | sed 's|^\./||g' | cat >> .gitignore; awk '!NF || !seen[$0]++' .gitignore

git add .

git status

git commit -a -m "Adding all new data files I generated in AWS to the repository."

git pull --no-edit

git push
```