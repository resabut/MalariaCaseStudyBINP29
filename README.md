Malaria case study
Author: Joan Escriva Font

# conda
```bash
conda create -n malaria proteinortho=6.0.33
conda create -n malaria_busco -c bioconda busco=5.4.5
```
# git-lfs
git lfs track *.genome *.gff *.gtf


### Adding genomes
In the server  
Copy the genome data
```bash
mkdir Data Data/Genomes
scp /home2/resources/binp29/Data/malaria/plasmodiumGenomes.tgz Data/Genomes/
tar zxvf Data/Genomes/plasmodiumGenomes.tgz -C Data/Genomes/
rm Data/Genomes/plasmodiumGenomes.tgz # remove the file
scp /home2/resources/binp29/Data/malaria/Haemoproteus_tartakovskyi.raw.genome.gz Data/Genomes/
gunzip Data/Genomes/Haemoproteus_tartakovskyi.raw.genome.gz
```
    Plasmodium_berghei.genome
    Plasmodium_cynomolgi.genome
    Plasmodium_faciparum.genome
    Plasmodium_knowlesi.genome
    Plasmodium_vivax.genome
    Plasmodium_yoelii.genome
    Toxoplasma_gondii.genome
### Add scripts
```bash
mkdir Scripts
scp /resources/binp29/Data/malaria/*.py Scripts
scp /home2/resources/binp28/Data/gffParse.pl Scripts
chmod +x Scripts/gffParse.pl
```
### Adding tax data
```bash
mkdir Data/tax_data
scp /home2/resources/binp29/Data/malaria/taxonomy.dat Data/tax_data

```


#### Gene prediction
Other students ran it for us.
```bash
gmes_petap.pl
```

I copy them from the shared folder.
```bash
mkdir Data/Predictions
scp /home2/resources/binp29/Data/malaria/Tg.gff.gz Data/Predictions
gunzip Data/Predictions/Tg.gff.gz
scp /tmp/Prediction/*.gtf Data/Predictions
```


## Processing of *Haemoproteus tartakovskyi* data
### Clean genome sequence
The GC threshold should be 35 for more astringency.

    Bird: 41 %  
    Parasite: 19-42%  

Average for  the whole genome before filtering
```bash
awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} \
  END{ printf "%.2f%%\n", (gc*100)/(gc+at) }' Data/Genomes/Haem*
```
27.40% GC

Filtering using the provided script.
```bash
mkdir Results Results/01_clean
python3 Scripts/removeScaffold.py \
  Data/Genomes/Haemoproteus_tartakovskyi.raw.genome \  # input assembly
  35 \  # GC threshold
  Results/01_clean/Ht.genome \  # output fltered file 
  3000  # minimum scaffold length
```

To check the number of scaffolds in the clean file
```bash
cat Results/01_clean/Ht.genome | grep -c \>
```

After the filtering
```bash
awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} \
  END{ printf "%.2f%%\n", (gc*100)/(gc+at) }' Results/01_clean/Ht.genome
```
25.92% GC

### Make a gene prediction
Firstly, GeneMark runs the GeneMark-ES self-training algorithm and predicts the genes.
```bash
mkdir Results/02_gmes Results/02_gmes/HT_genome_1st_run
nohup gmes_petap.pl \
 --ES \  # run gene self-training algorithm
 --sequence Results/01_clean/Ht.genome \  # genome
 --cores 4 \  # cores to run
 --work_dir Results/02_gmes/HT_genome_1st_run \  # output folder
 --min_contig 3000 & # same as the previous threshold
```
Teh output file is under ./Results/02_gmes/HT_genome_1st_run/genemark.gtf   

To convert the gtf files to gff, gffParse.pl is used
```bash
cat Results/02_gmes/HT_genome_1st_run/genemark.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Results/02_gmes/Ht2.gff
gffParse.pl \
    -g Results/02_gmes/Ht2.gff \
    -i Results/01_clean/Ht.genome \
    -b HT \
    -d Results/03_gffParse \
    -p   \
    -c
```
gffParse.pl \
    -g Results/02_gmes/Ht2.gff \
    -i Results/01_clean/Ht.genome # \
    -b HT \  # basename for files
    -d Results/03_gffParse \  # output directory
    -p \ # also outputs amino acid files
    -c  #  adjusts genes if by changing readframe, STOP codons disappear


Run BLASTP
```bash
mkdir Results/04_blastp
nohup blastp \
 -query Results/03_gffParse/HT.faa \
 -db SwissProt -evalue 1e-10 `#database and evalue cut-off` \
 -out Results/04_blastp/Ht2.blastp -num_descriptions 10 \
 -num_alignments 5 -num_threads 4 &
jobs
```

link the database files
```bash
ln -s /resources/binp29/Data/malaria/taxonomy.dat Data/Taxonomy/taxonomy.dat
ln -s /resources/binp29/Data/malaria/uniprot_sprot.dat Data/Taxonomy/uniprot_sprot.dat
```

Use the blast taxonomy parser provided
```bash
python Scripts/datParser.py Results/04_blastp/Ht2.blastp Results/03_gffParse/HT.faa Data/Taxonomy/taxonomy.dat \
    Data/Taxonomy/uniprot_sprot.dat > Results/04_blastp/bird_scaffolds.txt
```
```bash
mkdir Results/05_filter
```
To remove the bird contigs
```bash
bird_contig=false
touch Results/05_filter/Ht_filtered.genome
while read line; do
  while read contig; do
    if [[ $line =~ $contig ]] ; then
#      echo "Found bird contig"
#      echo $line
      bird_contig=true
    fi
  done < Results/04_blastp/bird_scaffolds.txt
  
  if $bird_contig; then
    # it is a bird contig - do not save it
    read line
  else
    # it is a parasite contig
    echo $line 
    read line # read genome line
    echo $line
  fi
  bird_contig=false 
done < Results/01_clean/Ht.genome > Results/05_filter/Ht_filtered.genome
```
Number of contigs before and after filtering
```bash
cat Results/01_clean/Ht.genome | grep -c \>
cat Results/05_filter/Ht_filtered.genome | grep -c \>
```

    2343
    2330

#### make a new gene prediction

```bash
mkdir Results/06_gmes2
```

```bash
nohup gmes_petap.pl \
 --ES \
 --sequence Results/05_filter/Ht_filtered.genome \
 --cores 10 \
 --work_dir Results/06_gmes2/ \
 --min_contig 3000 & 
```
Now we make gff and parse
```bash
cat Results/06_gmes2/genemark.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Results/06_gmes2/Ht_filtered.gff
gffParse.pl \
    -g Results/06_gmes2/Ht_filtered.gff \
    -i Results/05_filter/Ht_filtered.genome \
    -b Ht_filter \
    -d Results/07_gffparse2 \
    -p   \
    -c
```
## gff parsing for the rest of the predictions
```bash
for file in Data/Predictions/*.g*; do

    species=$(basename "${file%.*}")
    cat $file | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Results/06_gmes2/${species}.gff
    gffParse.pl \
        -g $file \
        -i Data/Genomes/${species}.genome \
        -b "$species" \
        -d Results/07_gffparse2/${species} \
        -p   \
        -c
done
```

## genome information

```bash
# genome length and GC
for file in Data/Genomes/*.genome; do
  echo "$file"
  cat $file | awk '!/^>/ {tot_len+=length}{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} \
  END {printf "GC %.2f \t Genome length  %d \n", (gc*100)/(gc+at), (gc+at)}' 
done
```

Number of genes
```bash
find Results/07_gffparse2/ -name *.fna | 
  while read line; do
    echo $(basename "${line%.*}")
    cat $line | grep -c \>
  done
```

## proteinortho
```bash
mkdir Results/08_proteinortho
```

```bash
find Results/07_gffparse2/ -name *.faa |
  while read line; do
    cat $line | sed s'/\t.*//g' > ${line/.faa/_cleanheader.faa}
  done
nohup proteinortho6.pl -project=prot_orth $(find Results/07_gffparse2/ -name *header.faa)  &
mv prot_orth* Results/08_proteinortho/
```
Shared orthologues
```bash
cat Results/08_proteinortho/prot_orth.proteinortho.tsv |
  awk '($1==8) {count++} END {print "Shared orthologues " count}'
```
144


## busco
```bash
mkdir Results/09_busco
```
```bash
conda activate malaria_busco
find Results/07_gffparse2/ -name *.faa | 
  while read file; do
    species=$(basename "${file%.*}")
    busco -i $file -o Results/09_busco/$species -m prot -l apicomplexa
done
```

To check the results
```bash
touch Results/09_busco/busco_summary.txt 
find Results/09_busco -name short_summary.*.txt |
  while read file; do
    grep -A9 "Results:" $file >> Results/09_busco/busco_summary.txt
  done



```
