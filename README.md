Malaria case study
Author: Joan Escriva Font

# conda
```bash
conda create -n malaria 
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

To convert the gtf files to gff, gffParse.pl is used
```bash
mkdir Results/03_gffParse
gffParse.pl \
 -i Results/01_clean/Ht.genome \  # input fasta
 -g Results/02_gmes/HT_genome_1st_run \  # input annotation file
 -b HT \  # basename for files
 -d Results/03_gffParse \  # output directory
 -p \  # also input amino acid files
 -c  #  adjusts genes if by changing readframe, STOP codons disappear
```
