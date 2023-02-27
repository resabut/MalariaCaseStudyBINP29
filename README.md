Malaria case study
Author: Joan Escriva Font

# conda
```bash
conda create -n malaria fastqc=
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
```
    Plasmodium_berghei.genome
    Plasmodium_cynomolgi.genome
    Plasmodium_faciparum.genome
    Plasmodium_knowlesi.genome
    Plasmodium_vivax.genome
    Plasmodium_yoelii.genome
    Toxoplasma_gondii.genome

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
scp /tmp/Prediction/*.gtf /Data/Predictions
```


## Processing of *Haemoproteus tartakovskyi* data
### Clean genome sequence
The GC threshold should be 35 for more astringency.

    Bird: 41 %  
    Parasite: 19-42%  

