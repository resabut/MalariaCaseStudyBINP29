#!/usr/bin/env bash
mkdir Results/09_busco/temp
rm Results/09_busco/temp/*
grep -v "^#" Results/09_busco/Ht_filter/run_apicomplexa_odb10/full_table.tsv | awk '($2=="Complete") {print $1 "\n"}' | sort > Results/09_busco/temp/temp_0
first=Results/09_busco/temp/temp_0
# awk '!/^#/ {print $1}' $first | sort > Results/09_busco/temp
# wc -l Results/09_busco/temp
counter=0
find Results/09_busco/ -name full_table.tsv | sort -R |
  while read file; do
      if [[ $file != Results/09_busco/Toxoplasma_gondii/run_apicomplexa_odb10/full_table.tsv ]]; then
        counter=$((counter+1))
        comm -12 <(sort $first) <(grep -v "^#" $file | awk ' ($2=="Complete"){print $1}' | sort) > Results/09_busco/temp/temp_$counter
        echo "Comparing" $first "and" $file "and writing to temp_$counter"
        first=Results/09_busco/temp/temp_$counter
        wc -l Results/09_busco/temp/temp_$counter
        wc  -l $file
      fi
done
# take the latest temp file
final=Results/09_busco/temp/temp_$(ls Results/09_busco/temp/* | cut -f 3 -d "_" | sort -rn | head -1)


mkdir Results/09_busco/fasta_genes
rm Results/09_busco/fasta_genes/* # overwrite the existing files
# for every common BUSCO ID, get the gene ID and retrieve the gene sequence
for busco_id in $(cat $final); do
  # echo $busco_id
  # echo $counter
  # counter+=1
  # for every species
  find Results/07_gffparse2/ -name "*.faa" | grep -v header | grep -v Toxo |
    while read sp; do
        # get the species name
      species=$(basename "${sp%.*}")
      # echo $species $sp
        # parse the hmmer output files to get gene id
        gene_id=$(cat Results/09_busco/$species/run_apicomplexa_odb10/hmmer_output/initial_run_results/$busco_id.out |
                    grep -v "^#"| tail -1 | cut -f 1 -d " ")
        # echo $gene_id
        # retrieve the gene sequences for each faa

        grep -w -A 1 "$gene_id" $sp | sed "s/>.*/>$species/" | sed 's/-//g' >> Results/09_busco/fasta_genes/$busco_id.faa

  done
done





