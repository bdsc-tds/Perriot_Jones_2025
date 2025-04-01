raw_data=../data/raw/GEX/bam/
out_dir=../data/processed/kir_genotype

for sample in `ls ${raw_data}`
do
  echo ${sample}
  run-t1k -b ${raw_data}${sample}/possorted_genome_bam.bam \
  -c ../kiridx/kiridx_rna_coord.fa \
  -o ${sample} \
  --od ${out_dir} \
  -f ../kiridx/kiridx_rna_seq.fa \
  --barcode CR
done


