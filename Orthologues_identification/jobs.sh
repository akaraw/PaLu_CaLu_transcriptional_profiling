python longest_isoform.py palecto.faa > palecto_long.faa #NCBI formatted prteome fasta
python longest_isoform.py human.faa > human_long.faa #NCBI formatted proteome fasta

mmseqs createdb palecto_long.faa batdb
mmseqs createdb human_long.faa mandb

mmseqs search batdb mandb bat_man_hits tmp --max-seqs 1000 --min-seq-id 0.5
mmseqs convertalis  batdb mandb bat_man_hits orthologs.m8 --format-output "query,target,evalue,pident"
cat orthologs.m8 | awk '$3 < 1e-05 {print $1"\t"$2}' > orthologs.tsv
python deduplicate.py orthologs.tsv final_orthologs.tsv
python deduplicate2.py final_orthologs.tsv orthololist_human_vs_bat.tsv

python replace_values_bat.py ../RNAseq_pipeline/palecto/htseq/ orthololist_human_vs_bat.tsv
python replace_values_human.py ../RNAseq_pipeline/human/htseq/ ../RNAseq_pipeline/palecto/htseq/T5ortho.tsv
