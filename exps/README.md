```
wget http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz -c > chr21.fa

wget http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz
zgrep -P "^21\t" Homo_sapiens.GRCh38.107.gtf.gz | grep "protein_coding" > chr21.gtf

/usr/bin/time -v python3 ~/code/eflux/eflux.py -t 4 chr21.fa chr21.gtf


mamba create -c bioconda -n eflux-exps samtools seaborn biopython gffutils suppa openjdk=8.0.332 flux-simulator star rmats snakemake-minimal

```