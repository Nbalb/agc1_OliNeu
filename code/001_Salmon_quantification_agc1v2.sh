mkdir -p agc1_OliNeu/data/salmon
cd data/salmon
curl -OL http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
cd ../../
# You need an environment with salmon installed, if you don't have it, run:
# conda create -n salmon -c bioconda salmon
conda activate salmon 
salmon index -t Mus_musculus.GRCm39.cdna.all.fa.gz -i agc1_index -k 31

# Always check if the gt file is supported by tximeta (https://bioconductor.org/packages/release/bioc/vignettes/tximeta/inst/doc/tximeta.html#Pre-computed_checksums) 
index="data/salmon/agc1_index"

for rname in data/fastq/*[1,2,3]_R1.fastq.gz
do

mkdir fastq
# add lines to download fastq
cd fastq
base=`basename $rname _R1.fastq.gz`
echo "Quantifying $base"
salmon quant -i $index \
-l A \
-p 5 \
-1 ${base}_R1.fastq.gz \
-2 ${base}_R2.fastq.gz \
--validateMappings \
--gcBias \
--seqBias \
-o ${base}_quant

done

conda deactivate
