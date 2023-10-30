mkdir -p data/008/salmon/ketogenic_diet/raw_fastq
cd data/008/salmon/ketogenic_diet/raw_fastq

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/003/SRR19764203/SRR19764203.fastq.gz -o SRR19764203_GSM6257088_Ketogenic_Diet_Oligodendrocytes_-_R4_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/008/SRR19764208/SRR19764208.fastq.gz -o SRR19764208_GSM6257083_Standard_Diet_Oligodendrocytes_-_R4_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/010/SRR19764210/SRR19764210.fastq.gz -o SRR19764210_GSM6257081_Standard_Diet_Oligodendrocytes_-_R2_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/007/SRR19764207/SRR19764207.fastq.gz -o SRR19764207_GSM6257084_Standard_Diet_Oligodendrocytes_-_R5_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/009/SRR19764209/SRR19764209.fastq.gz -o SRR19764209_GSM6257082_Standard_Diet_Oligodendrocytes_-_R3_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/004/SRR19764204/SRR19764204.fastq.gz -o SRR19764204_GSM6257087_Ketogenic_Diet_Oligodendrocytes_-_R3_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/002/SRR19764202/SRR19764202.fastq.gz -o SRR19764202_GSM6257089_Ketogenic_Diet_Oligodendrocytes_-_R5_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/006/SRR19764206/SRR19764206.fastq.gz -o SRR19764206_GSM6257085_Ketogenic_Diet_Oligodendrocytes_-_R1_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/011/SRR19764211/SRR19764211.fastq.gz -o SRR19764211_GSM6257080_Standard_Diet_Oligodendrocytes_-_R1_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/005/SRR19764205/SRR19764205.fastq.gz -o SRR19764205_GSM6257086_Ketogenic_Diet_Oligodendrocytes_-_R2_Mus_musculus_RNA-Seq.fastq.gz

# You need an environment with salmon installed, if you don't have it, run:
# conda create -n salmon -c bioconda salmon
conda activate salmon 

# Always check if the gt file is supported by tximeta (https://bioconductor.org/packages/release/bioc/vignettes/tximeta/inst/doc/tximeta.html#Pre-computed_checksums) 
index="data/salmon/agc1_index"

for rname in data/008/ketogenic_diet/raw_fastq/*.fastq.gz
do

base=`basename $rname | sed -E 's/.*_([A-Za-z]+)_Diet_[^-]+_-_([^_]+)_.*\.fastq\.gz/\1_\2/'`
echo "Quantifying $base"
salmon quant -i $index \
-l A \
-p 12 \
-r $rname \
--validateMappings \
--gcBias \
--seqBias \
-o data/008/salmon/${base}_quant

done

conda deactivate