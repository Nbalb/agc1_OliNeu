# install majiq folowing their guide: 
# https://biociphers.bitbucket.io/majiq-docs-academic/getting-started-guide/installing.html
data=/mnt/d/OneDrive/PhD/projects/agc1_OliNeu/data
raw_reads=$data/delivery_20200428/raw_reads
cd $data
mkdir data/isoforms

# Align reads
mkdir -p isoforms/idx/
curl -L ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz \
-o isoforms/genome.fa.gz
gzip -d isoforms/genome.fa.gz
hisat2-build -p 6 isoforms/genome.fa isoforms/idx/mm10

trimmed_fastq=$data/isoforms/trimmed_fastq
bams=$data/isoforms/bams
mkdir $trimmed_fastq
mkdir $bams

for i in {1..6}
do

r1=${raw_reads}/1_ID1603_$i-Olineu${i}_S1_L001_R1_001.fastq.gz
r2=${raw_reads}/1_ID1603_$i-Olineu${i}_S1_L001_R2_001.fastq.gz
out_trim1=$trimmed_fastq/Olineu_${i}_R1_trimmed.fastq.gz
out_trim2=$trimmed_fastq/Olineu_${i}_R2_trimmed.fastq.gz
report=$trimmed_fastq/fastp_sample${i}.html
echo "Trimming sample $i"
fastp -i $r1 -I $r2 -o $out_trim1 -O $out_trim2 -h $report -w 6

out_bam=$bams/Olineu_${i}.bam
sorted_bam=$bams/Olineu${i}.sorted.bam

echo "Aligning sample $i"
hisat2 -x isoforms/idx/mm10 -p 6 -1 $out_trim1 -2 $out_trim2 | samtools view -bS - > $out_bam

echo "Sorting sample $i"
samtools sort -m 3G -@ 6 -O BAM -o $sorted_bam $out_bam && rm $out_bam

echo "Indexing sample $i"
samtools index $sorted_bam
done


### Reads are forward-stranded
### Check settings.ini inside isoforms folder to verify settings and run builder
### Create Python environment and install htslib
curl -OL https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2
tar -xf htslib-1.13.tar.bz2  # extract archive
cd htslib-1.13  # change into htslib source directory
# configure, make, and install htslib to ~/install/htslib-1.13
./configure --prefix=$HOME/install/htslib-1.13
make
make install

cd $data/isoforms
curl -OL https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gff3.gz
mkdir $data/isoforms/build
python3 -m venv majiq_env
source majiq_env/bin/activate

majiq build gencode.vM25.annotation.gff3.gz -c settings.ini -j 6 -o $data/isoforms/build

### PSI analysis
cd $data/isoforms/build
majiq psi Olineu1.sorted.majiq Olineu2.sorted.majiq Olineu3.sorted.majiq -j 6 -o $data/isoforms/psi/ -n wt
majiq psi Olineu4.sorted.majiq Olineu5.sorted.majiq Olineu6.sorted.majiq -j 6 -o $data/isoforms/psi/ -n kd

### DeltaPSI analysis
majiq deltapsi \
-grp1 Olineu1.sorted.majiq Olineu2.sorted.majiq Olineu3.sorted.majiq \
-grp2 Olineu4.sorted.majiq Olineu5.sorted.majiq Olineu6.sorted.majiq \
-j 6 -o $data/isoforms/psi/ -n wt kd

### Voila visualization
voila tsv splicegraph.sql $data/isoforms/psi/wt-kd.deltapsi.voila -f voila_results
voila view -p 59408 splicegraph.sql $data/isoforms/psi/wt-kd.deltapsi.voila   

