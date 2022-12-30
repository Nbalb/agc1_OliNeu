#mamba create -n majiq_env -c http://majiq.biociphers.org/download/channel/linux-64 majiq
cd mnt/d/onedrive/linux/agc1

### Following (https://biociphers.bitbucket.io/majiq/quick.html#conf_file)
### Check if reads are stranded/rev_stranded/unstranded
### Use paramter -s 0(unstranded), 1(stranded), 2(rev stranded) and check the summary
featureCounts -T 6 -p -s 0 -a gencode.vM25.annotation.gff3 -o agc1feature_unstrand.counts /mnt/d/onedrive/linux/agc1/fastq/wt1.sorted.bam
featureCounts -T 6 -p -s 1 -a gencode.vM25.annotation.gff3 -o agc1feature_stranded.counts /mnt/d/onedrive/linux/agc1/fastq/wt1.sorted.bam
featureCounts -T 6 -p -s 2 -a gencode.vM25.annotation.gff3 -o agc1feature_rev_strand.counts /mnt/d/onedrive/linux/agc1/fastq/wt1.sorted.bam

### Reads are stranded
### Check settings.ini inside isoforms folder to verify settings and run builder
cd
source env/bin/activate
cd $linux/agc1/isoforms/
majiq build ../gencode.vM25.annotation.gff3 -c settings.ini -j 4 -o .

### PSI analysis
majiq psi wt1.sorted.majiq wt2.sorted.majiq wt3.sorted.majiq -j 4 -o psi/ -n wt
majiq psi kd1.sorted.majiq kd2.sorted.majiq kd3.sorted.majiq -j 4 -o psi/ -n kd

### DeltaPSI analysis
majiq deltapsi -grp1 wt1.sorted.majiq wt2.sorted.majiq wt3.sorted.majiq -grp2 kd1.sorted.majiq kd2.sorted.majiq kd3.sorted.majiq -j 4 -o psi/ -n wt kd

### Voila visualization
voila tsv splicegraph.sql psi/wt_kd.deltapsi.voila -f voila_results
voila view -p 59408 splicegraph.sql psi/wt_kd.deltapsi.voila   

