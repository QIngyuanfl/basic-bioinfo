echo "/hellogene/scgene01/bio/bin/bwa index $1" >bwa.sh
echo "/hellogene/scgene01/bio/bin/bwa mem -t 14 $1 $2 $3 | /hellogene/scgene01/bio/bin/samtools view -SbF 0x804 - | /hellogene/scgene01/bio/bin/samtools sort -@ 8 -  > genome.map.bam" >>bwa.sh
echo "/hellogene/scgene01/bio/bin/samtools stats genome.map.bam | grep ^COV | cut -f 2- > map.cov" >>bwa.sh
echo "/hellogene/scgene01/bio/bin/samtools depth genome.map.bam > map.dp" >>bwa.sh
echo "/hellogene/scgene01/bio/bin/samtools index genome.map.bam" >>bwa.sh
qsub -cwd -l vf=14g,cpu=32 bwa.sh
