#$ -S /bin/sh
if [ $# -eq 0 ]; then
        echo "Program: Sequence Alignment with Bowtie2"
        echo "Usage: $0 <genome.fasta> <fq1> <fq2>"
        exit 1
fi


ref=$1
fq1=$2
fq2=$3
echo "#$ -S /bin/sh" > bowtie.sh
echo "/hellogene/scgene01/bio/software/bowtie2-2.2.5/bowtie2-build ${ref} seq" >> bowtie.sh
echo "/hellogene/scgene01/bio/software/bowtie2-2.2.5/bowtie2 -x seq -1  ${fq1} -2 ${fq2} -p 25 | /hellogene/scgene01/bio/bin/samtools view -SbF 0x804 - | /hellogene/scgene01/bio/bin/samtools sort -@ 25 - > seq.sort.bam " >> bowtie.sh
echo "/hellogene/scgene01/bio/bin/samtools index seq.sort.bam"  >> bowtie.sh
qsub -cwd -l vf=20g,cpu=32 bowtie.sh

