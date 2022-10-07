# A tutorial for reprocessing data deposited on GEO
The SHARE-seqV2 alignment pipeline generates a pair of fastqs for each sample. These samples are ready to be submitted to GEO. Once the the fastqs are processed by GEO and they will be converted to SRA file with modified read header. This turorial demonstrates how to download deposited data and run the share-seqV2 pipeline.

We have deposited SHARE-seqV2 data on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207308)
1. As an example, we can download the species mixing ATAC SRA record by fasterq-dump -A SRR19912835 --split-files  -p
This will generate two fastq files: SRR19912835_1.fastq and SRR19912835_2.fastq

The read header looks like:\
head -4 SRR19912835_1.fastq\
@SRR19912835.1 A01389:111:H2Y5KDMXY:1:1101:1127:1000_R1.003,R2.032,R3.081,P1.06 length=50\
GGGCTACACAGAGAAACCCTGTCTCGAAAAACAAACAAAACAAAACAAAA\
+SRR19912835.1 A01389:111:H2Y5KDMXY:1:1101:1127:1000_R1.003,R2.032,R3.081,P1.06 length=50\
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:

2. Run these lines to modify the header and convert the format to be compatable with SHARE-seq pipeline. \
cat SRR19912835_1.fastq | awk '{if(NR%4==1) print "@"$2; else if(NR%4==2) print; else if(NR%4==3) print "+"; else if(NR%4==3) print $0}' | bgzip > speciesmix.ATAC.R1.fastq.gz\
cat SRR19912835_2.fastq | awk '{if(NR%4==1) print "@"$2; else if(NR%4==2) print; else if(NR%4==3) print "+"; else if(NR%4==3) print $0}' | bgzip > speciesmix.ATAC.R2.fastq.gz

The new read looks like:\
zless speciesmix.ATAC.R1.fastq.gz | head -4
@A01389:111:H2Y5KDMXY:1:1101:1127:1000_R1.003,R2.032,R3.081,P1.06\
GGGCTACACAGAGAAACCCTGTCTCGAAAAACAAACAAAACAAAACAAAA\
+\
@A01389:111:H2Y5KDMXY:1:1101:1832:1000_R1.083,R2.003,R3.037,P1.06

3. Move the resulting fastq.gz files into a new directoy (/fakepath/fastqfoler/).

4. Update the rawdir and Project in Share_seqV2_example.sh accordingly:\
rawdir=/fakepath/fastqfoler/ \
Project=(speciesmix.ATAC) \
Start=Demultiplexed # this optine makes the pipepline awares the starting file is the demultiplex fastqs.

5. Update the project name in config.example.yaml. 
6. Run the Share_seqV2_example.sh
