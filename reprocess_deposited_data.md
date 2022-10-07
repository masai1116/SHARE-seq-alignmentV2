# A tutorial for reprocessing data deposited on GEO
The SHARE-seqV2 alignment pipeline generates a pair of fastqs for each sample. These samples are ready to be submitted to GEO. Once the the fastqs are processed by GEO and they will be converted to SRA file with modified read header. This turorial demonstrates how to download deposited data and run the share-seqV2 pipeline.

We have deposited SHARE-seqV2 data on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207308)
As an example, we can download the species mixing ATAC SRA record by fasterq-dump -A SRR19912835 --split-files  -p
This will generate two fastq files: SRR19912835_1.fastq and SRR19912835_2.fastq
The read header looks like:\
head -4 SRR19912835_1.fastq\
@SRR19912835.1 A01389:111:H2Y5KDMXY:1:1101:1127:1000_R1.003,R2.032,R3.081,P1.06 length=50\
GGGCTACACAGAGAAACCCTGTCTCGAAAAACAAACAAAACAAAACAAAA\
+SRR19912835.1 A01389:111:H2Y5KDMXY:1:1101:1127:1000_R1.003,R2.032,R3.081,P1.06 length=50\
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:\

