#!/bin/bash
# Author: Sai Ma <sai@broadinstitute.org>
# Last modified date: 2021/05/25
# Designed for processing share-atac/rna/taps/dipC/cite/cellhashing/crop

# input file
# 1) yaml
# 2) BCL or fastq
# when there are 4 fastqs per lane, fastqs need to be named as "*S1_L001/2/3/4_R1/R2/I1/I2_001.fastq.gz" (when there are multiple lanes) or "_S1_R1/R2/I1/I2_001.fastq.gz" (when there is one lane), R1: bio read1, R2: bio read2, I1: index1, I2: index2
# when there are 2 fastqs per lane, fastqs need to be named as "*S1_L001/2/3/4_R1/2_001.fastq.gz" or "_S1_R1/2_001.fastq.gz". important to rename files as R1 and "R4"

rawdir=/mnt/users/sai/Script/Split-seq_Sai/example_fastq/
dir=/mnt/users/sai/Script/Split-seq_Sai/example_output/
yaml=/mnt/users/sai/Script/Split-seq_Sai/config.example.yaml
Project=(BMMC.RNA BMMC.ATAC)

Type=(RNA ATAC)
# ATAC RNA notrim (keep all information, in case of any spatial design)
# fill the same information in yaml file

Genomes=(hg19 hg19)
# both (hg19+mm10) mm10 hg19 hg38 noalign (skip alignment)
ReadsPerBarcode=(10 10)
# reads cutoff to barcodes: 100 for full run; 10 for QC run

Start=Fastq # Bcl or Fastq
Runtype=QC # QC or full,  QC only analyze 12M reads
chem=fwd # rev or fwd: nova 1.5 & nextseq use rev; nova1.0 uses fwd

###########################
## RNA-seq options, usually don't change these options
removeSingelReadUMI=F # default F; T if seq deep enough or slant is a big concern; def use F for crop-seq
keepIntron=T #default T
cores=16
genename=gene_name # default gene_name; gene_name (official gene symbol) or gene_id (ensemble gene name)
refgene=gencode # default gencode; gencode or genes; genes is UCSC genes; gencode also annotate ncRNA
mode=fast # fast or regular; default fast; fast: dedup with custom  script; regular: dedup with umitools
# fast mode gives more UMIs because taking genome position into account when dedup
# fast mode doesn't collapse UMIs map to different position
# in fast mode, the lib size estimation is not accurate

keepMultiMapping=(T T T T T T T T T T T T T T T T)
# default F; F for species mixing or cell lines, T for low yield tissues (only keep the primarily aligned reads), doesn't matter for ATAC
# allowing multi-mapping redas will increase the percent of mito reads

keepmito=(F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F)
## default F, remove mito reads

cleanup=F
# clean up fastq or not

#### do not change following

# define path
toolPATH='/mnt/Apps/JDB_tools/'
myPATH='/mnt/users/sai/Script/Split-seq_Sai/'
tssFilesPATH='/mnt/users/sai/Script/Split-seq_Sai/TSSfiles/' 
picardPATH='/mnt/bin/picard/picard.jar'
bowtieGenome='/mnt/users/sai/Script/Split-seq_Sai/refGenome/bowtie2/'
starGenome='/mnt/users/sai/Data/star/genome/'
genomeBed='/mnt/users/sai/Script/Split-seq_Sai/genomeBed/'
bismarkGenome='/mnt/users/sai/Script/Split-seq_Sai/refGenome/bismark/'
fastpPath='/mnt/users/sai/Package/fastp/'

# real code
export SHELL=$(type -p bash)
source ~/.bashrc
export LC_COLLATE=C
export LANG=C

echo "the number of projects is" ${#Project[@]}
echo "Running $Runtype pipeline"

if [ ! -d $dir ]; then mkdir $dir; fi
if [ ! -d $dir/fastqs ]; then mkdir $dir/fastqs ; fi
if [ ! -d $dir/temp ]; then mkdir $dir/temp ; fi

cp $myPATH/Split*.sh $dir/
cp $yaml $dir/
cd $dir 
if [ -f $dir/Run.log ]; then rm $dir/Run.log; fi

export PATH="/mnt/users/sai/miniconda2/bin:$PATH"

# check yaml file, make sure no duplicated P1.xx
p1=$(cat $yaml | grep P1. | sort | uniq -d)
if [ -z "$p1" ]; then
    # empty
    echo "Yaml looks good"
else
    echo "Yaml contains duplicated P1.xx, exiting..."
    exit
fi

# if start with bcl file
if [ "$Start" = Bcl ]; then
    echo "Bcl2fastq"
    if [ -f $rawdir/fastqs/Undetermined_S0_R1_001.fastq.gz ] || [ -f $rawdir/fastqs/Undetermined_S1_R1_001.fastq.gz ] || [ -f $rawdir/fastqs/Undetermined_S0_L001_R1_001.fastq.gz ]; then
	echo "Found Undetermined_S0_L001_I1_001.fastq.gz, skip Bcl2Fastq"
    else
	echo "Converting bcl to fastq"
	mkdir $rawdir/fastqs/
	bcl2fastq -p $cores -R $rawdir --mask-short-adapter-reads 0 -o $rawdir/fastqs/ --create-fastq-for-index-reads  2>>$dir/Run.log
	cd $rawdir/fastqs/    
    fi
    if [ -f $dir/fastqs/filesi1.xls ]; then rm $dir/fastqs/filler*; fi
    
    rawdir=$rawdir/fastqs/
    Start=Fastq
fi


if [ "$Start" = Fastq ]; then
    echo "Skip bcltofastq"
    if [ -f $dir/fastqs/filesi1.xls ]; then rm $dir/fastqs/filler*; fi
    cd $rawdir
    if ls *L001_R1_001.fastq.gz 1> /dev/null 2>&1; then
	temp=$(ls *_L001_R1_001.fastq.gz)
	Run=$(echo $temp | sed -e 's/\_S0\_L001\_R1\_001.fastq.gz//')
	singlelane=F
	temp=$(ls *L00*_R1_001.fastq.gz)
	VAR=( $temp )
	nolane=${#VAR[@]}
	echo "Detected $nolane lanes"
    elif ls *S1_R1_001.fastq.gz 1> /dev/null 2>&1; then
	echo "Detected single lane"
	temp=$(ls *S1_R1_001.fastq.gz)
	Run=$(echo $temp | sed -e 's/\_\S1\_\R1\_\001.fastq.gz//')
	singlelane=T
	nolane=1
    else
	echo "No fastq with matched naming format detected; exit..."
	exit
    fi
    
    echo "Run number is:" $Run

    # split fastqs
    mkdir $dir/smallfastqs/
    if [ -f $dir/smallfastqs/0001.1.$Run.R2.fastq.gz ]; then
        echo "Found 0001.$Run.R2.fastq, skip split fastqs"
    else
	if [ ! -f $dir/R1.1.fastq.gz ]; then
            echo "Link fastqs"
	    if [ $singlelane == T ]; then
		ln -s $rawdir/"$Run"_S1_R1_001.fastq.gz $dir/R1.1.fastq.gz
		ln -s $rawdir/"$Run"_S1_R2_001.fastq.gz $dir/R2.1.fastq.gz
		ln -s $rawdir/"$Run"_S1_I1_001.fastq.gz $dir/I1.1.fastq.gz
		ln -s $rawdir/"$Run"_S1_I2_001.fastq.gz $dir/I2.1.fastq.gz
	    else
		parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_R1_001.fastq.gz '$dir'/R1.{}.fastq.gz' ::: $(seq $nolane)
		parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_R2_001.fastq.gz '$dir'/R2.{}.fastq.gz' ::: $(seq $nolane)
		parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_I1_001.fastq.gz '$dir'/I1.{}.fastq.gz' ::: $(seq $nolane)
		parallel 'ln -s '$rawdir'/'$Run'_S0_L00{}_I2_001.fastq.gz '$dir'/I2.{}.fastq.gz' ::: $(seq $nolane)
	    fi
	fi
	
	if [ "$Runtype" = full ]; then
	    # Runing full pipeline
	    echo "Split fastqs to small files"
	    dosplitfull(){
                /mnt/users/sai/Package/fastp/fastp -i $2/R1.$1.fastq.gz -o $2/smallfastqs/$1.$3.R1.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                /mnt/users/sai/Package/fastp/fastp -i $2/R2.$1.fastq.gz -o $2/smallfastqs/$1.$3.R2.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                /mnt/users/sai/Package/fastp/fastp -i $2/I1.$1.fastq.gz -o $2/smallfastqs/$1.$3.I1.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
                /mnt/users/sai/Package/fastp/fastp -i $2/I2.$1.fastq.gz -o $2/smallfastqs/$1.$3.I2.fastq.gz -S 80000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		wait
            }
            export -f dosplitfull
            parallel --delay 1 dosplitfull {} $dir $Run ::: $(seq $nolane)
	elif [ "$Runtype" = QC ]; then
	    # Runing QC pipeline
	    echo "Split fastqs to small files"
	    dosplitQC(){
		/mnt/users/sai/Package/fastp/fastp -i $2/R1.$1.fastq.gz -o $2/smallfastqs/$1.$3.R1.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		/mnt/users/sai/Package/fastp/fastp -i $2/R2.$1.fastq.gz -o $2/smallfastqs/$1.$3.R2.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		/mnt/users/sai/Package/fastp/fastp -i $2/I1.$1.fastq.gz -o $2/smallfastqs/$1.$3.I1.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		/mnt/users/sai/Package/fastp/fastp -i $2/I2.$1.fastq.gz -o $2/smallfastqs/$1.$3.I2.fastq.gz --reads_to_process $4 -S 4000000 --thread 1 -d 4 -A -G -L -Q 2>>$2/split.log &
		wait
            }
            export -f dosplitQC
	    let reads=12100000/$nolane
	    parallel --delay 1 dosplitQC {} $dir $Run $reads ::: $(seq $nolane)
	else
            echo "Unknown sequencer type, exiting" && exit
        fi
    fi
    
    # trim and index fastq
    if [ -f $dir/fastp.json ]; then rm $dir/fastp.json $dir/fastp.html; fi
    ls $dir/smallfastqs | grep R1 > $dir/filesr1.xls
    ls $dir/smallfastqs | grep R2 > $dir/filesr2.xls
    ls $dir/smallfastqs | grep I1 > $dir/filesi1.xls
    ls $dir/smallfastqs | grep I2 > $dir/filesi2.xls
    cd $dir/
    if [ -f $dir/fastqs/Sub.0001.1.discard.R1.fq.gz ]; then
        echo "Found Sub.0001.1.discard.R1.fq.gz, skip updating index"
    else
        echo "Update index and trim fastqs"
	noreadfile=`ls $dir/smallfastqs | grep R1 2>/dev/null | wc -l`
	noindexfile=`ls $dir/smallfastqs | grep I1 2>/dev/null | wc -l`
	if [ $noreadfile == $noindexfile ]; then
	    paste filesr1.xls filesr2.xls filesi1.xls filesi2.xls | awk -v OFS='\t' '{print $1, $2, $3, $4, substr($1,1,7)}'> Filelist2.xls
	    parallel --jobs $cores --colsep '\t' 'if [ -f '$dir'/fastqs/Sub.{5}.discard.R1.fq.gz ]; then echo "found Sub.{5}.discard.R1.fq.gz"; \
	    	     	    	   	  else python3 '$myPATH'/fastq.process.py3.v0.8.py \
                                          -a '$dir'/smallfastqs/{1} -b '$dir'/smallfastqs/{2} \
                                          --c '$dir'/smallfastqs/{3} --d '$dir'/smallfastqs/{4} \
					  --out '$dir'/fastqs/Sub.{5} \
					  -t '$chem' -y '$yaml' && pigz --fast -p 4 '$dir'/fastqs/Sub.{5}*fq; fi' :::: Filelist2.xls
	else
	    paste filesr1.xls filesr2.xls | awk -v OFS='\t' '{print $1, $2, substr($1,1,7)}'> Filelist2.xls
	    parallel --jobs $cores --colsep '\t' 'if [ -f '$dir'/fastqs/Sub.{3}.discard.R1.fq.gz ]; then echo "found Sub.{3}.discard.R1.fq.gz"; \ 
	    	     	    	   	  else python3 '$myPATH'/fastq.process.py3.v0.8.py -a '$dir'/smallfastqs/{1} -b '$dir'/smallfastqs/{2} \
                                          --out '$dir'/fastqs/Sub.{3} \
					  -t '$chem' -y '$yaml' && pigz --fast -p 4 '$dir'/fastqs/Sub.{3}*fq; fi' :::: Filelist2.xls
	    # --qc # option is available to process even smaller number of reads
	fi
    fi
    rm filesr1.xls filesr2.xls filesi1.xls filesi2.xls
    if [ -f Filelist2.xls ]; then
	rm Filelist2.xls
    fi
fi 

# merge fastq
echo "Merge fastqs"
parallel --jobs 4 'if [ -f '$dir'/fastqs/{}.R1.fastq.gz ] || [ -d '$dir'/{} ]; then echo "found {}.R1.fastq.gz or {} folder"; \
	 	      	   else ls '$dir'/fastqs/Sub*{}*R1.fq.gz | xargs cat > '$dir'/fastqs/{}.R1.fastq.gz && echo "Generated {}.R1.fastq.gz"; fi' ::: ${Project[@]} discard 
parallel --jobs 4 'if [ -f '$dir'/fastqs/{}.R2.fastq.gz ] || [ -d '$dir'/{} ]; then echo "found {}.R2.fastq.gz or {} folder"; \
                           else ls '$dir'/fastqs/Sub*{}*R2.fq.gz | xargs cat > '$dir'/fastqs/{}.R2.fastq.gz && echo "Generated {}.R2.fastq.gz"; fi' ::: ${Project[@]} discard

## clean up small fastqs
if [ $cleanup == "T" ]; then
    rm -r $dir/smallfastqs/* $dir/*/Sub.*.fq.gz
    touch $dir/smallfastqs/0001.1.$Run.R2.fastq.gz
    touch $dir/fastqs/Sub.0001.1.discard.R1.fq.gz
    rm $dir/fastqs/discard.R1.fastq.gz
    rm $dir/fastqs/discard.R2.fastq.gz
    touch $dir/fastqs/discard.R1.fastq.gz
    touch $dir/fastqs/discard.R2.fastq.gz
fi

# align
index=0
for Name in ${Project[@]}; do
    echo "project $index : $Name"
    
    if [ ${Type[$index]} == "notrim" ]; then
	if [ ! -d $dir/$Name ]; then
	    mkdir $dir/$Name
	fi
    fi
    if [ -d $dir/$Name ]; then
	echo "Found $Name dir, skip this project"
    else
	if [ ${Genomes[$index]} == "noalign" ]; then
	    mkdir $dir/$Name
	    mv $dir/fastq/$Name*.fastq.gz $dir/$Name
	    ## this will skip aligning this project
	fi
	if [ ${Genomes[$index]} == "both" ]; then
	    if [ ${Type[$index]} == "ATAC" ]; then
		Genome1=(hg19 mm10) # Genome1 for alignment, Genome2 for other steps
	    else
		Genome1=(both) # for RNA and cellhash
	    fi
            Genome2=(hg19 mm10)
	else
	    Genome1=${Genomes[$index]}
	    Genome2=${Genomes[$index]}
	fi
	
	for Species in ${Genome1[@]}; do
	    if [ -d $dir/$Name ]; then
		echo "Found $Name folder, skip alignment"
	    else
		echo "Align $Name to $Species ${Type[$index]} library"
		cd $dir/fastqs/
		if [ -f $dir/fastqs/$Name.$Species.bam ] || [ -f $dir/fastqs/$Name.$Species.st.bam ]; then
		    echo "Found $Name.$Species.bam, skip alignment"
		elif [ ${Type[$index]} == "ATAC" ]; then
		    (bowtie2 -X2000 -p $cores --rg-id $Name \
			     -x $bowtieGenome/$Species/$Species \
			     -1 $dir/fastqs/$Name.R1.fastq.gz \
			     -2 $dir/fastqs/$Name.R2.fastq.gz | \
			 samtools view -bS -@ $cores - -o $Name.$Species.bam) 2>$Name.$Species.align.log
		elif [ ${Type[$index]} == "RNA" ]; then
		    if [ -f $dir/fastqs/$Name.$Species.align.log ]; then
			echo "found $dir/fastqs/$Name.$Species.align.log, skip alignment for UMI reads"
		    else
			echo "Align UMI reads"
			STAR --chimOutType WithinBAM \
			     --runThreadN $cores \
			     --genomeDir $starGenome/$Species/ \
			     --readFilesIn $dir/fastqs/$Name.R1.fastq.gz  \
			     --outFileNamePrefix $dir/fastqs/$Name.$Species. \
			     --outFilterMultimapNmax 20 \
			     --outFilterScoreMinOverLread 0.3 \
			     --outFilterMatchNminOverLread 0.3 \
			     --outSAMattributes NH HI AS nM MD \
			     --limitOutSJcollapsed 5000000 \
			     --outSAMtype BAM Unsorted \
			     --limitIObufferSize 400000000 \
			     --outReadsUnmapped None \
			     --readFilesCommand zcat
			    
			mv $dir/fastqs/$Name.$Species.Aligned.out.bam $dir/fastqs/$Name.$Species.bam
			rm -r *_STARtmp *Log.progress.out *SJ.out.tab
			mv $dir/fastqs/$Name.$Species.Log.final.out $dir/fastqs/$Name.$Species.align.log
		    fi
		else
		    echo "Unknown parameter"
		    exit
		fi
		if [ -f $dir/fastqs/$Name.$Species.st.bam ]; then
		    echo "Found $Name.$Species.st.bam, skip sorting bam"
		else
		    echo "Sort $Name.$Species.bam"
		    cd $dir/fastqs/
		    samtools sort -@ 8 -m 8G $Name.$Species.bam > $Name.$Species.st.bam
		    samtools index -@ $cores $Name.$Species.st.bam
		    rm $Name.$Species.bam
		fi
		# exit		
		# Update RGID's and Headers
		if [ -f $dir/fastqs/$Name.$Species.rigid.reheader.st.bam.bai ]; then
		    echo "Found $Name.$Species.rigid.reheader.st.bam, skip update RGID"
		elif [ ${Type[$index]} == "RNA" ]; then
		    if [ -f $Name.$Species.rigid.reheader.unique.st.bam.bai ]; then
                        echo "Found $Name.$Species.rigid.reheader.st.bam, skip update RG tag"
		    else
			echo "Update RGID for $Name.$Species.st.bam"
			echo "Remove low quality reads and unwanted chrs"
			samtools view -H $Name.$Species.st.bam | sed 's/chrMT/chrM/g' > $Name.$Species.st.header.sam
			if [ ${keepMultiMapping[$index]} == "T" ] && [ ${Type[$index]} == "RNA" ]; then
			    # keep primary aligned reads only
			    # modify chrMT to chrM
			    cat $Name.$Species.st.header.sam <(samtools view -@ $cores $Name.$Species.st.bam | \
				   sed 's/chrMT/chrM/g') | \
				   samtools view -@ $cores -bS -F 256 > $Name.$Species.rigid.reheader.unique.st.bam
			else
                            cat $Name.$Species.st.header.sam <(samtools view -@ $cores $Name.$Species.st.bam | \
				   sed 's/chrMT/chrM/g')	| \
				   samtools view -@ $cores -bS -q 30  > $Name.$Species.rigid.reheader.unique.st.bam
			fi
			samtools index -@ $cores $Name.$Species.rigid.reheader.unique.st.bam
			rm $Name.$Species.st.header.sam
		    fi
		fi
		
		# ATAC processing convert bam to bed
		if [ ${Type[$index]} == "ATAC" ]; then
		    if [ -f $Name.$Species.rmdup.bam ]; then
			echo "Skip processing ATAC bam"
		    else
			if [ -f $Name.$Species.namesort.bam ]; then
		            echo "Skip sort bam on name"
			else
			    echo "Sort bam on name"
			    echo "Remove low quality reads, unwanted chrs & namesort" $Name.$Species.st.bam
			    if [ ${keepmito[$index]} == "F" ]; then
				chrs=`samtools view -H $Name.$Species.st.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v Y | awk '{if(length($0)<6)print}'`
			    else
				chrs=`samtools view -H $Name.$Species.st.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v Y | awk '{if(length($0)<6)print}'`
			    fi
			    samtools view -b -q 30 -f 0x2 $Name.$Species.st.bam  `echo $chrs` | samtools sort -@ $cores -m 3G -n -o $Name.$Species.namesort.bam
			fi
			if [ -f $Name.$Species.bed.gz ]; then
			    echo "Skip converting bam to bed.gz"
			else
                            echo "Convert namesort.bam to bed.gz & mark duplicates"
			    bedtools bamtobed -i $Name.$Species.namesort.bam -bedpe | \
				sed 's/_/\t/g' | \
				awk -v OFS="\t" '{if($10=="+"){print $1,$2+4,$6-5,$8}else if($10=="-"){print $1,$2-5,$6+4,$8}}' |\
				sort --parallel=$cores -S 40G  -k4,4 -k1,1 -k2,2n -k3,3n | \
				uniq -c | \
				awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | \
				pigz --fast -p $cores > $Name.$Species.bed.gz
			    rm $Name.$Species.namesort.bam
			    ## when the DNA is taged twice, only use R2 for dedup..
			fi
			# convert to a bam file for QC
                        bedToBam -i <(zcat $Name.$Species.bed.gz) -g $genomeBed/$Species.chrom.sizes | samtools sort -@ $cores -m 2G - > $Name.$Species.rmdup.bam
                        samtools index -@ $cores $Name.$Species.rmdup.bam
                    fi
		fi
				       
		# count reads for ATAC
		if [ ${Type[$index]} == "ATAC" ] ; then
		    if [ -f $Name.$Species.filtered.counts.csv ] || [ -f $Name.$Species.filtered.counts.csv.gz ]; then
			echo "found $Name.$Species.filtered.counts.csv"
		    else
			echo "count unfiltered reads"
			zcat $Name.$Species.bed.gz | \
			    awk -v OFS='\t' '{a[$4] += $5} END{for (i in a) print a[i], i}'| \
			    awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }'> $Name.$Species.wdup.RG.freq.bed
			/usr/local/bin/Rscript $myPATH/sum_reads_v2.R $dir/fastqs/ $Name.$Species.wdup.RG.freq.bed --save
			mv $Name.$Species.wdup.RG.freq.bed.csv $Name.$Species.unfiltered.counts.csv

			echo "count filtered reads"
			zcat $Name.$Species.bed.gz | cut -f4 | uniq -c | \
			    awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }' > $Name.$Species.rmdup.RG.freq.bed
			/usr/local/bin/Rscript $myPATH/sum_reads_v2.R $dir/fastqs/ $Name.$Species.rmdup.RG.freq.bed --save
			mv $Name.$Species.rmdup.RG.freq.bed.csv $Name.$Species.filtered.counts.csv

			rm $Name.$Species.wdup.RG.freq.bed $Name.$Species.rmdup.RG.freq.bed 
		    fi
		fi

		# remove barcode with low counts from the fragment file for ATAC
                if [ ${Type[$index]} == "ATAC" ] ; then
		    if [ -f $Name.$Species.fragments.tsv.gz ]; then
			echo "Skip removing low counts barcode combination"
                    else
			echo "Remove low counts barcode combination"
			sed -e 's/,/\t/g' $Name.$Species.filtered.counts.csv | \
                            awk -v OFS=',' 'NR>=2 {if($5 >= '${ReadsPerBarcode[$index]}') print $1,$2,$3,$4} '  > $Name.$Species.barcodes.txt
			grep -wFf $Name.$Species.barcodes.txt  <(zcat $Name.$Species.bed.gz) | \
			    sort --parallel=$cores -S 40G -k1,1 -k2,2n -k3,3n -k4,4 | bgzip  > $Name.$Species.fragments.tsv.gz
			tabix -p bed $Name.$Species.fragments.tsv.gz
		    fi
		fi
		
		# RNA processing
		if [ ${Type[$index]} == "RNA" ] && [ ! -f $Name.mm10.rigid.reheader.unique.st.bam ]; then
		    cd $dir/fastqs/
		    if [ ${Genomes[$index]} == "both" ]; then
			echo "Split into hg and mm"
			chrs1=`samtools view -H $Name.$Species.rigid.reheader.unique.st.bam | grep hg19 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
			chrs2=`samtools view -H $Name.$Species.rigid.reheader.unique.st.bam | grep mm10 | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<8)print}'`
			samtools view -@ $cores -b $Name.$Species.rigid.reheader.unique.st.bam -o $Name.$Species.temp1.bam `echo ${chrs1[@]}`
			samtools view -@ $cores -b $Name.$Species.rigid.reheader.unique.st.bam -o $Name.$Species.temp2.bam `echo ${chrs2[@]}`
			samtools view -@ $cores -h $Name.$Species.temp1.bam | sed 's/hg19_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ $cores -b -o $Name.hg19.rigid.reheader.unique.st.bam
			samtools view -@ $cores -h $Name.$Species.temp2.bam | sed 's/mm10_/chr/g' | sed 's/chrMT/chrM/g'| samtools view -@ $cores -b -o $Name.mm10.rigid.reheader.unique.st.bam

			samtools index -@ $cores $Name.hg19.rigid.reheader.unique.st.bam &
			samtools index -@ $cores $Name.mm10.rigid.reheader.unique.st.bam &
			wait
			rm $Name.$Species.temp1.bam $Name.$Species.temp2.bam
		    else
			echo "Single species is aligned"
		    fi		
		fi

		# assign feature to reads
		cd  $dir/fastqs/
		if [ ${Type[$index]} == "RNA" ]; then
		    for Species2 in ${Genome2[@]}; do
			if [ -f $dir/fastqs/$Name.$Species2.exon.featureCounts.bam ] || [ -f $dir/fastqs/$Name.$Species2.wdup.bam.bai ]; then
			    echo "Skip exon feasure count"
			else
		    	    # excliude multimapping, uniquely mapped reads only, -Q 30, for real sample, might consider include multi-mapping
			    echo "Feature counting on exons"
			    # count exon
                            if [ ${keepMultiMapping[$index]} == "T" ]; then
				featureCounts -T $cores -Q 0 -M -a $myPATH/gtf/$Species2.$refgene.gtf -t exon -g $genename \
                                    -o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.rigid.reheader.unique.st.bam >>$dir/Run.log
                            else
				featureCounts -T $cores -Q 30 -a $myPATH/gtf/$Species2.$refgene.gtf -t exon -g $genename \
                                    -o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.rigid.reheader.unique.st.bam >>$dir/Run.log
                            fi
			    # Extract reads that assigned to genes
			    mv $Name.$Species2.rigid.reheader.unique.st.bam.featureCounts.bam $Name.$Species2.exon.featureCounts.bam
			fi
			if [ -f $dir/fastqs/$Name.$Species2.wdup.bam.bai ]; then
			    echo "Skip intron and exon feasure count"
			else
			    # count both intron and exon
			    echo "Count feature on both intron and exon"
			    if [ $keepIntron == "T" ]; then
				if [ ${keepMultiMapping[$index]} == "T" ]; then
				    featureCounts -T $cores -Q 0 -M -a $myPATH/gtf/$Species2.$refgene.gtf -t gene -g $genename \
						  -o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.exon.featureCounts.bam >>$dir/Run.log
				# then for reads mapped to multiple genes,  only keep if all its alignments are within a single gene
				else
				    featureCounts -T $cores -Q 30 -a $myPATH/gtf/$Species2.$refgene.gtf -t gene -g $genename \
					-o $Name.$Species2.feature.count.txt -R BAM $Name.$Species2.exon.featureCounts.bam >>$dir/Run.log
				fi
				samtools sort -@ $cores -m 2G -o $Name.$Species2.wdup.bam  $Name.$Species2.exon.featureCounts.bam.featureCounts.bam
				rm $Name.$Species2.exon.featureCounts.bam.featureCounts.bam
			    else
				samtools sort -@ $cores -m 2G -o $Name.$Species2.wdup.bam $Name.$Species2.exon.featureCounts.bam
			    fi
			    samtools index -@ $cores $Name.$Species2.wdup.bam
			fi
		    done
		fi
		
		# group UMIs
		if [ ${Type[$index]} == "RNA" ]; then
		    for Species2 in ${Genome2[@]}; do
			if [ -f $dir/fastqs/$Name.$Species2.grouped.bam ] || [ -f $dir/fastqs/$Name.$Species2.groups.tsv.gz ] || [ -f $dir/fastqs/$Name.$Species2.groups.tsv ]; then
			    echo "Found $Name.$Species.groups.bam, skip grouping UMIs"
			else
			    echo "Group reads to unique UMIs"
			    if [ $mode == "regular" ]; then
				## this is previously used. Seems to get more slant and fewer UMIs, but get accurate lib size estimation
				umi_tools group --extract-umi-method=read_id \
                                          --per-gene --gene-tag=XT --per-cell \
                                          -I $dir/fastqs/$Name.$Species2.wdup.bam \
                                          --output-bam -S $dir/fastqs/$Name.$Species2.grouped.bam \
                                          --group-out=$dir/fastqs/$Name.$Species2.groups.tsv --skip-tags-regex=Unassigned >>$dir/Run.log
			    else
				## own UMI dedup by matching bc-umi-align position
                                samtools view -@ $cores $Name.$Species2.wdup.bam | grep XT:Z: | \
                                    sed 's/Unassigned_Ambiguity/discard/g' | \
                                    sed 's/Unassigned_MappingQuality/discard/g' | \
                                    awk 'gsub(/[_]/,"\t", $1)' | \
                                    awk -v OFS='\t' '{if($NF ~/discard/){$NF=$(NF-1)} print $5, $6, $2, $3, $NF}' | \
                                    sed 's/XT:Z://g' > $Name.$Species2.wdup.bed
                                # remove dup reads
                                python3 $myPATH/rm_dup_barcode_UMI_v3.py \
                                        -i $Name.$Species2.wdup.bed -o $Name.$Species2.groups.tsv --m 1
#                                pigz --fast -p $cores $dir/fastqs/$Name.$Species2.groups.tsv
                                rm $Name.$Species2.wdup.bed
			    fi
			fi
			# convert groupped UMI to bed file
			if [ -f $dir/fastqs/$Name.$Species2.bed.gz ]; then
			    echo "Skip removing single-read UMIs and GGGG"
			else
			    if [ $removeSingelReadUMI == "T" ]; then
				echo "Filter UMI that has only 1 read and convert to bed file"
				if [ $mode == "regular" ]; then
				    # in regular mode, duplicate UMI reads are kept in the tsv, so I remove those reads first (awk 'NR==1{ id="N"} {if(id != $11 ) {id = $11; print $2, $6, $10}}')
				    less $dir/fastqs/$Name.$Species2.groups.tsv | \
					awk 'gsub(/[_]/,"\t", $1)' | \
                                       	awk 'FNR>1 {if($10 >1){if($9 != "GGGGGGGGGG"){print}}}' | \
                                        awk 'NR==1{ id="N"} {if(id != $11 ) {id = $11; print $2, $6, $10}}' | \
                                        awk -v OFS="\t" 'NR==1{ t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} \
                                            else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}}' | \
                                        pigz --fast -p $cores > $Name.$Species2.bed.gz
				else
				    less $dir/fastqs/$Name.$Species2.groups.tsv | \
					awk -v OFS="\t" '{if($3 >1){print}}' | \
					sort  --parallel=$cores -S 24G -k1,1 -k2,2 | \
					awk -v OFS="\t" 'NR==1 { t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}} END {print t1, t2, umisum, readsum}' | \
                                        pigz --fast -p $cores > $Name.$Species2.bed.gz
				fi
			    else
				echo "Convert to bed file"
				if [ $mode == "regular" ]; then
				    ## count reads per gene
				    less $dir/fastqs/$Name.$Species2.groups.tsv | \
					awk 'gsub(/[_]/,"\t", $1)' | \
					awk 'FNR>1 {if($9 != "GGGGGGGGGG"){print}}' | \
					awk 'NR==1 { id="N"} {if(id != $11 ) {id = $11; print $2, $6, $10}}' | \
					awk -v OFS="\t" 'NR==1 { t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}}' | \
					pigz --fast -p $cores > $Name.$Species2.bed.gz
					
					## 3rd column is umi, 4th column is read
					## note: difference between these two versions of groups.tsv
					## 1) umitools output keep all the barcode-UMI and don't collapse them
					## 2) my script already collpsed them at alignment position level
				else
				    less $dir/fastqs/$Name.$Species2.groups.tsv | \
					sort --parallel=$cores -S 24G -k1,1 -k2,2 | \
					awk -v OFS="\t" 'NR==1 { t1=$1;t2=$2;readsum=0; umisum=0} {if(t1==$1 && t2==$2) {readsum+=$3; umisum+=1} else {print t1, t2, umisum, readsum; t1=$1;t2=$2;umisum=1;readsum=$3}} END {print t1, t2, umisum, readsum}' | \
					pigz --fast -p $cores > $Name.$Species2.bed.gz
				fi
			    fi
			fi

			# count reads for RNA
			if [ -f $Name.$Species2.filtered.counts.csv ] || [ -f $Name.$Species.filtered.counts.csv.gz ]; then
                            echo "found $Name.$Species.filtered.counts.csv"
			else
                            echo "count unfiltered reads"
			    zcat $Name.$Species2.bed.gz | \
				awk -v OFS='\t' '{a[$1] += $4} END{for (i in a) print a[i], i}' | \
				awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }'> $Name.$Species2.wdup.RG.freq.bed
                            /usr/local/bin/Rscript $myPATH/sum_reads_v2.R $dir/fastqs/ $Name.$Species2.wdup.RG.freq.bed --save
                            mv $Name.$Species2.wdup.RG.freq.bed.csv $Name.$Species2.unfiltered.counts.csv

                            echo "count filtered reads"
                            zcat $Name.$Species2.bed.gz | \
				awk -v OFS='\t' '{a[$1] += $3} END{for (i in a) print a[i], i}' | \
				awk -v OFS='\t' '{if($1 >= '${ReadsPerBarcode[$index]}') print }' > $Name.$Species2.rmdup.RG.freq.bed
                            /usr/local/bin/Rscript $myPATH/sum_reads_v2.R $dir/fastqs/ $Name.$Species2.rmdup.RG.freq.bed --save
                            mv $Name.$Species2.rmdup.RG.freq.bed.csv $Name.$Species2.filtered.counts.csv

                            rm $Name.$Species2.wdup.RG.freq.bed $Name.$Species2.rmdup.RG.freq.bed
			fi		
			
			# remove barcode combination that has less then N reads
			if [ -f $Name.$Species2.cutoff.bed.gz ]; then
			    echo "Skip removing low counts barcode combination"
			else
			    echo "Remove low counts barcode combination"
			    sed -e 's/,/\t/g' $Name.$Species2.filtered.counts.csv | \
				awk -v OFS=',' 'NR>=2 {if($5 >= '${ReadsPerBarcode[$index]}') print $1,$2,$3,$4} '  > $Name.$Species2.barcodes.txt
				grep -wFf $Name.$Species2.barcodes.txt <(zcat $Name.$Species2.bed.gz) | pigz --fast -p $cores > $Name.$Species2.cutoff.bed.gz
			fi
			# deduplicats for non-UMI reads
		    done
		fi
		
		# Gene body coverage and reads distribution
		if [ ${Type[$index]} == "RNA" ]; then
		    # split bam to hg and mm
		    cd $dir/fastqs/
		    for Species2 in ${Genome2[@]}; do
			if [ -f $dir/fastqs/$Name.$Species2.read_distribution.txt ]; then
			    echo "Skip calculate read disbution"
			else
			    echo "Calculate read distribution"
			    if [ $Runtype = QC ]; then
				read_distribution.py -i $dir/fastqs/$Name.$Species2.wdup.bam -r $genomeBed/$Species2.UCSC_RefSeq.bed > $Name.$Species2.read_distribution.txt 2>>$dir/Run.log
			    else
				# only use 1% of reads
				samtools view -s 0.01 -o $dir/fastqs/temp.bam $dir/fastqs/$Name.$Species2.wdup.bam
				read_distribution.py -i $dir/fastqs/temp.bam -r $genomeBed/$Species2.UCSC_RefSeq.bed > $Name.$Species2.read_distribution.txt 2>>$dir/Run.log
				rm $dir/fastqs/temp.bam
			    fi
			    # plot reads disbution
			    tail -n +5 $Name.$Species2.read_distribution.txt | head -n -1 > temp1.txt
			    head -n 3  $Name.$Species2.read_distribution.txt | grep Assigned | sed 's/Total Assigned Tags//g' | sed 's/ //g' > temp2.txt
			    /usr/local/bin/Rscript $myPATH/Read_distribution.R $dir/fastqs/ $Name.$Species2 --save
			    rm temp1.txt temp2.txt
			fi
		    done
		fi 

		if [ ${Type[$index]} == "ATAC" ]; then
		    # get final quality stats
		    echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Raw" > $Name.$Species.stats.log
		    samtools idxstats $Name.$Species.st.bam >> $Name.$Species.stats.log
		    echo -e "Chromosome\tLength\tProperPairs\tBadPairs:Filtered" >> $Name.$Species.stats.log
		    samtools idxstats $Name.$Species.rmdup.bam >> $Name.$Species.stats.log
		    
		    # get insert-sizes
		    if [ -f $Name.$Species.rmdup.hist_data.pdf ]; then
			echo "Skip checking insert size"
		    else
			echo '' > $Name.$Species.rmdup.hist_data.log
		    fi
		    if [ ! -f $Name.$Species.RefSeqTSS ]; then		
			# make TSS pileup fig
			echo "Create TSS pileup"; set +e
			$toolPATH/pyMakeVplot.py -a $Name.$Species.rmdup.bam -b $tssFilesPATH/$Species.TSS.bed -e 2000 -p ends -v -u -o $Name.$Species.RefSeqTSS
		    fi
		fi		
       		if [ -d "$dir/fastqs/tmp" ]; then rm -r $dir/fastqs/tmp; fi
	    fi
	done
	
	for Species2 in ${Genome2[@]}; do
	    if [ ${Type[$index]} == "RNA" ] || [ ${Type[$index]} == "crop" ] || [ ${Type[$index]} == "cite" ] || [ ${Type[$index]} == "cellhash" ]; then
	  	# plot UMI/cell or gene/cell
		echo "Convert bed to UMI count matrix, plot UMI_gene_perCell, generate h5"
		/usr/local/bin/Rscript $myPATH/UMI_gene_perCell_plot_v3.R $dir/fastqs/ $Name.$Species2 --save
	    fi
	done


	# estimate lib size
	if [ -f $Name.counts.csv ]; then
	    echo "Found $Name.counts.csv, skip calculate lib size"
	else
	    echo "Estimate lib size"
            if [ ${Genomes[$index]} == "both" ]; then
		if [ -f $Name.hg19.unfiltered.counts.csv ] && [ ! -f $Name.hg19.filtered.counts.csv ]; then
		    cp $Name.hg19.unfiltered.counts.csv $Name.hg19.filtered.counts.csv
		    echo "Error: Could locate $Name.hg19.filtered.counts.csv"
		fi
		if [ -f $Name.mm10.unfiltered.counts.csv ] && [ ! -f $Name.mm10.filtered.counts.csv ]; then
		    cp $Name.mm10.unfiltered.counts.csv $Name.mm10.filtered.counts.csv
		    echo "Error: Could locate $Name.mm10.filtered.counts.csv"
		fi
		if [ -f $Name.hg19.unfiltered.counts.csv ] && [ -f $Name.mm10.unfiltered.counts.csv ]; then
                    echo "Calcuating library size for $Name"
		    /usr/local/bin/Rscript $myPATH/lib_size_sc_V5_species_mixing.R ./ $Name ${ReadsPerBarcode[$index]} ${Type[$index]} --save
		fi
	    else
		if [ -f $Name.${Genomes[$index]}.unfiltered.counts.csv ] && [ ! -f $Name.${Genomes[$index]}.filtered.counts.csv ]; then
		    cp $Name.${Genomes[$index]}.unfiltered.counts.csv $Name.${Genomes[$index]}.filtered.counts.csv
		else
		    /usr/local/bin/Rscript $myPATH/lib_size_sc_V5_single_species.R ./ $Name ${ReadsPerBarcode[$index]} ${Genomes[$index]} ${Type[$index]} --save
		fi
	    fi
	fi
	if [ ! -d $dir/$Name/ ]; then 
	    mkdir $dir/$Name/ && mv $dir/fastqs/$Name.* $dir/$Name
	fi
    fi
    index=$((index + 1))
done

cd $dir
# get stats for ATAC
echo "Name" > Names.atac.xls
index=0 && count=0
if [ -d $dir/temp/ ]; then rm -r $dir/temp/; fi
mkdir $dir/temp/
for Name in ${Project[@]}; do
    if [ ${Type[$index]} == "ATAC" ]; then
	if [ ${Genomes[$index]} == "both" ]; then
	    Genome2=(hg19 mm10)
	    for Species in ${Genome2[@]}; do
		echo $Name.$Species >> Names.atac.xls
		mkdir $dir/temp/$Name.$Species/
		cp $dir/$Name/$Name.$Species.dups.log $dir/temp/$Name.$Species/
		cp $dir/$Name/$Name.$Species.align.log $dir/temp/$Name.$Species/
		cp $dir/$Name/$Name.$Species.stats.log $dir/temp/$Name.$Species/
		cp $dir/$Name/$Name.$Species.rmdup.hist_data.log $dir/temp/$Name.$Species/
		cp $dir/$Name/$Name.$Species.RefSeqTSS $dir/temp/$Name.$Species/
		let "count=count+1"
            done
        else
            Species=${Genomes[$index]}
	    echo $Name.$Species >> Names.atac.xls
	    mkdir $dir/temp/$Name.$Species/
	    cp $dir/$Name/$Name.$Species.dups.log $dir/temp/$Name.$Species/
	    cp $dir/$Name/$Name.$Species.align.log $dir/temp/$Name.$Species/
	    cp $dir/$Name/$Name.$Species.stats.log $dir/temp/$Name.$Species/
	    cp $dir/$Name/$Name.$Species.rmdup.hist_data.log $dir/temp/$Name.$Species/
	    cp $dir/$Name/$Name.$Species.RefSeqTSS $dir/temp/$Name.$Species/
	fi
    fi
    index=$((index + 1))
done
cp -r $dir/temp/$Name.$Species/ $dir/temp/$Name.$Species.temp/

if [ $count -eq 0 ]; then
    rm Names.atac.xls
else
    cd $dir/temp/
    cp $dir/Names.atac.xls $dir/temp/
    /usr/bin/python $myPATH/pySinglesGetAlnStats_sai.py -i ./Names.atac.xls -d ./
    mv $dir/temp/Names.merged.xls $dir/Names.atac.merged.xls
    mv Names.merged.xls.iSize.txt.mat Names.atac.iSize.txt.mat
fi

# gather useful files
cd $dir
mkdir $dir/Useful
cp $dir/*/*.png $dir/Useful
cp $dir/*/*.pdf $dir/Useful

# generate bigwig, for TAPS-seq and RNA-seq
let i=0 
unset Subproject
for (( Index=0; Index<${#Type[@]}; Index++ ))
do
    if [ ${Type[$Index]} == "RNA" ] || [ ${Type[$Index]} == "ATAC" ]; then
	Subproject[$i]=${Project[$Index]}
	let "i+=1"
    fi
done

count=`ls -1 $dir/Useful/*bw 2>/dev/null | wc -l`
Genome=(hg19 mm10 hg38) # define genome

if [ ! -z $Subproject ]; then
    if [ $count != 0 ]; then
	echo "Found bigwig files, skip converting bam to bigwig"
    else
	echo "Genarate bigwigs"
	for Species in ${Genome[@]}; do
	    parallel --jobs 6 --delay 1 'igvtools count -w 25 -f mean -e 250 '$dir'/{}/{}.'$Species'.rmdup.bam '$dir'/Useful/{}.'$Species'.rmdup.wig '$genomeBed'/'$Species'.chrom.sizes >>'$dir'/Run.log' ::: `echo ${Subproject[@]}`
	    cd $dir/Useful
	    parallel --jobs 6 --delay 1 'wigToBigWig {} '$genomeBed'/'$Species'genome.bed {.}.bw 2>>'$dir'/Run.log' ::: *wig
	    rm *wig
	done
    fi
fi

# clean up
echo "Cleaning up folder"
cp $dir/*/*.csv $dir/Useful/
rm $dir/Useful/*filtered.counts.csv 
rm -r $dir/temp/

pigz --fast -p $cores $dir/*/*.csv
pigz --fast -p $cores $dir/*/*.groups.tsv


rm $dir/*/*wdup.all.bam* $dir/*/*namesort.bam $dir/igv.log $dir/*/*exon.featureCounts.bam 
rm $dir/*/*grouped.bam $dir/*/*rigid.reheader.unique.st.bam*

if [ $cleanup == "T" ]; then
    rm $dir/fastqs/discard.R1.fastq.gz touch $dir/fastqs/discard.R1.fastq.gz
    rm $dir/fastqs/discard.R2.fastq.gz touch $dir/fastqs/discard.R2.fastq.gz
    $dir/*/*.groups.tsv.gz
fi

echo "The pipeline is completed!! Author: Sai Ma <sai@broadinstitute.org>"

exit
