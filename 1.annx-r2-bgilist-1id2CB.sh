###get file path
OUTDIR=/path/outdir
mkdir -p $OUTDIR
INDIR=/path/datadir
###get all chip files
files=`ls $INDIR/*.fastq.gz`
echo $files
###loop for dealing with each chip
for file in $files
    do
    filename=`basename $file _R1.fastq.gz`
    samplename=$filename
    DIR=$OUTDIR/$samplename
    mkdir -p $DIR
    cd $DIR

    mkdir -p CellBarcodes
    outfile=$DIR/$samplename\_mappingratio.txt
    touch $outfile
    mv $INDIR/$filename\_R1.fastq.gz $INDIR/$filename\_R1.fastq
    mv $INDIR/$filename\_R2.fastq.gz $INDIR/$filename\_R2.fastq
    outdir=$DIR 
    fq1=$INDIR/$filename\_R1.fastq
    fq2=$INDIR/$filename\_R2.fastq 
    idxfile=$INDIR/anx5-barcodes.txt

    ###remove adapters
    fastp -n 10 -q 35 -w 40 --detect_adapter_for_pe -i $fq1 -I $fq2 -o $samplename.1.clean.fq.gz -O $samplename.2.clean.fq.gz -j $samplename.json -h $samplename.html -R $samplename
    rm *.clean.*

    ###get fastq containing reads with ANXV barcode
    fastq-multx -B $idxfile -m 4 -b $fq1 $fq2 -o %.R1.fastq -o %.R2.fastq
    a=`zcat $fq1 |wc -l`
    echo "RawReads $((a/4))" > $outfile
    d=0
    e=0
    echo "CellBarcodes	Count	Sample" > CBlist_$samplename.txt
    pre="ANX5-labeled"
    nfq1=$pre.R1.fastq
    nfq2=$pre.R2.fastq
    b=`wc -l < $pre.R1.fastq`
    echo "$pre.labled $((b/4))" >> $outfile
    ###get ANXV labeled cell barcodes
    python3.7 $INDIR/whitelist.py --stdin $nfq1 --read2-in=$nfq2 --error-correct-threshold=0 --method=umis --extract-method=regex --bc-pattern2="(?P<discard_1>CTACGATCCGACTTTCTGCG){s<=2}(?P<cell_1>.{10})(?P<discard_2>CCTTCC){s<=1}(?P<cell_2>.{10})(?P<discard_3>CGATG){s<=1}(?P<umi_1>.{10})TT.*" --plot-prefix=$pre\_expect_whitelist --log2stderr --expect-count=0 --knee-method=None > $pre\_whitelist.txt
    ###filter with BGI cell barcode whitelist
    join $pre\_whitelist.txt $INDIR/BGI_droplet_scRNA_readStructureV2_cDNA_T1-2.sort.txt |awk '{print $1"\t\t"$2"\t"}' > comm.$pre\_whitelist.txt
    ###get ANXV labeled cell fastq
    umi_tools extract -I $nfq1 --read2-in=$nfq2 --extract-method=regex --bc-pattern2="(?P<discard_1>CTACGATCCGACTTTCTGCG){s<=2}(?P<cell_1>.{10})(?P<discard_2>CCTTCC){s<=1}(?P<cell_2>.{10})(?P<discard_3>CGATG){s<=1}(?P<umi_1>.{10})TT.*" --stdout=truecell.$pre.1.fastq --read2-out=truecell.$pre.2.fastq --filtered-out=Nocell.$pre.1.fastq --filtered-out2=Nocell.$pre.2.fastq --whitelist=comm.$pre\_whitelist.txt
    c=`wc -l < truecell.$pre.1.fastq`
    echo "$pre.CB $((c/4))" >> $outfile
    ((d=d+c))
    ((e=e+b))
    n=`wc -l < $pre\_whitelist.txt`
    echo "raw.$pre.CBtypes $n" >> $outfile
    n=`wc -l < comm.$pre\_whitelist.txt`
    echo "$pre.CBtypes $n" >> $outfile
    ###final ANXV CB list
    awk -v a=$pre 'BEGIN { OFS = "\t"} {print $0,a} ' comm.$pre\_whitelist.txt >> CBlist_$samplename.txt

    echo "withCBTotal   $((d/4))" >> $outfile
    echo "withSLTotal   $((e/4))" >> $outfile
    ###clean folder
    mv *whitelist.txt CellBarcodes/
    mv *.png CellBarcodes/
    mv *.txt CellBarcodes/
    rm *.fastq
done
