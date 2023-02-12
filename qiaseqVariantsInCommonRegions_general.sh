ml load bedtools/2.27.0;

DIRIN=$1;
DIR=$DIRIN'/target';
mkdir -p $DIR;

NAME=$2;

DIROUT=$3;
mkdir -p $DIROUT;

echo 'Keeping variants within target regions...';
for vcf in $(ls $DIRIN/*.vcf | grep -v ori | grep -v target);
do
  echo $vcf;
  #$5 new 8
  outfile=$(awk -F "/" '{print $5}' <<< $vcf);
  echo $outfile;
  grep '##' $vcf > tmp.vcf
  grep '#CHROM' $vcf >> tmp.vcf
  bedtools intersect -a $vcf -b ../../gene-panels/DHS-001Z.primers-150bp_GRCh38.chrNames.bed | sort -u >> tmp.vcf;

  mv tmp.vcf $DIR'/'$outfile;

done;

DIR=$DIRIN'/filtered';
mkdir -p $DIR;

echo 'basic filtering...';
# # the third field with '-' as delimiter is kept as id by default
for vcf in $(ls $DIRIN'/target/'*'.vcf' | grep -v ori);
do
  echo $vcf;
  # $6 $5
  outfile=$(awk -F "/" '{print $6}' <<< $vcf | awk -F "-" '{print $5}' | cut -d '.' -f 1);
  echo $outfile;
  # only keep SNVs in COSMIC
  grep -v "INDEL" $vcf | grep 'COSV' > $DIR'/'$outfile'_snps.vcf';

  python qiaseqBasicFiltering.py $DIR'/'$outfile'_snps.vcf';
  cp $DIR'/'$outfile'_snps_filtered.vcf' $DIROUT'/'$outfile'.vcf';

done;

echo 'preparing data for extracting genomic context features...';
# get a list of all mutations
cat $DIR/*filtered.vcf | grep -v '##' | grep -v '#CH' | cut -f 1,2 | sort -u > ../../all-qiaseq-calls-$NAME.txt;

mkdir -p ../../mutation-contexts-$NAME;

python prepareGenomicContext.py -$NAME;


for chrom in '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y' 'M';
do
  echo $chrom
  if [ $chrom == 'M' ];
  then
    chromName='MT';
  else
    chromName=$chrom;
  fi;

  bedtools getfasta -fi '../../../human-genome/fasta/chr'$chrom'.fa' -bed '../../mutation-contexts-'$NAME'/locations_chr'$chromName'.bed' -tab > '../../mutation-contexts-'$NAME'/sequences_chr'$chromName'.txt';
done;


echo 'extracting genomic context features...';
python extractGenomicContext.py -$NAME;

echo 'generating and saving final feature vectors...';
python extractFeatures.py $DIROUT'/' '' $NAME'.pkl' '' '-'$NAME;


mkdir -p ./.tmp

echo 'classifying variants with model...'
# if third argument is given and is not 0, model filter is turned off
# example:
# python outerLoop_general.py $DIROUT'/' $NAME'.pkl' 1 
python outerLoop_general.py $DIROUT'/' $NAME'.pkl'

echo 'annotating variants with VEP...'
# replace this with your own vep installation path/container
singularity run --bind /tudelft.net/staff-umbrella/liquidbiopsy/software/vep:/mnt,/tudelft.net/staff-umbrella/liquidbiopsy/SONIA/scripts/qiaseq-filtering/.tmp/:/home/nfs/stavrosmakrodi ~/vep/ensembl-vep_latest.sif vep -i ~/calls.vcf -o ~/calls_annotated.txt --most_severe --cache --dir_cache /mnt --cache_version 107 --offline --force_overwrite

tail -n+42 ./.tmp/calls_annotated.txt > ./.tmp/calls_annotated.tsv

echo 'post-processing variants...'
python postprocess_general.py

#echo 'cleaning up...'
#rm -rf ./.tmp
