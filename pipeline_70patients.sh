#! /bin/sh

VCFDIR='../../vcf/Qiaseq/';
VCFFILTEREDDIR='../../vcf/Qiaseq-filtered/';
DIRc='../../vcf/Qiaseq-common-regions/';
DIRu='../../vcf/Qiaseq-unique-regions/';
QIASEQTARGETSFILE='../../gene-panels/DHS-001Z.primers-150bp_GRCh38.chrNames.bed';
COMMONREGIONSFILE='../../gene-panels/common_regions_sorted_merged_38.chrNames.bed';
ALLCALLSFILE='../../all-qiaseq-calls.txt';
MUTATIONCONTEXTDIR='../../mutation-contexts/';
HUMANGENOMEPATH='../../../human-genome/fasta/';
SAVELOCATIONFEATURESCOMMON='../../featureVectorsCommon.pkl';
SAVELOCATIONFEATURESUNIQUE='../../featureVectorsUnique.pkl';
ONCOMINEDIR='../../vcf/Oncomine/';
LOPOCALLSFILE='../../common_regions_loo_calls.csv';
LOPORESULTSFILE='./result_reducedFeats_cosmic.pkl';
BASELINERESULTSFILE='./result_baseline.pkl';
UNIQUECALLSFILE='../allqiaseqCalls_filtered_reducedFeats_cosmic.csv';
MODELSAVEFILE='trained_model_reduced_cosmic.pkl';
VEPFILESPATH='/tudelft.net/staff-umbrella/liquidbiopsy/software/vep';
VEPCONTAINERPATH='~/vep/ensembl-vep_latest.sif';
ONCOMINECALLS='../../vcf/Oncomine-108-samples.csv'
ml load bedtools/2.27;


mkdir -p $VCFFILTEREDDIR;

# basic filtering
for vcf in $(ls $VCFDIR*.vcf.gz | grep -v run1 | grep -v ori);
do

  if grep -q 'run2' <<< $vcf;
  then
    outfile=$(awk -F "/" '{print $5}' <<< $vcf | awk -F "-" '{print $5}');
  else
    outfile=$(awk -F "/" '{print $5}' <<< $vcf | awk -F "-" '{print $4}');
  fi

  zcat $vcf | grep -v "INDEL" > $VCFFILTEREDDIR$outfile'_snps.vcf';

  python qiaseqBasicFiltering.py $VCFFILTEREDDIR$outfile'_snps.vcf';

done;


mkdir -p $DIRc;
mkdir -p $DIRu;


# split into regions common with oncomine and qiaseq only
for vcf in $VCFFILTEREDDIR*filtered.vcf;
do

  outfile=$( awk -F "/" '{print $5}' <<< $vcf | awk -F "_" '{print $1".vcf"}');
  grep '##' $vcf > tmp.vcf
  grep '#CHROM' $vcf >> tmp.vcf

  # remove calls outside the Qiaseq target regions
  bedtools intersect -a $vcf -b $QIASEQTARGETSFILE | sort -u >> tmp.vcf;

  # split to common and unique regions
  bedtools intersect -a tmp.vcf -b $COMMONREGIONSFILE > $DIRc$outfile;
  bedtools intersect -a tmp.vcf -b $COMMONREGIONSFILE -v > $DIRu$outfile;
  rm tmp.vcf;
done;

# get a list of all mutations
cat $VCFFILTEREDDIR*filtered.vcf | grep -v '##' | grep -v '#CH' | cut -f 1,2 | sort -u > $ALLCALLSFILE;


mkdir -p $MUTATIONCONTEXTDIR;

# this script will be affected by changes in variables ALLCALLSFILE and MUTATIONCONTEXTDIR
python prepareGenomicContext.py;


for chrom in '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y' 'M';
do
  echo $chrom
  if [ $chrom == 'M' ];
  then
    chromName='MT';
  else
    chromName=$chrom;
  fi;

  bedtools getfasta -fi $HUMANGENOMEPATH'chr'$chrom'.fa' -bed $MUTATIONCONTEXTDIR'locations_chr'$chromName'.bed' -tab > $MUTATIONCONTEXTDIR'sequences_chr'$chromName'.txt';
done;

# this script will be affected by changes in variables ALLCALLSFILE and MUTATIONCONTEXTDIR
python extractGenomicContext.py;

# finalize feature extraction
python extractFeatures.py $DIRc $DIRu $SAVELOCATIONFEATURESCOMMON $SAVELOCATIONFEATURESUNIQUE '';

mkdir -p './.tmp/';
# run leave one patient out to evaluate model (inner cross validation loop)
python inner_loop_lopo.py 0 $ONCOMINEDIR $DIRc $SAVELOCATIONFEATURESCOMMON $LOPOCALLSFILE $LOPORESULTSFILE;

cp $LOPOCALLSFILE '.tmp/commonCalled.csv';
# replace this with your own vep installation path/container
# this command contains some absolute paths, see VEP documentation for how to set VEP up
singularity run --bind /tudelft.net/staff-umbrella/liquidbiopsy/software/vep:/mnt,/tudelft.net/staff-umbrella/liquidbiopsy/SONIA/scripts/qiaseq-filtering/.tmp/:/home/nfs/stavrosmakrodi ~/vep/ensembl-vep_latest.sif vep -i ~/calls.vcf -o ~/calls_annotated.txt --most_severe --cache --dir_cache /mnt --cache_version 107 --offline --force_overwrite

tail -n+42 ./.tmp/calls_annotated.txt > ./.tmp/calls_annotated.tsv

python postprocess_general.py '../../finalCommonRegionCallsFiltered.csv';


# run baseline model with vaf, quality etc.
python baseline.py $ONCOMINEDIR $ONCOMINECALLS $DIRc $SAVELOCATIONFEATURESCOMMON $BASELINERESULTSFILE;

# run outer loop to filter variants in qiaseq-specific regions
# first this scripts uses all variants from the common regions
# to optimize model parameters
python outerLoop_reducedFeats_cosmic.py $ONCOMINEDIR $DIRc $SAVELOCATIONFEATURESCOMMON $MODELSAVEFILE $SAVELOCATIONFEATURESUNIQUE $UNIQUECALLSFILE $DIRu

# replace this with your own vep installation path/container
# this command contains some absolute paths, see VEP documentation for how to set VEP up
singularity run --bind /tudelft.net/staff-umbrella/liquidbiopsy/software/vep:/mnt,/tudelft.net/staff-umbrella/liquidbiopsy/SONIA/scripts/qiaseq-filtering/.tmp/:/home/nfs/stavrosmakrodi ~/vep/ensembl-vep_latest.sif vep -i ~/calls.vcf -o ~/calls_annotated.txt --most_severe --cache --dir_cache /mnt --cache_version 107 --offline --force_overwrite

tail -n+42 ./.tmp/calls_annotated.txt > ./.tmp/calls_annotated.tsv;

mv .tmp/uniqueCalled.csv .tmp/commonCalled.csv;
python postprocess_general.py '../../finalUniqueRegionCallsFiltered.csv';
