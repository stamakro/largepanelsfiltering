
ml load bedtools;

# get a list of all mutations
cat ../../vcf/Qiaseq-filtered/*filtered.vcf | grep -v '##' | grep -v '#CH' | cut -f 1,2 | sort -u > ../../all-qiaseq-calls.txt;


mkdir -p ../../mutation-contexts;

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

  bedtools getfasta -fi '../../../human-genome/fasta/chr'$chrom'.fa' -bed '../../mutation-contexts/locations_chr'$chromName'.bed' -tab > '../../mutation-contexts/sequences_chr'$chromName'.txt';
done;

# get a list of all mutations
cat ../../vcf/Qiaseq-new-patients/filtered/*filtered.vcf | grep -v '##' | grep -v '#CH' | cut -f 1,2 | sort -u > ../../all-qiaseq-calls-new.txt;


mkdir -p ../../mutation-contexts-new;

python prepareGenomicContext.py '-new';


for chrom in '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y' 'M';
do
  echo $chrom
  if [ $chrom == 'M' ];
  then
    chromName='MT';
  else
    chromName=$chrom;
  fi;

  bedtools getfasta -fi '../../../human-genome/fasta/chr'$chrom'.fa' -bed '../../mutation-contexts-new/locations_chr'$chromName'.bed' -tab > '../../mutation-contexts-new/sequences_chr'$chromName'.txt';
done;
