#! /bin/sh
#SBATCH --partition=general --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=50000
#SBATCH --time=00:45:00
#SBATCH --job-name=sonia
#SBATCH --mail-user=stavrosmakrodi
#SBATCH --mail-type=ALL


# load bedtools (v2.28 gives a weird segmentation fault for some files)
ml use /opt/insy/modules/DBL/modulefiles;
ml load bedtools/2.27.0;
ml load samtools;


ROOT='/tudelft.net/staff-umbrella/liquidbiopsy/SONIA/';
BEDPATH=$ROOT'gene-panels/';
BAMPATH=$ROOT'bam/';


# find common regions
bedtools intersect -a $	BEDPATH'Oncomine_Breast_cfDNA_v2.08212017.Designed.bed' -b $BEDPATH'DHS-001Z.primers-150bp.bed' > $BEDPATH'common_regions.bed';
#
bedtools sort -i $BEDPATH'common_regions.bed' > $BEDPATH'common_regions_sorted.bed'
bedtools merge -i $BEDPATH'common_regions_sorted.bed' | sort -V > $BEDPATH'common_regions_sorted_merged.bed'

#! !! the qiaseq reads are mapped to grch38! for now manually liftover common_regions.bed and save it as common_regions_sorted_merged_38.bed in the same dir  !!!
# !!! let's how many more times we have to do this
#
# remove the 'chr' from the qiaseq bed file because it is not in the bam file
cat $BEDPATH'DHS-001Z.primers-150bp_GRCh38.bed' | tr -d "chr" > $BEDPATH'DHS-001Z.primers-150bp_GRCh38.chrNames.bed'
cat $BEDPATH'common_regions_sorted_merged_38.bed' | tr -d "chr" > $BEDPATH'common_regions_sorted_merged_38.chrNames.bed'

# measure bp overlap in each location
bedtools intersect -a $BEDPATH'Oncomine_Breast_cfDNA_v2.08212017.Designed.bed' -b $BEDPATH'DHS-001Z.primers-150bp.bed' -wao > $BEDPATH'overlap.tsv';


COVERAGEPATH=$BEDPATH'coverage/'
mkdir -p $COVERAGEPATH;

# check that all oncomine samples have reads on all of the desired regions of the oncomine assay
mkdir -p $COVERAGEPATH'oncomine-on-oncomine';
#
for f in $BAMPATH'Oncomine/'*.bam;
do
	echo $f;
	f2=$(awk -F "/" '{print $NF}' <<< $f);
	fout=$(sed 's/bam/tsv/g' <<< $f2);

	nl=$(wc -l $COVERAGEPATH'oncomine-on-oncomine/'$fout | awk -F " " '{print $1}');

	if ! (( $nl == 37 ));
 	then
 		bedtools coverage -a ../gene-panels/Oncomine_Breast_cfDNA_v2.08212017.Designed.bed -b $f > $COVERAGEPATH'oncomine-on-oncomine/'$fout;
 	else
 		echo $fout;
 	fi;

done;

for f in $BAMPATH'Oncomine/Oncomine-files-added-020621/'*.bam;
do
	echo $f;
 	f2=$(awk -F "/" '{print $NF}' <<< $f);
 	fout=$(sed 's/bam/tsv/g' <<< $f2);

 	nl=$(wc -l $COVERAGEPATH'oncomine-on-oncomine/'$fout | awk -F " " '{print $1}');

 	if ! (( $nl == 37 ));
 	then
 		bedtools coverage -a ../gene-panels/Oncomine_Breast_cfDNA_v2.08212017.Designed.bed -b $f > $COVERAGEPATH'oncomine-on-oncomine/'$fout;
 	else
 		echo $fout;
 	fi;

done;


# 'manually' change the qiaseq bed file to replace 'chr1' with '1'
# also the qiaseq reads are mapped to GRCH38!!
cat $BEDPATH'DHS-001Z.primers-150bp_GRCh38.bed' | tr -d "chr" > $BEDPATH'DHS-001Z.primers-150bp_GRCh38.chrNames.bed'

mkdir -p $COVERAGEPATH'qiaseq-on-qiaseq';

for f in $BAMPATH'Qiaseq/'*.bam;
do
	echo $f;
	f2=$(awk -F "/" '{print $NF}' <<< $f);
	fout=$(sed 's/bam/tsv/g' <<< $f2);

	nl=$(wc -l $COVERAGEPATH'qiaseq-on-qiaseq/'$fout | awk -F " " '{print $1}');

	if ! (( $nl == 4831 ));
	then
		bedtools coverage -a ../gene-panels/DHS-001Z.primers-150bp_GRCh38.chrNames.bed -b $f > $COVERAGEPATH'qiaseq-on-qiaseq/'$fout;
	else
		echo $fout;
	fi;

done;

exit;

# check that all oncomine samples have reads on all of the regions covered by both amplicon sets
mkdir -p $COVERAGEPATH'oncomine-on-common';

for f in $BAMPATH'Oncomine/'*.bam;
do
	echo $f;
	f2=$(awk -F "/" '{print $NF}' <<< $f);
	fout=$(sed 's/bam/tsv/g' <<< $f2);

	nl=$(wc -l $COVERAGEPATH'oncomine-on-common/'$fout | awk -F " " '{print $1}');

	if ! (( $nl == 36 ));
	then
		bedtools coverage -a ../gene-panels/common_regions_sorted_merged.bed -b $f > $COVERAGEPATH'oncomine-on-common/'$fout;
	else
		echo $fout;
	fi;
done;

for f in $BAMPATH'Oncomine/Oncomine-files-added-020621/'*.bam;
do
	echo $f;
	f2=$(awk -F "/" '{print $NF}' <<< $f);
	fout=$(sed 's/bam/tsv/g' <<< $f2);

	nl=$(wc -l $COVERAGEPATH'oncomine-on-common/'$fout | awk -F " " '{print $1}');

	if ! (( $nl == 36 ));
	then
		bedtools coverage -a ../gene-panels/common_regions_sorted_merged.bed -b $f > $COVERAGEPATH'oncomine-on-common/'$fout;
	else
		echo $fout;
	fi;

done;

mkdir -p $COVERAGEPATH'qiaseq-on-common';

for f in $BAMPATH'Qiaseq/'*.bam;
do
	echo $f;
	f2=$(awk -F "/" '{print $NF}' <<< $f);
	fout=$(sed 's/bam/tsv/g' <<< $f2);

	nl=$(wc -l $COVERAGEPATH'qiaseq-on-common/'$fout | awk -F " " '{print $1}');

	if ! (( $nl == 36 ));
	then
		bedtools coverage -a ../gene-panels/common_regions_sorted_merged_38.chrNames.bed -b $f > $COVERAGEPATH'qiaseq-on-common/'$fout;
	else
		echo $fout;
	fi;

done;

exit;
# generate plots that compare the two panels
python readCoverageStats.py;

# liftover hotspots file (ucsc genome browser, manually) and remove the 'chr' so that it matches the bam file
cat $BEDPATH'Oncomine_Breast_cfDNA_v2.08212017.Hotspots_38.bed' | tr -d "chr" > $BEDPATH'Oncomine_Breast_cfDNA_v2.08212017.Hotspots_38.chrNames.bed'
