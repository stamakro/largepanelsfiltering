
python filterOncomine_preprocess.py


singularity run --bind /tudelft.net/staff-umbrella/liquidbiopsy/software/vep:/mnt,/tudelft.net/staff-umbrella/liquidbiopsy/SONIA/scripts/qiaseq-filtering/.tmp/:/home/nfs/stavrosmakrodi ~/vep/ensembl-vep_latest.sif vep -i ~/allcalls.vcf -o ~/calls_annotated.txt --most_severe --cache --dir_cache /mnt --cache_version 107 --offline --force_overwrite

tail -n+42 ./.tmp/calls_annotated.txt > ./.tmp/calls_annotated.tsv

echo 'post-processing variants...'
python filterOncomine_postprocess.py
