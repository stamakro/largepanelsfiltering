import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

directory = '../../vcf/Oncomine/'
files = os.listdir(directory)


# delete 2 samples that were run twice
tobedel = []
for i, f in enumerate(files):
    if f[:5] == 'L7063' and '16_R' in f or f[:5] == 'L7041' and '15_R' in f:
        tobedel.append(i)

for ii in sorted(tobedel)[::-1]:
    del files[ii]

samples = [f[:5] for f in files]


# first look at all 108 samples
data = pd.read_csv('../../vcf/Oncomine-108-samples-CADD.csv')
print('start:\t%d variants' % data.shape[0])

data = data[data['Sample'].isin(samples)]
print('common qiaseq/oncomine samples:\t%d' % data.shape[0])

allSubjects = data['subject'].unique()

data = data[data['Type'] == 'SNP']
print('only snps:\t%d' % data.shape[0])

# data = data[data['Mol.Coverage'] > 100]
    # print('coverage:\t%d' % data.shape[0])

data.drop(['Filter.1', 'Gene.ID', 'Region.Name',
       'Subset.Of', 'VCF.Position', 'VCF.Ref', 'VCF.Variant', 'Read.Cov', 'Filter',
       'Allele.Read.Cov', 'Allele.Read.Freq', 'Filter.2',
	'Filter.3', 'Allele.Mol.Freq', 'Filter.4',
       'Strand.Bias', 'Filter.5', 'Common.Signal.Shift', 'Filter.6',
       'Reference.Signal.Shift', 'Filter.7', 'Variant.Signal.Shift',
       'Filter.8', 'Relative.Read.Quality', 'Filter.9', 'HP.Length',
       'Filter.10', 'Context.Error.', 'Filter.11', 'Context.Error..1',
       'Filter.12', 'Context.Strand.Bias', 'Filter.13', 'Sample.Name',
       'Barcode', 'Run.Name', 'Location'], axis=1, inplace=True)




dataHeterozygous = data[data['Allele.Call'] == 'Heterozygous']
print('remove no-calls:\t%d' % dataHeterozygous.shape[0])

dataQuality = data[data['Quality'] > 20]
print('remove low quality:\t%d' % dataQuality.shape[0])

dataQualityHeterozygous = dataQuality[dataQuality['Allele.Call'] == 'Heterozygous']
print('Q>20 & heterozygous:\t%d' % dataQualityHeterozygous.shape[0])




data3reads = data[data['Allele.Mol.Cov'] > 3]
print('3 var reads: %d' % data3reads.shape[0])




new = set(dataQualityHeterozygous[['Sample', 'Chrom', 'Position', 'Ref', 'Variant']].astype(str).apply(lambda x: ':'.join(x), axis=1))
old = set(data3reads[['Sample', 'Chrom', 'Position', 'Ref', 'Variant']].astype(str).apply(lambda x: ':'.join(x), axis=1))


data3reads.to_csv('../../vcf/oncomine-filter-variations/oldway-snv-molcov100-allelemolcov3.csv')
dataQualityHeterozygous.to_csv('../../vcf/oncomine-filter-variations/newway-snv-molcov100-qual20-callheterozygous.csv')





sys.exit(0)

data.reset_index(inplace=True)
#remove call in homopolymer region
data.drop(index=16, inplace=True)
print(data.shape[0])

assert data.groupby('Allele.Call').ngroups == 1
