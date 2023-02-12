import os
import pandas as pd
import numpy as np



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
data = pd.read_csv('../../vcf/Oncomine-108-samples-CADD.csv', index_col=0)

data = data[data['Sample'].isin(samples)]

allSubjects = data['subject'].unique()

data = data[data['Type'] == 'SNP']


data = data[data['Quality'] > 20]

data = data[data['Allele.Call'] == 'Heterozygous']


data = data[data['Allele.Mol.Cov'] > 3]

data.reset_index(inplace=True)
#remove call in homopolymer region
data.drop(index=16, inplace=True)
print(data.shape[0])

assert data.groupby('Allele.Call').ngroups == 1

data.drop(['Filter.1', 'LOD', 'Gene.ID', 'Region.Name',
       'Subset.Of', 'VCF.Position', 'VCF.Ref', 'VCF.Variant', 'Read.Cov', 'Filter', 'Allele.Call',
       'Allele.Read.Cov', 'Allele.Read.Freq', 'Filter.2',
        'Filter.3', 'Allele.Mol.Freq', 'Filter.4',
       'Strand.Bias', 'Filter.5', 'Common.Signal.Shift', 'Filter.6',
       'Reference.Signal.Shift', 'Filter.7', 'Variant.Signal.Shift',
       'Filter.8', 'Relative.Read.Quality', 'Filter.9', 'HP.Length',
       'Filter.10', 'Context.Error.', 'Filter.11', 'Context.Error..1',
       'Filter.12', 'Context.Strand.Bias', 'Filter.13', 'Sample.Name',
       'Barcode', 'Run.Name', 'Location'], axis=1, inplace=True)


vep = pd.read_table('.tmp/calls_annotated.tsv')

# 1 consequence per variant
assert vep.shape[0] == vep['#Uploaded_variation'].unique().shape[0]

cosmic2cons = vep[['#Uploaded_variation', 'Consequence']].set_index('#Uploaded_variation').to_dict()['Consequence']

data['VEP_consequence'] = data['Cosmicvariant'].map(cosmic2cons)

ind = np.where([data['Position'] == 7573897])[1]
data.iloc[ind, -1] = 'intron_variant'

ind = np.where([data['Position'] == 7578210])[1]
data.iloc[ind, -1] = 'synonymous_variant'

ind = np.where([data['Position'] == 7578263])[1]
data.iloc[ind, -1] = 'stop_gained'

ind = np.where([data['Position'] == 7578393])[1]
data.iloc[ind, -1] = 'missense_variant'

desiredConsequences = {'missense_variant', 'start_lost', 'stop_lost', 'stop_gained', 'splice_acceptor_variant', 'splice_donor_variant', 'splice_donor_5th_base_variant'}
dataFiltered = data[data['VEP_consequence'].isin(desiredConsequences)]


icgcInfo = pd.read_table('/tudelft.net/staff-umbrella/liquidbiopsy/human-genome/icgc-brca-somatic-snvs/cosmic_in_qiaseq_header.vcf')
icgcInfo.set_index('ID', inplace=True)
icgcInfo.drop(['#chrom', 'position', 'ref', 'alt', 'Mutation strand'], axis=1, inplace=True)


# ---
frc = []
frp = []
euc = []
eup = []
usc = []
usp = []
totc = []
totp = []
status = []

for c in dataFiltered['Cosmicvariant']:
    try:
        row = icgcInfo.loc[c]
        stats = list(row)
    except KeyError:
        stats = ['Unknown', 0, 0, 0, 0, 0, 0, 0, 0]

    status.append(stats[0])
    frc.append(stats[1])
    frp.append(stats[2])
    euc.append(stats[3])
    eup.append(stats[4])
    usc.append(stats[5])
    usp.append(stats[6])
    totc.append(stats[7])
    totp.append(stats[8])



dataFiltered['ICGC-status'] = status
dataFiltered['ICGC-FR-count'] = frc
dataFiltered['ICGC-FR-%'] = frp
dataFiltered['ICGC-EU-count'] = euc
dataFiltered['ICGC-EU-%'] = eup
dataFiltered['ICGC-US-count'] = usc
dataFiltered['ICGC-US-%'] = usp
dataFiltered['ICGC-TOTAL-count'] = totc
dataFiltered['ICGC-TOTAL-%'] = totp

# manually fix for the two variants not in cosmic, one is not present
dataFiltered.loc[23,'ICGC-TOTAL-count'] = 16

dataFiltered = dataFiltered[dataFiltered['ICGC-TOTAL-count'] > 0]
dataFiltered.to_csv('oncomine_calls_filtered.csv')
