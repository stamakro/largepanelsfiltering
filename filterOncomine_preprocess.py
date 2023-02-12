import os
import pandas as pd
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


muts = set()
for i in range(data.shape[0]):
    muts.add('%s:%d' % (data.iloc[i]['Chrom'], data.iloc[i]['Position']))

muts = sorted(muts)


liftedOver = ['chr12:25245350', 'chr12:56088557', 'chr14:104780214', 'chr17:39723967', 'chr17:7670579', 'chr17:7673776',
'chr17:7673781', 'chr17:7673803', 'chr17:7674214', 'chr17:7674232', 'chr17:7674892', 'chr17:7674945', 'chr17:7674953', 'chr17:7675075', 'chr17:7675160',
'chr17:7675217', 'chr2:197402110', 'chr3:179203765', 'chr3:179210291', 'chr3:179218294','chr3:179218303', 'chr3:179218306', 'chr3:179221146',
'chr3:179234297', 'chr3:179234302', 'chr6:152011697', 'chr6:152094402', 'chr6:152098791']




vcf = data.iloc[:, [3,4,-2,5,6]]
vcf['QUAL'] = data['Quality']
vcf['FILTER'] = '.'
vcf['INFO'] = '.'
vcf.rename({'Cosmicvariant': 'ID', 'Chrom': '#CHROM', 'Position': 'POS', 'Ref': 'REF', 'Variant': 'ALT'}, axis=1, inplace=True)

for i in range(vcf.shape[0]):
    vcf.iloc[i,1] = liftedOver[muts.index(vcf.iloc[i,0] + ':' + str(vcf.iloc[i,1]))].split(':')[1]



vcfUnique = vcf.drop_duplicates('ID')
vcfUnique.to_csv('.tmp/allcalls.vcf', sep='\t', index=False)
