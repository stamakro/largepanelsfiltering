import os
import pickle
import numpy as np
import pandas as pd
from scipy.stats import beta
import sys



with open('trained_model_reduced_cosmic.pkl', 'rb') as f:
    modelDict = pickle.load(f)

clf = modelDict['model']
scaler = modelDict['scaler']


vcfLocation = sys.argv[1]
allfiles = os.listdir(vcfLocation)

featureVectorLocation = sys.argv[2]

# if this is true, then the model filter is not applied
# by default false
try:
    turnOffModel = bool(int(sys.argv[3]))
except:
    turnOffModel = False

with open(featureVectorLocation, 'rb') as f:
    fdict = pickle.load(f)

fd = fdict['mutDict']
featNames = fdict['featNames']

allLcodes = sorted(fd.keys())



sampleID = []
chromosome = []
position = []
reference = []
variant = []
xs = []

for sID in allLcodes:
    for m in fd[sID]:
        sampleID.append(sID)

        chrom, pos, ref, alt = m.split('-')

        chromosome.append(chrom)
        position.append(pos)
        reference.append(ref)
        variant.append(alt)

        # remove unused features
        featureVector = fd[sID][m]

        # apply centering and scaling to features
        xs.append(featureVector)


X = np.array(xs)
tokeep = [0, 2, 5, 6, 7, 8, 9, 11, 13, 14, 15, 16, 17, 22, 23, 24, 25, 31, 32, 33, 34, 36, 37]
X = X[:, tokeep]
# P = np.array(ps)

# allpatientIDs = np.unique(P)

featNames = np.array(['COSMIC-variant', '#var-alleles', 'QUAL', 'DP', 'VDB', 'SGB', 'RPB', 'MQB', 'MQSB', 'BQB', 'MQ0F', 'AC', 'AN', 'DP4_RF',
'DP4_RR', 'DP4_AF', 'DP4_AR', 'MQ', 'REFinGT_1', 'AD_1', 'DP_1', 'VAF_1', 'REFinGT_2', 'AD_2', 'DP_2', 'VAF_2', 'REFinGT_3',
'AD_3', 'DP_3', 'VAF_3', 'repeat_at_position', 'numAs', 'numCs', 'numGs', 'numTs', 'numNs', 'GC%', 'repeat%'])

featNames = featNames[tokeep]

X[np.isnan(X)] = 0.5


cosmicInd = np.where(X[:, 0] == 1)[0]
X = np.delete(X, 0, axis=1)
X = X[cosmicInd]
sampleID = np.array(sampleID)[cosmicInd]
chromosome = np.array(chromosome)[cosmicInd]
position = np.array(position)[cosmicInd]
reference = np.array(reference)[cosmicInd]
variant = np.array(variant)[cosmicInd]



featNames = featNames[1:]

commonRegionCalls = pd.DataFrame({'patient': sampleID,
'chr': chromosome,
'pos': position,
'ref': reference,
'alt': variant})


for i, ff in enumerate(featNames):
    commonRegionCalls[ff] = X[:, i]



predictions = clf.decision_function(scaler.transform(X))

commonRegionCalls['model_score'] = predictions
commonRegionCalls['model_score_positive'] = (predictions > 0).astype(int)

if turnOffModel:
    commonCalled = commonRegionCalls
else:
    commonCalled = commonRegionCalls[commonRegionCalls['model_score_positive'] == 1]

# --------------------

for root, _, files in os.walk(vcfLocation):
    pass
pos2cosmic = dict()
duplicates = dict()

for vcf in files:
    if len(vcf) > 1:
        with open(root + '/' + vcf) as f:
            for line in f:
                if line[0] == '#':
                    continue

                fields = line.split('\t')

                ID = fields[2]

                if ID[:4] == 'COSV':
                    pos = fields[0] + '-' + str(fields[1]) + '-' + fields[3] + '-' + fields[4]

                    if pos in pos2cosmic:
                        if pos2cosmic[pos] != ID:
                            pos2cosmic[pos] += ':' + ID
                    else:
                        pos2cosmic[pos] = ID


cosmicIDs = []
for i in range(commonCalled.shape[0]):
    p = commonCalled.iloc[i, 1] + '-' + str(commonCalled.iloc[i, 2]) + '-' + commonCalled.iloc[i, 3] +  '-' + commonCalled.iloc[i, 4]

    cosmicIDs.append(pos2cosmic[p])

assert len([k for k in cosmicIDs if ':' in k]) == 0

commonCalled['cosmicID'] = cosmicIDs

commonCalled.to_csv('./.tmp/commonCalled.csv')

allmuts = commonCalled.drop_duplicates(subset='cosmicID', inplace=False)

# chrom, pos, id, ref, alt, quality
allmuts = allmuts.iloc[:, [1,2,29,3,4,5]]
allmuts.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL']
allmuts['FILTER'] = '.'
allmuts['INFO'] = '.'

allmutsAutosome = allmuts[allmuts['#CHROM'] != 'X']
allmutsSex = allmuts[allmuts['#CHROM'] == 'X']

allmutsAutosome.sort_values(by=['#CHROM', 'POS'], inplace=True)
allmutsSex.sort_values(by='POS', inplace=True)

allmuts = pd.concat((allmutsAutosome, allmutsSex), 0)
allmuts.to_csv('./.tmp/calls.vcf', sep='\t', index=False)


sys.exit(0)
# --------------------
# error, homozygous, somatic, heterozygous
alphas = np.array([2.02388244, 19.2159724 , 3.68778916, 20.03851025])
betas = np.array([23.12416727, 6.60997908, 9.63767649, 20.72082979])
ps = np.array([0.10561397, 0.06878159, 0.29823147, 0.52737297])

zz = np.zeros((commonCalled.shape[0], 4))


for i in range(4):
    zz[:, i] = ps[i] * beta.pdf(np.array(commonCalled['VAF_2']), alphas[i], betas[i])


totals = np.sum(zz, 1)

zz = (zz.T / totals).T


commonCalled['prob_error'] = zz[:,0]
commonCalled['prob_g_heterozygous'] = zz[:, 3]
commonCalled['prob_g_homozygous'] = zz[:, 1]
commonCalled['prob_somatic'] = zz[:, 2]

variantClasses = np.array(['error', 'germ.hom.', 'somatic', 'germ.het.'])
commonCalled['label'] = variantClasses[np.argmax(zz, axis=1)]


cos2db = dict()
weirdChrom = set()
with open('../../matchedCosmicDBsnp.txt') as f:
    for line in f:
        fields = line.split('\t')
        chrom = int(fields[0].split('_')[1].split('.')[0])

        if chrom > 0 and chrom < 23:
            chrom = str(chrom)
        elif chrom == 23:
            chrom = 'X'
        else:
            weirdChrom.add(fields[0])
            continue

        db = fields[2]
        stuff = fields[-1].split('=')
        assert stuff[-1][:4] == 'COSV'
        assert stuff[-1][-1] == '\n'
        cosmic = stuff[-1][:-1]
        if cosmic in cos2db:
            assert cos2db[cosmic] == db
        else:
            cos2db[cosmic] = db

        pos = chrom + '-' + fields[1] + '-' + fields[3] + '-' + fields[4]
        pos2cosmic[pos] = cosmic

extra = dict()
for k in cos2db:
    if ':' in k:
        for kk in k.split(':'):
            extra[kk] = cos2db[k]

for k in extra:
    cos2db[k] = extra[k]



dbSNPids = []
for k in cosmicIDs:
    currentID = '.'
    if ':' in k:
        for kk in k.split(':'):
            if kk in cos2db:
                currentID = cos2db[kk]

    else:
        if k in cos2db:
            currentID = cos2db[k]

    dbSNPids.append(currentID)



commonCalled['dbSNP'] = dbSNPids


commonCalledFiltered = commonCalled[commonCalled['dbSNP'] == '.']


# --------------------------



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

for c in commonCalledFiltered['cosmicID']:
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


commonCalledFiltered['ICGC-status'] = status
commonCalledFiltered['ICGC-FR-count'] = frc
commonCalledFiltered['ICGC-FR-%'] = frp
commonCalledFiltered['ICGC-EU-count'] = euc
commonCalledFiltered['ICGC-EU-%'] = eup
commonCalledFiltered['ICGC-US-count'] = usc
commonCalledFiltered['ICGC-US-%'] = usp
commonCalledFiltered['ICGC-TOTAL-count'] = totc
commonCalledFiltered['ICGC-TOTAL-%'] = totp

import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
xx = np.linspace(0,1,150)

dens_ori = gaussian_kde(X[:,15])

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(xx, dens_ori(xx), color='k', label='start')

dens_model = gaussian_kde(commonCalled['VAF_2'])
ax.plot(xx, dens_model(xx), color='r', label='model')

dens_filter = gaussian_kde(commonCalledFiltered['VAF_2'])
ax.plot(xx, dens_filter(xx), color='g', label='model+filter')

plt.legend()

ax.set_xlabel('VAF')
ax.set_ylabel('Density')

plt.show()

print('!!!!Result not saved!!!')
