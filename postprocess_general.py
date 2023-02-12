import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta
import sys

def compareImpact(x):
    x = set(x)
    if 'HIGH' in x:
        return 'HIGH'
    if 'MODERATE' in x:
        return 'MODERATE'
    if 'LOW' in x:
        return 'LOW'

    assert 'MODIFIER' in x
    return 'MODIFIER'

commonCalled = pd.read_csv('./.tmp/commonCalled.csv', index_col=0)

vep = pd.read_table('.tmp/calls_annotated.tsv')

# 1 consequence per variant
assert vep.shape[0] == vep['#Uploaded_variation'].unique().shape[0]

cosmic2cons = vep[['#Uploaded_variation', 'Consequence']].set_index('#Uploaded_variation').to_dict()['Consequence']

commonCalled['VEP_consequence'] = commonCalled['cosmicID'].map(cosmic2cons)


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

        # pos = chrom + '-' + fields[1] + '-' + fields[3] + '-' + fields[4]
        # pos2cosmic[pos] = cosmic

extra = dict()
for k in cos2db:
    if ':' in k:
        for kk in k.split(':'):
            extra[kk] = cos2db[k]

for k in extra:
    cos2db[k] = extra[k]


commonCalled['dbSNP'] = commonCalled['cosmicID'].map(cos2db)



# dbSNPids = []
# for k in cosmicIDs:
#     currentID = '.'
#     if ':' in k:
#         for kk in k.split(':'):
#             if kk in cos2db:
#                 currentID = cos2db[kk]
#
#     else:
#         if k in cos2db:
#             currentID = cos2db[k]
#
#     dbSNPids.append(currentID)
#
#
#
# commonCalled['dbSNP'] = dbSNPids


# --------------------------


# commonCalledFiltered = commonCalled[commonCalled['VEP_impact'] != 'MODIFIER']


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

for c in commonCalled['cosmicID']:
    try:
        row = icgcInfo.loc[c]
        stats = list(row)
    except KeyError:
        if ':' not in c:
            stats = ['Unknown', 0, 0, 0, 0, 0, 0, 0, 0]
        else:
            stats = ['Unknown', 0, 0, 0, 0, 0, 0, 0, 0]
            for cc in c.split(':'):
                if cc in icgcInfo.index:
                    row = icgcInfo.loc[cc]
                    stats = list(row)
                    break



    status.append(stats[0])
    frc.append(stats[1])
    frp.append(stats[2])
    euc.append(stats[3])
    eup.append(stats[4])
    usc.append(stats[5])
    usp.append(stats[6])
    totc.append(stats[7])
    totp.append(stats[8])


commonCalled['ICGC-status'] = status
commonCalled['ICGC-FR-count'] = frc
commonCalled['ICGC-FR-%'] = frp
commonCalled['ICGC-EU-count'] = euc
commonCalled['ICGC-EU-%'] = eup
commonCalled['ICGC-US-count'] = usc
commonCalled['ICGC-US-%'] = usp
commonCalled['ICGC-TOTAL-count'] = totc
commonCalled['ICGC-TOTAL-%'] = totp

desiredConsequences = {'missense_variant', 'start_lost', 'stop_lost', 'stop_gained', 'splice_acceptor_variant', 'splice_donor_variant', 'splice_donor_5th_base_variant'}

commonCalledFiltered = commonCalled[commonCalled['VEP_consequence'].isin(desiredConsequences)]

commonCalledFiltered = commonCalledFiltered[commonCalledFiltered['ICGC-TOTAL-count'] > 0]


# large enough vaf + dbsnp entry
s = commonCalledFiltered[commonCalledFiltered['label'] == 'germ.het.']['dbSNP'].notna()
if len(s[s].index):
    commonCalledFiltered.drop(s[s].index, inplace=True)

# large enough vaf + dbsnp entry
s = commonCalledFiltered[commonCalledFiltered['label'] == 'germ.hom.']['dbSNP'].notna()
if len(s[s].index):
    commonCalledFiltered.drop(s[s].index, inplace=True)



# import matplotlib.pyplot as plt
# from scipy.stats import gaussian_kde
# xx = np.linspace(0,1,150)
#
# dens_ori = gaussian_kde(commonCalled['VAF_2'])
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
#
# ax.plot(xx, dens_ori(xx), color='k', label='start')
#
# dens_model = gaussian_kde(commonCalled['VAF_2'])
# ax.plot(xx, dens_model(xx), color='r', label='model')
#
# dens_filter = gaussian_kde(commonCalledFiltered['VAF_2'])
# ax.plot(xx, dens_filter(xx), color='g', label='model+filter')
#
# plt.legend()
#
# ax.set_xlabel('VAF')
# ax.set_ylabel('Density')
#
# plt.show()

try:
    fileLocation = sys.argv[1]
    commonCalledFiltered.to_csv(fileLocation)
except FileNotFoundError:
    print('!!!!Result not saved!!!')
