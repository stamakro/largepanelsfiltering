import os
import numpy as np
import pandas as pd
import sys
from sklearn.metrics import precision_recall_curve, average_precision_score, roc_curve, precision_recall_fscore_support, recall_score, accuracy_score, matthews_corrcoef
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC, LinearSVC
from sklearn.model_selection import StratifiedKFold
import pickle
import matplotlib.pyplot as plt
from scipy.stats import beta

def methodOfMoments(m, v):
    # estimate a, b params of beta distribution
    # using the method of moments
    # m: mean, v: variance of sample
    # returns estimated a, b
    pp = (m * (1 - m) / v) - 1

    a = m * pp
    b = (1 - m) * pp

    assert a > 0 and b > 0

    return a, b


def weightedMeanAndVariance(x, w):
    # calculate the weighted sample mean and variance of data x
    # with weights given by w
    m = np.average(x, weights=w)

    v = np.average( (x - m) ** 2, weights=w)
    n = x.shape[0]

    v = n * v / (n - 1)

    return m, v



def mixtureOfBetasEM(x, maxiter=5000, eps=1e-4):
    # components initialized as follows:
    # component 1: errors
    # component 2: germline homozygous
    # component 3: somatic heterozygous
    # component 4: germline heterozygous

    N = x.shape[0]
    Ncomp = 4
    alphasOld = np.array([1.0, 20.0, 8.0, 25.])
    betasOld = np.array([30.0, 1.0, 20.0, 25.])
    #pOld = np.array([1./3, 1./3, 1./3])
    pOld = np.ones(Ncomp) / Ncomp

    LLold = - np.inf
    converged = False

    zz = np.zeros((N, Ncomp))

    for counter in range(maxiter):

        # E step
        for i in range(Ncomp):
            zz[:, i] = pOld[i] * beta.pdf(x, alphasOld[i], betasOld[i])


        totals = np.sum(zz, 1)

        zz = (zz.T / totals).T


        # M step
        alphasNew = np.zeros(alphasOld.shape)
        betasNew = np.zeros(betasOld.shape)
        pNew = np.mean(zz, 0)

        for i in range(Ncomp):
            mi, vi = weightedMeanAndVariance(x, zz[:, i])

            alphasNew[i], betasNew[i] = methodOfMoments(mi, vi)

        zz_thr = np.zeros(zz.shape)
        for ss in range(N):
            zz_thr[ss, np.argmax(zz[ss])] = 1.

        LLnew = 0.
        for ii in range(Ncomp):
            LLnew += np.sum(zz_thr[:, i] * (beta.logpdf(x, alphasNew[i], betasNew[i]) + np.log(pNew[i])))

        converged = True
        if LLnew - LLold > eps:
            converged = False

        alphasOld = alphasNew
        betasOld = betasNew
        pOld = pNew
        LLold = LLnew

        print('Iteration %4d, LL = %.4f' % (counter, LLnew))
        if converged:
            break

    return alphasNew, betasNew, pNew, zz


def plotMixture(x, alphas, betas, ps, histogramBins=30, ax=None, linestyle='--', showHistogram=True):
    Ncomp = alphas.shape[0]

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if showHistogram:
        ax.hist(x, bins=histogramBins, color='C0', edgecolor='k', density=True)

    xx = np.linspace(0,1,500)
    for i in range(Ncomp):
        ax.plot(xx, ps[i] * beta.pdf(xx, alphas[i], betas[i]), color='C'+str(i+1), linestyle=linestyle)


    return ax







def informedness(Ytrue, Ypred):
    # aka youden's J
    return np.sum(recall_score(Ytrue, Ypred, average=None)) - 1


# oncomine
directory = sys.argv[1]
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
    muts.add('%s:%d-%d' % (data.iloc[i]['Chrom'], data.iloc[i]['Position'], data.iloc[i]['Position']))

muts = sorted(muts)

# for m in muts:
#     print(m)

liftedOver = ['chr12:25245350', 'chr12:56088557', 'chr14:104780214', 'chr17:39723967', 'chr17:7670579', 'chr17:7673776',
'chr17:7673781', 'chr17:7673803', 'chr17:7674214', 'chr17:7674232', 'chr17:7674892', 'chr17:7674945', 'chr17:7674953', 'chr17:7675075', 'chr17:7675160',
'chr17:7675217', 'chr2:197402110', 'chr3:179203765', 'chr3:179210291', 'chr3:179218294','chr3:179218303', 'chr3:179218306', 'chr3:179221146',
'chr3:179234297', 'chr3:179234302', 'chr6:152011697', 'chr6:152094402', 'chr6:152098791']


liftoverDict = {}
for m19, m38 in zip(muts, liftedOver):
    liftoverDict[m19.split('-')[0]] = m38

positions = []
for i in range(data.shape[0]):
    m = data.iloc[i]['Chrom'] + ':' + str(data.iloc[i]['Position'])
    m38 = liftoverDict[m]
    positions.append(int(m38.split(':')[1]))


data['Position_grch38'] = positions


patientInfo = pd.read_csv('../../Lcode-Oncomine-Qiaseq-per-subject.csv')
subject2QiaseqL = pd.Series(patientInfo['Lcode qiaseq'].values, index=patientInfo['Subject']).to_dict()
qiaseqL2subject = dict()
for s in subject2QiaseqL:
    qiaseqL2subject[subject2QiaseqL[s]] = s

# for sub in subject2QiaseqL:
#     assert os.path.exists('../../vcf/Qiaseq-common-regions/'+ subject2QiaseqL[sub] + '.vcf')



mutationDict = {}
for sub in allSubjects:
    mutationDict[sub] = set()

for i in range(data.shape[0]):
    subject = data.iloc[i]['subject']
    chr = data.iloc[i]['Chrom'][3:]
    pos = data.iloc[i]['Position_grch38']

    mutationDict[subject].add(chr + ':' + str(pos))



vcfLocation = sys.argv[2]
allfiles = os.listdir(vcfLocation)

featureVectorsCommonPath = sys.argv[3]
with open(featureVectorsCommonPath, 'rb'  ) as f:
    fdict = pickle.load(f)

fd = fdict['mutDict']
featNames = fdict['featNames']

allLcodes = sorted(fd.keys())



qiaseqLs = list(patientInfo['Lcode qiaseq'])
tobedel = []
for i, ff in enumerate(allLcodes):
    index = qiaseqLs.index(ff)
    if patientInfo.iloc[index]['Lcode oncomine'] == 'x':
        tobedel.append(i)

for tt in sorted(tobedel, reverse=True):
    del allLcodes[tt]


oncomineCalledMutations = dict()


ys = []
xs = []
ps = []


for i, lcode in enumerate(allLcodes):

    oncomineCalledMutations[qiaseqL2subject[lcode]] = set()

    qiaseqCalls = fd[lcode]

    for call in qiaseqCalls:
        chr, pos, ref, alt = call.split('-')

        if chr + ':' + pos in mutationDict[qiaseqL2subject[lcode]]:
            ys.append(1)
            oncomineCalledMutations[qiaseqL2subject[lcode]].add(chr + ':' + pos)

        else:
            ys.append(0)

        ps.append(lcode)
        xs.append(qiaseqCalls[call])

print(np.sum(ys))

#dd = pd.DataFrame({'patient': ps, 'call': ys})
#dd.to_csv('../../calls_all.csv')

differenceDict = dict()
for k in mutationDict:
    if len(mutationDict[k]) > 0:
        differenceDict[k] = mutationDict[k].difference(oncomineCalledMutations[k])

calledByQ = []
missedByQ = []

for k in mutationDict:
    ddd = data[data['subject'] == k]
    if len(mutationDict[k]) > 0:
        for m in mutationDict[k]:
            chrom, pos = m.split(':')

            chrom = 'chr' + str(chrom)
            pos = int(pos)
            ddd2 = ddd[ddd['Chrom'] == chrom]

            ddd3 = ddd2[ddd2['Position_grch38'] == pos]
            assert ddd3.shape[0] == 1

            myVAF = ddd3.iloc[0]['Frequency']

            if m in differenceDict[k]:
                missedByQ.append(myVAF)
                #print('%s\t%s\t%f' % (k, m, myVAF))
            else:
                calledByQ.append(myVAF)


X = np.array(xs)
y = np.array(ys)
P = np.array(ps)

allpatientIDs = np.unique(P)

featNames = np.array(['COSMIC-variant', '#var-alleles', 'QUAL', 'DP', 'VDB', 'SGB', 'RPB', 'MQB', 'MQSB', 'BQB', 'MQ0F', 'AC', 'AN', 'DP4_RF',
'DP4_RR', 'DP4_AF', 'DP4_AR', 'MQ', 'REFinGT_1', 'AD_1', 'DP_1', 'VAF_1', 'REFinGT_2', 'AD_2', 'DP_2', 'VAF_2', 'REFinGT_3',
'AD_3', 'DP_3', 'VAF_3', 'repeat_at_position', 'numAs', 'numCs', 'numGs', 'numTs', 'numNs', 'GC%', 'repeat%'])

X[np.isnan(X)] = 0.5

keep = np.ones(X.shape[1], int)
keep[np.where(np.std(X, axis=0) == 0)[0]] = 0

keep[np.where(featNames == 'VDB')[0][0]] = 0

keep[np.where(featNames == 'DP')[0][0]] = 0


ind = [i for i, fn in enumerate(featNames) if '_1' in fn or '_3' in fn]
keep[ind] = 0

tokeep = np.where(keep)[0]

print(keep.shape)
print(tokeep)
#sys.exit(0)
X = X[:, tokeep]
featNames = featNames[tokeep]


assert P.shape[0] == X.shape[0]

y[y == 0] = -1

nFolds = 4

paramValues = np.array([1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 5.0, 10., 20., 50.])

innerPrecision = np.zeros((nFolds, paramValues.shape[0]))
innerRecall = np.zeros(innerPrecision.shape)
innerF1 = np.zeros(innerPrecision.shape)
innerMCC = np.zeros(innerPrecision.shape)
innerAccuracy = np.zeros(innerPrecision.shape)
innerInformedness = np.zeros(innerPrecision.shape)

innerPrecision_RBF = np.zeros((nFolds, paramValues.shape[0], paramValues.shape[0]))
innerRecall_RBF = np.zeros(innerPrecision_RBF.shape)
innerF1_RBF = np.zeros(innerPrecision_RBF.shape)
innerMCC_RBF = np.zeros(innerPrecision_RBF.shape)
innerAccuracy_RBF = np.zeros(innerPrecision_RBF.shape)
innerInformedness_RBF = np.zeros(innerPrecision_RBF.shape)

cosmicInd = np.where(X[:, 0] == 1)[0]
X = X[cosmicInd]
y = y[cosmicInd]

X = np.delete(X, 0, axis=1)

featNames = featNames[1:]

cv = StratifiedKFold(n_splits=nFolds, shuffle=True, random_state=42)

for fold, (trn, val) in enumerate(cv.split(X, y)):
    print('fold %d' % fold)
    Xtrain = X[trn]
    Xval = X[val]
    ytrain = y[trn]
    yval = y[val]

    ss = StandardScaler()
    Xss = ss.fit_transform(Xtrain)

    Xssval = ss.transform(Xval)


    for i, C in enumerate(paramValues):

        clf = LinearSVC(C=C, penalty='l2', loss='hinge', max_iter=500000)
        clf.fit(Xss, ytrain)


        predictions = clf.predict(Xssval)

        innerPrecision[fold, i], innerRecall[fold, i], innerF1[fold, i], _ = precision_recall_fscore_support(yval, predictions, average='binary')
        innerAccuracy[fold, i] = accuracy_score(yval, predictions)
        innerMCC[fold, i] = matthews_corrcoef(yval, predictions)
        innerInformedness[fold, i] = informedness(yval, predictions)

    for i, C in enumerate(paramValues):
        for j, g in enumerate(paramValues):
            clf = SVC(C=C, kernel='rbf', gamma=g)
            clf.fit(Xss, ytrain)


            predictions = clf.predict(Xssval)

            innerPrecision_RBF[fold, i, j], innerRecall_RBF[fold, i, j], innerF1_RBF[fold, i, j], _ = precision_recall_fscore_support(yval, predictions, average='binary')
            innerAccuracy_RBF[fold, i, j] = accuracy_score(yval, predictions)
            innerMCC_RBF[fold, i, j] = matthews_corrcoef(yval, predictions)
            innerInformedness_RBF[fold, i, j] = informedness(yval, predictions)



perfLinear = np.mean(innerInformedness, axis=0)
perfRBF = np.mean(innerInformedness_RBF, axis=0)


bestLinear = np.max(perfLinear)
bestRBF = np.max(perfRBF)

bestModelInd = np.argmax([bestLinear, bestRBF])

if bestModelInd == 0:
    # lin
    ii = np.argmax(perfLinear)
    clf = LinearSVC(C=paramValues[ii], penalty='l2', loss='hinge', max_iter=500000)
else:
    # rbf
    ii = np.argmax(perfRBF)
    ind1, ind2 = np.unravel_index(ii, (paramValues.shape[0], paramValues.shape[0]))

    clf = SVC(C=paramValues[ind1], kernel='rbf', gamma=paramValues[ind2])

sys.exit(0)
ss2 = StandardScaler()
XCS = ss2.fit_transform(X)

clf.fit(XCS, y)

modelSaveLocation = sys.argv[4]
with open(modelSaveLocation, 'wb') as f:
    pickle.dump({'model': clf, 'scaler': ss2}, f)


# show weights of best linear model
ii = np.argmax(perfLinear)
clf_lin = LinearSVC(C=paramValues[ii], penalty='l2', loss='hinge', max_iter=500000)

clf_lin.fit(XCS, y)
iii = np.argsort(clf_lin.coef_[0])

for ff, ww in zip(featNames[iii], clf_lin.coef_[0][iii]):
    print('%20s\t%1.4f' % (ff, ww))


########################## start application to other regions ##########################
featureVectorsUniquePath = sys.argv[5]
with open(featureVectorsUniquePath, 'rb'  ) as f:
    fdict = pickle.load(f)


fd = fdict['mutDict']
featNamesTest = fdict['featNames']

allLcodes = sorted(fd.keys())

sampleID = []
chromosome = []
position = []
reference = []
variant = []
xx = []

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
        xx.append(featureVector[tokeep])


xx = np.array(xx)

cosmicInd = np.where(xx[:, 0] == 1)[0]

xx = np.delete(xx, 0, axis=1)
xx = xx[cosmicInd]
sampleID = np.array(sampleID)[cosmicInd]
chromosome = np.array(chromosome)[cosmicInd]
position = np.array(position)[cosmicInd]
reference = np.array(reference)[cosmicInd]
variant = np.array(variant)[cosmicInd]

uniqueRegionCalls = pd.DataFrame({'patient': sampleID,
'chr': chromosome,
'pos': position,
'ref': reference,
'alt': variant})


for i, ff in enumerate(featNames):
    uniqueRegionCalls[ff] = xx[:, i]

print(np.sum(np.isnan(xx)))

xx[np.isnan(xx)] = 0.5
predictions = clf.decision_function(ss2.transform(xx))

uniqueRegionCalls['model_score'] = predictions
uniqueRegionCalls['model_score_positive'] = (predictions > 0).astype(int)

uniqueRegionCallsPath = sys.argv[6]
uniqueRegionCalls.to_csv(uniqueRegionCallsPath)


called = uniqueRegionCalls[uniqueRegionCalls['model_score'] > 0]

import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import gaussian_kde, norm

vafs = np.array(called['VAF_2'])

vafDensity = gaussian_kde(vafs)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

vv = np.linspace(0,1,1000)

ax.hist(called['VAF_2'], bins=30, edgecolor='k', color='C0', density=True)
ax.plot(vv, vafDensity(vv), color='C1')

mix = GaussianMixture(n_components=3)
mix.fit(vafs.reshape(-1,1))

for i, (mm, v, ww) in enumerate(zip(mix.means_, mix.covariances_, mix.weights_)):
    sigma = np.sqrt(v[0])
    print('%d: %.3f' % (i, mm))

    ax.plot(vv, norm.pdf(vv, mm, sigma) * ww, color='C'+str(i+2))

somPeak = np.argmin(mix.means_)

predictedComponent = mix.predict(vafs.reshape(-1,1))


#called['somatic'] = (predictedComponent == somPeak).astype(int)

#calledAndSomatic = called[called['somatic'] == 1]



alphasNew, betasNew, pNew, zz = mixtureOfBetasEM(np.array(called['VAF_2']))
ax = plotMixture(np.array(called['VAF_2']), alphasNew, betasNew, pNew, histogramBins=50)

alphasOld = np.array([1.0, 20.0, 8.0, 25.])
betasOld = np.array([30.0, 1.0, 20.0, 25.])
pOld = np.ones(4) / 4.

plotMixture(np.array(called['VAF_2']), alphasOld, betasOld, pOld, histogramBins=50, ax=ax, linestyle='-', showHistogram=False)


variantClasses = np.array(['error', 'germ.hom.', 'somatic', 'germ.het.'])
called['prob_error'] = zz[:,0]
called['prob_g_heterozygous'] = zz[:, 3]
called['prob_g_homozygous'] = zz[:, 1]
called['prob_somatic'] = zz[:, 2]

called['label'] = variantClasses[np.argmax(zz, axis=1)]

vcfLocationUnique = sys.argv[7]
for root, _, files in os.walk(vcfLocationUnique):
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
                            print('duplicate! %s' % pos)
                            pos2cosmic[pos] += ':' + ID
                    else:
                        pos2cosmic[pos] = ID




cosmicIDs = []
for i in range(called.shape[0]):
    p = called.iloc[i, 1] + '-' + str(called.iloc[i, 2]) + '-' +called.iloc[i, 3] + '-' + called.iloc[i, 4]

    cosmicIDs.append(pos2cosmic[p])

called['cosmicID'] = cosmicIDs

called.to_csv('./.tmp/uniqueCalled.csv')

allmuts = called.drop_duplicates(subset='cosmicID', inplace=False)

# chrom, pos, id, ref, alt, quality
allmuts = allmuts.iloc[:, [1,2,34,3,4,5]]
allmuts.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL']
allmuts['FILTER'] = '.'
allmuts['INFO'] = '.'

allmutsAutosome = allmuts[allmuts['#CHROM'] != 'X']
allmutsSex = allmuts[allmuts['#CHROM'] == 'X']

allmutsAutosome.sort_values(by=['#CHROM', 'POS'], inplace=True)
allmutsSex.sort_values(by='POS', inplace=True)

allmuts = pd.concat((allmutsAutosome, allmutsSex), 0)
allmuts.to_csv('./.tmp/calls.vcf', sep='\t', index=False)


#
# fig = plt.figure()
# ax = fig.add_subplot(111)
#
# ax.boxplot((Xtest[predictions == 1, 20], Xtest[predictions == -1, 20]))
#
# ax.set_ylabel('VAF', fontsize=25)
# ax.set_xticklabels(['KEPT\n'+str(np.sum(predictions==1)), 'FILTERED OUT\n'+str(np.sum(predictions==-1))])
#
#
# plt.show()
