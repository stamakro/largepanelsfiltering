import os
import numpy as np
import pandas as pd
import sys
from scipy.stats import mannwhitneyu, ttest_ind, spearmanr, gaussian_kde
from sklearn.metrics import precision_recall_curve, average_precision_score, roc_curve, precision_recall_fscore_support, recall_score, accuracy_score, matthews_corrcoef
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC, LinearSVC
from sklearn.model_selection import StratifiedKFold
import pickle

def informedness(Ytrue, Ypred):
    # aka youden's J
    return np.sum(recall_score(Ytrue, Ypred, average=None)) - 1



thresholdTypeInd = int(sys.argv[1])
# 0: 0, 1: max rc s.t. pr = 1
assert thresholdTypeInd > -1 and thresholdTypeInd < 3

print('!!!!!!! %d !!!!!!!' % thresholdTypeInd)

# oncomine
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

for sub in subject2QiaseqL:
    assert os.path.exists('../../vcf/Qiaseq-common-regions/'+ subject2QiaseqL[sub] + '.vcf')



mutationDict = {}
for sub in allSubjects:
    mutationDict[sub] = set()

for i in range(data.shape[0]):
    subject = data.iloc[i]['subject']
    chr = data.iloc[i]['Chrom'][3:]
    pos = data.iloc[i]['Position_grch38']

    mutationDict[subject].add(chr + ':' + str(pos))


vcfLocation = '../../vcf/Qiaseq-common-regions/'
allfiles = os.listdir(vcfLocation)

with open('../../featureVectorsCommon.pkl', 'rb'  ) as f:
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

dd = pd.DataFrame({'patient': ps, 'call': ys})
dd.to_csv('../../calls_all.csv')

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

featNames = ['COSMIC-variant', '#var-alleles', 'QUAL', 'DP', 'VDB', 'SGB', 'RPB', 'MQB', 'MQSB', 'BQB', 'MQ0F', 'AC', 'AN', 'DP4_RF',
'DP4_RR', 'DP4_AF', 'DP4_AR', 'MQ', 'REFinGT_1', 'AD_1', 'DP_1', 'VAF_1', 'REFinGT_2', 'AD_2', 'DP_2', 'VAF_2', 'REFinGT_3',
'AD_3', 'DP_3', 'VAF_3', 'repeat_at_position', 'numAs', 'numCs', 'numGs', 'numTs', 'numNs', 'GC%', 'repeat%']

X[np.isnan(X)] = 0.5


tobedelFeats = np.where(np.std(X, axis=0) == 0)[0]
X = np.delete(X, tobedelFeats, axis=1)

featNames = np.delete(featNames, tobedelFeats)

ind = np.where(featNames == 'VDB')[0][0]

X = np.delete(X, ind, axis=1)
featNames = np.delete(featNames, ind)

ind = [i for i, fn in enumerate(featNames) if '_1' in fn or '_3' in fn]

X = np.delete(X, ind, axis=1)
featNames = np.delete(featNames, ind)

pvalsU = np.zeros(featNames.shape[0])
pvalsT = np.zeros(featNames.shape[0])

for i in range(X.shape[1]):
    pvalsU[i] = mannwhitneyu(X[y==0, i], X[y==1, i])[1]
    pvalsT[i] = ttest_ind(X[y==0, i], X[y==1, i])[1]


assert P.shape[0] == X.shape[0]


y[y == 0] = -1

nInnerFolds = 4


paramValues = np.array([1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 5.0, 10., 20., 50.])


innerPrecision = np.zeros((allpatientIDs.shape[0], nInnerFolds, paramValues.shape[0]))
innerRecall = np.zeros(innerPrecision.shape)
innerF1 = np.zeros(innerPrecision.shape)
innerMCC = np.zeros(innerPrecision.shape)
innerAccuracy = np.zeros(innerPrecision.shape)
innerInformedness = np.zeros(innerPrecision.shape)
threshold = np.zeros(innerPrecision.shape)

innerPrecision_RBF = np.zeros((allpatientIDs.shape[0], nInnerFolds, paramValues.shape[0], paramValues.shape[0]))
innerRecall_RBF = np.zeros(innerPrecision_RBF.shape)
innerF1_RBF = np.zeros(innerPrecision_RBF.shape)
innerMCC_RBF = np.zeros(innerPrecision_RBF.shape)
innerAccuracy_RBF = np.zeros(innerPrecision_RBF.shape)
innerInformedness_RBF = np.zeros(innerPrecision_RBF.shape)
threshold_RBF = np.zeros(innerPrecision_RBF.shape)

trainingSetSize = np.zeros(innerPrecision.shape[0])
testSetSize = np.zeros(innerPrecision.shape[0])
nMutations = np.zeros(innerPrecision.shape[0])
nMutationsCalled = np.zeros(innerPrecision.shape[0])
bestModel = np.zeros(innerPrecision.shape[0])
precision = np.zeros(innerPrecision.shape[0])
recall = np.zeros(innerPrecision.shape[0])
f1 = np.zeros(innerPrecision.shape[0])
mcc = np.zeros(innerPrecision.shape[0])
acc = np.zeros(innerPrecision.shape[0])
informednessYouden = np.zeros(innerPrecision.shape[0])

sys.exit(0)

for pInd, patient in enumerate(allpatientIDs):

    testInd = np.where(P == patient)[0]
    trainInd = np.where(P != patient)[0]

    trainingSetSize[pInd] = trainInd.shape[0]
    testSetSize[pInd] = testInd.shape[0]

    XtrainOuter = X[trainInd]
    ytrainOuter = y[trainInd]

    Xtest = X[testInd]
    ytest = y[testInd]

    nMutations[pInd] = np.sum(ytest == 1)

    cv = StratifiedKFold(n_splits=nInnerFolds, shuffle=True, random_state=42)
    for fold, (trn, val) in enumerate(cv.split(XtrainOuter, ytrainOuter)):

        print('Patient %d/%d, fold %d/%d' % (pInd, allpatientIDs.shape[0], fold, nInnerFolds))

        Xtrain = XtrainOuter[trn]
        Xval = XtrainOuter[val]
        ytrain = ytrainOuter[trn]
        yval = ytrainOuter[val]

        ss = StandardScaler()
        Xss = ss.fit_transform(Xtrain)

        Xssval = ss.transform(Xval)

        for i, C in enumerate(paramValues):
            clf = LinearSVC(C=C, penalty='l2', loss='hinge', max_iter=500000)
            clf.fit(Xss, ytrain)

            if thresholdTypeInd == 0:
                predictions = clf.predict(Xssval)
            else:
                pr, rc, thres = precision_recall_curve(ytrain, clf.decision_function(Xss))

                candidateIndices = np.where(pr == 1.)[0]

                ii = candidateIndices[np.argmax(rc[candidateIndices])]

                try:
                    chosenThreshold = thres[ii]
                except IndexError:
                    assert ii == thres.shape[0]
                    chosenThreshold = thres[-1] + 0.01


                predictions = (clf.decision_function(Xssval) > chosenThreshold).astype(int)
                predictions[predictions == 0] = -1

                threshold[pInd, fold, i] = chosenThreshold

            innerPrecision[pInd, fold, i], innerRecall[pInd, fold, i], innerF1[pInd, fold, i], _ = precision_recall_fscore_support(yval, predictions, average='binary')
            innerAccuracy[pInd, fold, i] = accuracy_score(yval, predictions)
            innerMCC[pInd, fold, i] = matthews_corrcoef(yval, predictions)
            innerInformedness[pInd, fold, i] = informedness(yval, predictions)


        for i, C in enumerate(paramValues):
            for j, g in enumerate(paramValues):

                clf = SVC(C=C, kernel='rbf', gamma=g)
                clf.fit(Xss, ytrain)

                if thresholdTypeInd == 0:
                    predictions = clf.predict(Xssval)
                else:
                    pr, rc, thres = precision_recall_curve(ytrain, clf.decision_function(Xss))

                    candidateIndices = np.where(pr == 1.)[0]

                    ii = candidateIndices[np.argmax(rc[candidateIndices])]

                    try:
                        chosenThreshold = thres[ii]
                    except IndexError:
                        assert ii == thres.shape[0]
                        chosenThreshold = thres[-1] + 0.01


                    predictions = (clf.decision_function(Xssval) > chosenThreshold).astype(int)
                    predictions[predictions == 0] = -1

                    threshold_RBF[pInd, fold, i, j] = chosenThreshold

                innerPrecision_RBF[pInd, fold, i, j], innerRecall_RBF[pInd, fold, i, j], innerF1_RBF[pInd, fold, i, j], _ = precision_recall_fscore_support(yval, predictions, average='binary')
                innerAccuracy_RBF[pInd, fold, i, j] = accuracy_score(yval, predictions)
                innerMCC_RBF[pInd, fold, i, j] = matthews_corrcoef(yval, predictions)
                innerInformedness_RBF[pInd, fold, i, j] = informedness(yval, predictions)


    if thresholdTypeInd == 0:
        perfLinear = np.mean(innerInformedness[pInd], axis=0)
        perfRBF = np.mean(innerInformedness_RBF[pInd], axis=0)

    else:
        perfLinear = np.mean(innerRecall[pInd], axis=0)
        perfRBF = np.mean(innerRecall_RBF[pInd], axis=0)


    bestLinear = np.max(perfLinear)
    bestRBF = np.max(perfRBF)

    bestModelInd = np.argmax([bestLinear, bestRBF])
    bestModel[pInd] = bestModelInd

    if bestModelInd == 0:
        # lin
        ii = np.argmax(perfLinear)
        clf = LinearSVC(C=paramValues[ii], penalty='l2', loss='hinge', max_iter=500000)
    else:
        # rbf
        ii = np.argmax(perfRBF)
        ind1, ind2 = np.unravel_index(ii, (paramValues.shape[0], paramValues.shape[0]))

        clf = SVC(C=paramValues[ind1], kernel='rbf', gamma=paramValues[ind2])

    ss2 = StandardScaler()
    XtrainOuterCS = ss2.fit_transform(XtrainOuter)
    XtestCS = ss2.transform(Xtest)

    clf.fit(XtrainOuterCS, ytrainOuter)

    if thresholdTypeInd == 0:
        predictions = clf.predict(XtestCS)
    else:
        pr, rc, thres = precision_recall_curve(ytrainOuter, clf.decision_function(XtrainOuterCS))

        candidateIndices = np.where(pr == 1.)[0]

        ii = candidateIndices[np.argmax(rc[candidateIndices])]

        try:
            chosenThreshold = thres[ii]
        except IndexError:
            assert ii == thres.shape[0]
            chosenThreshold = thres[-1] + 0.01

        predictions = (clf.decision_function(XtestCS) > chosenThreshold).astype(int)
        predictions[predictions == 0] = -1


    nMutationsCalled[pInd] = np.sum(predictions == 1)
    precision[pInd], recall[pInd], f1[pInd], _ = precision_recall_fscore_support(ytest, predictions, average='binary')
    acc[pInd] = accuracy_score(ytest, predictions)
    mcc[pInd] = matthews_corrcoef(ytest, predictions)
    informednessYouden[pInd] = informedness(ytest, predictions)


resDict = dict()
resDict['acc'] = acc
resDict['mcc'] = mcc
resDict['rc'] = recall
resDict['f1'] = f1
resDict['pr'] = precision
resDict['inf'] = informednessYouden
resDict['n_muts'] = nMutations
resDict['n_muts_called'] = nMutationsCalled
resDict['test_size'] = testSetSize
resDict['train_size'] = trainingSetSize
resDict['patients'] = allpatientIDs
resDict['best_model'] = bestModel

names = ['0', 'maxRC']

with open('../result_reducedFeats.pkl', 'wb') as f:
    pickle.dump(resDict, f)
