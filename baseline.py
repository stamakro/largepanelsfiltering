import os
import numpy as np
import pandas as pd
import sys
from sklearn.metrics import precision_recall_curve, average_precision_score, roc_curve, precision_recall_fscore_support, recall_score, accuracy_score, matthews_corrcoef
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
import pickle

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


oncominePath = sys.argv[2]
# first look at all 108 samples
data = pd.read_csv(oncominePath)

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


vcfLocation = sys.argv[3]
allfiles = os.listdir(vcfLocation)

# load feature vectors
with open(sys.argv[4], 'rb'  ) as f:
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
sys.exit(0)
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


tobedelFeats = np.where(np.std(X, axis=0) == 0)[0]
X = np.delete(X, tobedelFeats, axis=1)

featNames = np.delete(featNames, tobedelFeats)


assert P.shape[0] == X.shape[0]


depth = X[:, np.where(featNames=='DP_2')[0][0]]
vaf = X[:, np.where(featNames=='VAF_2')[0][0]]
quality = X[:, np.where(featNames=='QUAL')[0][0]]
known = X[:, np.where(featNames=='COSMIC-variant')[0]]

vafValues = np.array([0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.20])
# depthValues = np.linspace(np.min(depth), np.max(depth), vafValues.shape[0])
# qualityValues = np.linspace(np.min(quality), np.max(quality), vafValues.shape[0])
depthValues = np.array([10, 20, 50, 75, 100, 150, 200, 500, 1000])
qualityValues = np.array([1., 2., 5., 10., 15., 20., 25., 30., 40., 50])

innerPrecision = np.zeros((allpatientIDs.shape[0], vafValues.shape[0], depthValues.shape[0], qualityValues.shape[0]))
innerRecall= np.zeros(innerPrecision.shape)
innerF1 = np.zeros(innerPrecision.shape)
innerMCC = np.zeros(innerPrecision.shape)
innerAccuracy = np.zeros(innerPrecision.shape)
innerInformedness = np.zeros(innerPrecision.shape)


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


for pInd, patient in enumerate(allpatientIDs):

    print('%d/70' % pInd)

    testInd = np.where(P == patient)[0]
    trainInd = np.where(P != patient)[0]

    trainingSetSize[pInd] = trainInd.shape[0]
    testSetSize[pInd] = testInd.shape[0]

    XtrainOuter = X[trainInd]
    ytrainOuter = y[trainInd]

    depthTrain = XtrainOuter[:, np.where(featNames=='DP_2')[0][0]]
    vafTrain = XtrainOuter[:, np.where(featNames=='VAF_2')[0][0]]
    qualityTrain = XtrainOuter[:, np.where(featNames=='QUAL')[0][0]]
    knownTrain = XtrainOuter[:, np.where(featNames=='COSMIC-variant')[0][0]]

    Xtest = X[testInd]
    ytest = y[testInd]

    depthTest = Xtest[:, np.where(featNames=='DP_2')[0][0]]
    vafTest = Xtest[:, np.where(featNames=='VAF_2')[0][0]]
    qualityTest = Xtest[:, np.where(featNames=='QUAL')[0][0]]
    knownTest = Xtest[:, np.where(featNames=='COSMIC-variant')[0][0]]




    nMutations[pInd] = np.sum(ytest == 1)


    for ii, v in enumerate(vafValues):
        for jj, d in enumerate(depthValues):
            for kk, q in enumerate(qualityValues):


                ypred_val = np.logical_and(np.logical_and(depthTrain > d, vafTrain > v), qualityTrain > q)
                ypred_val[knownTrain == 0] = 0

                innerInformedness[pInd, ii, jj, kk] = informedness(ytrainOuter, ypred_val)

    qind, vind, dind = np.unravel_index(np.argmax(innerInformedness[pInd]), innerInformedness.shape[1:])

    qstar = qualityValues[qind]
    vstar = vafValues[vind]
    dstar = depthValues[dind]

    predictions = np.logical_and(np.logical_and(depthTest > dstar, vafTest > vstar), qualityTest > qstar)
    predictions[knownTest == 0] = 0


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


baselineResultPath = sys.argv[5]
with open(baselineResultPath, 'wb') as f:
    pickle.dump(resDict, f)


# plt.show()
