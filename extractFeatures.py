import numpy as np
import sys
import pandas as pd
import os

def extractFeatures(file, genomicContextDict=None, returnGenomicContext=False, ext=''):

    #
    #  0: is it a COSMIC variant
    #  1: number of different alt alleles
    #  2: variant quality
    #  3: depth
    #  4: VDB
    #  5: SGB
    #  6: RPB
    #  7: MQB
    #  8: MQSB
    #  9: BQB
    # 10: MQ0F
    # 11: AC
    # 12: AN
    # 13: DP4_RF
    # 14: DP4_RR
    # 15: DP4_AF
    # 16: DP4_AR
    # 17: MQ
    # 18: #ref alleleles
    # 19: AD
    # 20: DP
    # 21: VAF
    # 22: #ref alleleles
    # 23: AD
    # 24: DP
    # 25: VAF
    # 26: #ref alleleles
    # 27: AD
    # 28: DP
    # 29: VAF
    # 30: position flagged by repeatmasker
    # 31: #As in 60bp region
    # 32: #Cs in 60bp region
    # 33: #Gs in 60bp region
    # 34: #Ts in 60bp region
    # 35: #Ns in 60bp region
    # 36: GC fraction in 60bp region
    # 37: repeat fraction in 60bp region


    if genomicContextDict is None:

        genomicContext = pd.read_csv('../../gene-panels/qiaseq-panel-details/genomic-context-muts' + ext + '.csv', index_col=0)
        keys = np.array(genomicContext['chr'] + ':' + genomicContext['pos'].astype(str))

        genomicContextDict = dict()
        for i, mut in enumerate(keys):
            assert mut not in genomicContextDict
            genomicContextDict[mut] = np.array(genomicContext.iloc[i, 6:]).astype(float)


    feats = dict()
    featNames = ['COSMIC-variant', '#var-alleles', 'QUAL', 'DP', 'VDB', 'SGB', 'RPB', 'MQB', 'MQSB', 'BQB', 'MQ0F', 'AC', 'AN', 'DP4_RF',
    'DP4_RR', 'DP4_AF', 'DP4_AR', 'MQ', 'REFinGT_1', 'AD_1', 'DP_1', 'VAF_1', 'REFinGT_2', 'AD_2', 'DP_2', 'VAF_2', 'REFinGT_3',
    'AD_3', 'DP_3', 'VAF_3', 'repeat_at_position', 'numAs', 'numCs', 'numGs', 'numTs', 'numNs', 'GC%', 'repeat%']

    assert len(featNames) == 38

    header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'I21-1041-08', 'I21-1041-08-consensus', 'I21-1041-08-consensus-5']



    weirdKeys = set()

    vcf = pd.read_table(file, names=header)

    for i in range(vcf.shape[0]):
        # print(i)
        x = np.zeros((38,), float)

        mut = vcf.iloc[i]
        ref = mut['REF']

        assert len(ref) == 1


        altAlleles = mut['ALT'].split(',')
        assert len(altAlleles) == 1

        if min([len(a) for a in altAlleles]) > 1:
            assert mut['INFO'].split(';')[0] == 'INDEL'
            print('!')
            continue

        alt = altAlleles[0]
        # if len(altAlleles) > 1:
        #     assert max([len(a) for a in altAlleles]) == 1


        chr = str(mut['CHROM'])
        pos = str(mut['POS'])

        numberAllelesCalled = len(altAlleles)

        x[1] = numberAllelesCalled


        ID = mut['ID']

        assert 'rs' not in ID
        if ID[:4] == 'COSV':
            x[0] = 1
        else:
            assert ID == '.'

        x[2] = float(mut['QUAL'])

        info = mut['INFO'].split(';')
        infoDict = {}

        for field in info:
            assert '=' in field
            key, val = field.split('=')

            if key == 'DP4':
                vals = [float(vv) for vv in val.split(',')]
                infoDict['DP4_RF'] = vals[0]
                infoDict['DP4_RR'] = vals[1]
                infoDict['DP4_AF'] = vals[2]
                infoDict['DP4_AR'] = vals[3]
            else:
                try:
                    infoDict[key] = float(val)
                except ValueError:

                    weirdKeys.add(key)
                    infoDict[key] = val

        x[3] = infoDict['DP']

        # AC 12 / VAF 4
        try:
            x[11] = infoDict['AC']
        except ValueError:
            print('!!!')
            x[11] = max([float(ac) for ac in infoDict['AC'].split(',')])


        x[4] = infoDict['VDB']
        x[5] = infoDict['SGB']

        x[6] = infoDict['RPB']
        x[7] = infoDict['MQB']
        x[9] = infoDict['BQB']


        try:
            x[8] = infoDict['MQSB']
        except KeyError:
            x[8] = np.nan
        x[10] = infoDict['MQ0F']

        x[12] = infoDict['AN']
        x[13] = infoDict['DP4_RF']
        x[14] = infoDict['DP4_RR']
        x[15] = infoDict['DP4_AF']
        x[16] = infoDict['DP4_AR']
        x[17] = infoDict['MQ']


        for jj in range(3):
            gt, pl, ad, dp, af = mut[-3 + jj].split(':')
            assert len(ad.split(',')) == 2

            allGTs = set(gt.split('/'))
            if '0' in allGTs:
                if len(allGTs) == 1:
                    x[18 + 4*jj] = 2
                else:
                    x[18 + 4*jj] = 1

            x[19 + 4*jj] = float(ad.split(',')[1])
            x[20 + 4*jj] = float(dp)
            if x[19 + 4*jj] != 0:
                x[21 + 4*jj] = float(af.split(',')[1])
            else:
                x[21 + 4*jj] = 0.

        for jj, kk in enumerate(genomicContextDict['chr' + chr  + ':' +  pos]):
            x[30 + jj] = kk


        kk = chr + '-' + pos + '-' + ref + '-' + alt
        feats[kk] = x

    # print(weirdKeys)
    assert len(featNames) == len(x)

    if returnGenomicContext:
        return feats, featNames, genomicContextDict

    return feats, featNames


try:
    commonPath = sys.argv[1]
    uniquePath = sys.argv[2]
    saveCommon = sys.argv[3]
    saveUnique = sys.argv[4]
    ext = sys.argv[5]
except:
    commonPath = '../../vcf/Qiaseq-common-regions/'
    uniquePath = '../../vcf/Qiaseq-unique-regions/'
    saveCommon = '../../featureVectorsCommon.pkl'
    saveUnique = '../../featureVectorsUnique.pkl'
    ext = ''


for root, _, vcfs in os.walk(commonPath):
    pass

common = dict()
for i, vcf in enumerate(vcfs):
    print('%d/%d' % (i, len(vcfs)))
    if i == 0:
        common[vcf.split('.')[0]], featureNames, gcd = extractFeatures(root + vcf, returnGenomicContext=True, ext=ext)
    else:
        common[vcf.split('.')[0]], _ = extractFeatures(root + vcf, genomicContextDict=gcd, ext=ext)

if len(uniquePath) > 0:
    for root, _, vcfs in os.walk(uniquePath):
        pass

    unique = dict()
    for i, vcf in enumerate(vcfs):
        print('%d/%d' % (i, len(vcfs)))

        unique[vcf.split('.')[0]], _ = extractFeatures(root + vcf, genomicContextDict=gcd, ext=ext)


import pickle
with open(saveCommon, 'wb') as f:
    pickle.dump({'mutDict': common, 'featNames': featureNames}, f)

if len(uniquePath) > 0:
    with open(saveUnique, 'wb') as f:
        pickle.dump({'mutDict': unique, 'featNames': featureNames}, f)
