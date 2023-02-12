import numpy as np
import matplotlib.pyplot as plt
#import myplots
import pickle
import pandas as pd

def drawFigure(nmuts, tps, fps, ax, title='', minorticks=True, horizontal=True):
    if minorticks:
        ax.minorticks_on()

    if horizontal:
        nmuts = nmuts[::-1]
        tps = tps[::-1]
        fps = fps[::-1]

        ax.barh(np.arange(nmuts.shape[0]), -nmuts, color='C0', label='Oncomine calls')
        ax.barh(np.arange(nmuts.shape[0]), tps, color='C1', label='Oncomine & Qiaseq')
        ax.barh(np.arange(nmuts.shape[0]), fps, left=tps, color='C2', label='Only Qiaseq')

        ax.axvline(0, color='k', linewidth=3)

        ax.grid(True, 'both', 'x')
        # ax.set_yticks([0, 1, 2, 3])

        ax.set_ylabel('Patients', fontsize=12)
        ax.set_xlabel('# mutations', fontsize=12)

        ticks = ax.get_xticks()
        newticks = [str(np.abs(tt).astype(int)) for tt in ticks]
        ax.set_xticklabels(newticks)


        plt.setp(ax.get_xticklabels(), fontsize=10)
        plt.setp(ax.get_yticklabels(), fontsize=10)


        ax.set_ylim(-1,70)


    else:
        ax.bar(np.arange(nmuts.shape[0]), -nmuts, color='C0', label='Oncomine calls')
        ax.bar(np.arange(nmuts.shape[0]), tps, color='C1', label='Oncomine & Qiaseq')
        ax.bar(np.arange(nmuts.shape[0]), fps, bottom=tps, color='C2', label='Only Qiaseq')

        ax.axhline(0, color='k', linewidth=3)

        ax.grid(True, 'both', 'y')
        # ax.set_yticks([0, 1, 2, 3])

        ax.set_xlabel('Patients', fontsize=14)
        ax.set_ylabel('# mutations', fontsize=14)

        plt.setp(ax.get_xticklabels(), fontsize=10)
        plt.setp(ax.get_yticklabels(), fontsize=10)


        ax.set_xlim(-1,70)

    plt.legend(prop={'size': 6})

    ax.set_title(title)
    plt.tight_layout()

    return ax


def printPerformance(result, show_medians=False):
    print('Precision: %.3f' % np.mean(result['pr']))
    print('Recall: %.3f' % np.mean(result['rc']))
    print('F1: %.3f' % np.mean(result['f1']))
    print('MCC: %.3f' % np.mean(result['mcc']))
    print('Informedness: %.3f' % np.mean(result['inf']))

    print('\n')
    if show_medians:
        print('medians\n\n')
        print('Precision: %.3f' % np.median(result['pr']))
        print('Recall: %.3f' % np.median(result['rc']))
        print('F1: %.3f' % np.median(result['f1']))
        print('MCC: %.3f' % np.median(result['mcc']))
        print('Informedness: %.3f' % np.median(result['inf']))
        print('\n')

    print('\n')


plt.close('all')


modelName = '_reducedFeats_cosmic'

with open('numberOfMutationsDict.pkl', 'rb') as f:
    nmutsTotalDict = pickle.load(f)


patientInfo = pd.read_csv('../../Lcode-Oncomine-Qiaseq-per-subject.csv')
oncomine2QiaseqL = pd.Series(patientInfo['Lcode qiaseq'].values, index=patientInfo['Lcode oncomine']).to_dict()
qiaseqL2oncomine = dict()
for s in oncomine2QiaseqL:
    qiaseqL2oncomine[oncomine2QiaseqL[s]] = s




with open('result%s.pkl' % modelName, 'rb') as f:
    result = pickle.load(f)

print('SVM %s' % modelName)
printPerformance(result, show_medians=False)


ind = np.argsort(result['n_muts'])

tps = result['n_muts'] * result['rc']
fns = result['n_muts'] - tps
fps = result['n_muts_called'] - tps

patients = result['patients']

tps = result['n_muts'] * result['rc']
fns = result['n_muts'] - tps
fps = result['n_muts_called'] - tps



nmutsTotal = []
for p in patients:
    assert p in qiaseqL2oncomine
    pOnco = qiaseqL2oncomine[p]

    if pOnco in nmutsTotalDict:
        nmutsTotal.append(nmutsTotalDict[pOnco])
    else:
        nmutsTotal.append(0)

nmutsTotal = np.array(nmutsTotal)

ind = np.argsort(nmutsTotal)
nmutsTotal = nmutsTotal[ind]
tps = tps[ind]
fns = fns[ind]
fps = fps[ind]

tst = pd.DataFrame({'total': nmutsTotal, 'tp': tps, 'fp': fps})
tst.sort_values(['total', 'tp', 'fp'], ascending=True, inplace=True)


fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(1,3,3)
ax = drawFigure(tst['total'], tst['tp'], tst['fp'], ax, title='(C) SVM-based filtering', minorticks=False)

#fig.savefig('../../figures/result_model_filtering.png', dpi=600)
#fig.savefig('../../figures/result_model_filtering.eps', dpi=1200)


# std qiaseq pipeline - cosmic filter
dd = pd.read_csv('../../calls_cosmic.csv', index_col=0)
dd['call-binary'] = dd.apply(lambda x: 1 if x['call'] == 1 else 0, axis=1)



ncalls = dd.groupby('patient')['call-binary'].count().to_dict()
tpcalls = dd.groupby('patient')['call-binary'].sum().to_dict()

ncalls2 = []
ncalls2Total = []
tpcalls2 = []

for p in patients:
    pOnco = qiaseqL2oncomine[p]
    if pOnco in nmutsTotalDict:
        ncalls2Total.append(nmutsTotalDict[pOnco])
    else:
        ncalls2Total.append(0)
    ncalls2.append(ncalls[p])
    tpcalls2.append(tpcalls[p])


ncalls2 = np.array(ncalls2)
tpcalls2 = np.array(tpcalls2)
ncalls2Total = np.array(ncalls2Total)


fpcalls2 = ncalls2 - tpcalls2

precision = tpcalls2 / (tpcalls2 + fpcalls2)
f1 = 2 * precision / (precision + 1)

print('Precision: %.2f' % np.mean(precision))
print('F1: %.2f' % np.mean(f1))

print('\n\n')

tst = pd.DataFrame({'total': ncalls2Total, 'tp': tpcalls2, 'fp': fpcalls2})
tst.sort_values(['total', 'tp', 'fp'], ascending=True, inplace=True)


#fig = plt.figure()
ax = fig.add_subplot(1,3,1)
ax = drawFigure(tst['total'], tst['tp'], tst['fp'], ax, title='(A) No filtering - COSMIC variants')

# fig.savefig('../../figures/result_std_pipeline_cosmic.png', dpi=600)
# fig.savefig('../../figures/result_std_pipeline_cosmic.eps', dpi=1200)

####### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########
with open('result_baseline.pkl', 'rb') as f:
    result2 = pickle.load(f)

print('baseline')
printPerformance(result2, show_medians=False)
tps = result2['n_muts'] * result['rc']
fns = result2['n_muts'] - tps
fps = result2['n_muts_called'] - tps

tps = tps[ind]
fns = fns[ind]
fps = fps[ind]

patients = result['patients']


tst = pd.DataFrame({'total': nmutsTotal, 'tp': tps, 'fp': fps})
tst.sort_values(['total', 'tp', 'fp'], ascending=True, inplace=True)



# fig = plt.figure()
ax = fig.add_subplot(1,3,2)
ax = drawFigure(tst['total'], tst['tp'], tst['fp'], ax, title='(B) manual filtering')


#fig.tight_layout()
plt.subplots_adjust(left=0.07, right=0.979, top=0.925, bottom=0.125, hspace=0.2, wspace=0.272)
# fig.savefig('../../figures/result_basic_filtering.png', dpi=600)
# fig.savefig('../../figures/result_basic_filtering.eps', dpi=1200)
fig.savefig('../../figures/fig2.png', dpi=300)
fig.savefig('../../figures/fig2.svg', dpi=300)

#
plt.show()
