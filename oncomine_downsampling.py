from Bio.Sequencing.Applications import SamtoolsViewCommandline
import subprocess
import numpy as np
import os
import pandas as pd
import sys

path = '../../gene-panels/'

# count base pairs in the 2 panels
basesOncomine = 0
with open(path + 'Oncomine_Breast_cfDNA_v2.08212017.Designed.bed') as f:
    for i, line in enumerate(f):
        if i > 0:
            fields = line.split('\t')

            basesOncomine += int(fields[2]) - int(fields[1])


basesQiaseq = 0
with open(path + 'DHS-001Z.primers-150bp.bed') as f:
    for i, line in enumerate(f):
        if i > 0:
            fields = line.split('\t')

            basesQiaseq += int(fields[2]) - int(fields[1])


# # delete 2 samples that were run twice
# tobedel = []
# for i, f in enumerate(files):
#     if f[:5] == 'L7063' and '16_R' in f or f[:5] == 'L7041' and '15_R' in f:
#         tobedel.append(i)
#
# for ii in sorted(tobedel)[::-1]:
#     del files[ii]
#
# samples = [f[:5] for f in files]



patientInfo = pd.read_csv('../../Lcode-Oncomine-Qiaseq-per-subject.csv')
oncomine2QiaseqL = pd.Series(patientInfo['Lcode qiaseq'].values, index=patientInfo['Lcode oncomine']).to_dict()
del oncomine2QiaseqL['x']
qiaseqL2oncomine = dict()
for s in oncomine2QiaseqL:
    qiaseqL2oncomine[oncomine2QiaseqL[s]] = s


# start the loop
oncomineBamPath = '/tudelft.net/staff-umbrella/liquidbiopsy/SONIA/bam/Oncomine/'
qiaseqBamPath = '/tudelft.net/staff-umbrella/liquidbiopsy/SONIA/bam/Qiaseq/'
subsampledBamPath = '/tudelft.net/staff-umbrella/liquidbiopsy/SONIA/bam/Subsampled-Oncomine/'
if not os.path.exists(subsampledBamPath):
    os.mkdir(subsampledBamPath)

qiaseqBamFiles = os.listdir(qiaseqBamPath)
qiaseqBams = [q for q in qiaseqBamFiles if q[-4:] == '.bam']

patients  = 0
for bamfile in os.listdir(oncomineBamPath):
    if bamfile[-4:] != '.bam':
        continue


    codeO = bamfile.split('_')[0]

    if codeO == 'L7063' and '16_R' in bamfile or codeO == 'L7041' and '15_R' in bamfile:
        continue

    try:
        codeQ = oncomine2QiaseqL[codeO]
    except:
        continue


    patients += 1
    print(patients, flush=True)

    # if os.path.exists(subsampledBamPath + codeO + '.bam'):
    #     continue

    qiaseqBam = [q for q in qiaseqBams if codeQ in q]
    assert len(qiaseqBam) == 1
    bamfileQiaseq = qiaseqBam[0]
    print(bamfileQiaseq)

    cmd = "bedtools coverage -a ../../gene-panels/DHS-001Z.primers-150bp_GRCh38.chrNames.bed -b %s/%s -d | awk '{ total += $5; count++ } END { print total/count }'" % (qiaseqBamPath, bamfileQiaseq)

    returnValue = subprocess.run([cmd], capture_output=True, shell=True)
    coverageQ = float(returnValue.stdout)

    samtools_view_cmd = SamtoolsViewCommandline(input_file=qiaseqBamPath + bamfileQiaseq)
    samtools_view_cmd.c = True

    returnValue = subprocess.run([str(samtools_view_cmd)], capture_output=True, shell=True)
    readsQ = int(returnValue.stdout)

    readsPerBasePairQ = readsQ / basesQiaseq

    print(bamfile)
    cmd = "bedtools coverage -a ../../gene-panels/Oncomine_Breast_cfDNA_v2.08212017.Designed.bed -b %s/%s -d | awk '{ total += $10; count++ } END { print total/count }'" % (oncomineBamPath, bamfile)

    returnValue = subprocess.run([cmd], capture_output=True, shell=True)
    coverageO = float(returnValue.stdout)

    samtools_view_cmd = SamtoolsViewCommandline(input_file=oncomineBamPath + bamfile)
    samtools_view_cmd.c = True

    returnValue = subprocess.run([str(samtools_view_cmd)], capture_output=True, shell=True)
    readsO = int(returnValue.stdout)

    readsPerBasePairO = readsO / basesOncomine

    # fraction = coveragePerBasePairQ / coveragePerBasePairO
    fraction = coverageQ / coverageO

    print('Total read depth:\t%d\t%d\t%.4f' % (readsO, readsQ, readsQ / readsO))
    print('Read depth per bp:\t%d\t%d\t%.4f' % (readsPerBasePairO, readsPerBasePairQ, readsPerBasePairQ / readsPerBasePairO))
    print('Mean coverage:\t\t%.2f\t%.2f\t%.4f' % (coverageO, coverageQ, fraction))


    if fraction > 1:
        print('Subsampling not required for %s' % codeO)
    else:
        fraction = np.round(fraction, 5)
        # -s INT.FRAC decimal part is the fraction to keep; INT part is the seed, here set to 104
        cmd = 'samtools view -s %3.5f -b %s -o %s' % (104. + fraction, oncomineBamPath + bamfile, subsampledBamPath + codeO + '.bam')
        returnValue = subprocess.run(cmd, capture_output=True, shell=True)

        assert returnValue.returncode == 0

        cmd = "bedtools coverage -a ../../gene-panels/Oncomine_Breast_cfDNA_v2.08212017.Designed.bed -b %s/%s.bam -d | awk '{ total += $10; count++ } END { print total/count }'" % (subsampledBamPath, codeO)

        returnValue = subprocess.run([cmd], capture_output=True, shell=True)
        coverageO2 = float(returnValue.stdout)

        samtools_view_cmd = SamtoolsViewCommandline(input_file=subsampledBamPath + codeO + '.bam')
        samtools_view_cmd.c = True

        returnValue = subprocess.run([str(samtools_view_cmd)], capture_output=True, shell=True)
        readsO2 = int(returnValue.stdout)

        readsPerBasePairO2 = readsO2 / basesOncomine


        print('\nAfter subsampling:')
        print('Total read depth:\t%d' % readsO2)
        print('Read depth per bp:\t%d' % readsPerBasePairO2)
        print('Mean coverage:\t\t%.2f' % coverageO2)


for bamfile in os.listdir(oncomineBamPath + '/Oncomine-files-added-020621/'):
    if bamfile[-4:] != '.bam':
        continue


    codeO = bamfile.split('_')[0]

    if codeO == 'L7063' and '16_R' in bamfile or codeO == 'L7041' and '15_R' in bamfile:
        continue

    try:
        codeQ = oncomine2QiaseqL[codeO]
    except:
        continue


    patients += 1
    print(patients, flush=True)

    # if os.path.exists(subsampledBamPath + codeO + '.bam'):
    #     continue

    qiaseqBam = [q for q in qiaseqBams if codeQ in q]
    assert len(qiaseqBam) == 1
    bamfileQiaseq = qiaseqBam[0]
    print(bamfileQiaseq)

    cmd = "bedtools coverage -a ../../gene-panels/DHS-001Z.primers-150bp_GRCh38.chrNames.bed -b %s/%s -d | awk '{ total += $5; count++ } END { print total/count }'" % (qiaseqBamPath, bamfileQiaseq)

    returnValue = subprocess.run([cmd], capture_output=True, shell=True)
    coverageQ = float(returnValue.stdout)

    samtools_view_cmd = SamtoolsViewCommandline(input_file=qiaseqBamPath + bamfileQiaseq)
    samtools_view_cmd.c = True

    returnValue = subprocess.run([str(samtools_view_cmd)], capture_output=True, shell=True)
    readsQ = int(returnValue.stdout)

    readsPerBasePairQ = readsQ / basesQiaseq

    print(bamfile)
    cmd = "bedtools coverage -a ../../gene-panels/Oncomine_Breast_cfDNA_v2.08212017.Designed.bed -b %s/%s -d | awk '{ total += $10; count++ } END { print total/count }'" % (oncomineBamPath + '/Oncomine-files-added-020621/', bamfile)

    returnValue = subprocess.run([cmd], capture_output=True, shell=True)
    coverageO = float(returnValue.stdout)

    samtools_view_cmd = SamtoolsViewCommandline(input_file=oncomineBamPath + '/Oncomine-files-added-020621/' + bamfile)
    samtools_view_cmd.c = True

    returnValue = subprocess.run([str(samtools_view_cmd)], capture_output=True, shell=True)
    readsO = int(returnValue.stdout)

    readsPerBasePairO = readsO / basesOncomine

    # fraction = coveragePerBasePairQ / coveragePerBasePairO
    fraction = coverageQ / coverageO

    print('Total read depth:\t%d\t%d\t%.4f' % (readsO, readsQ, readsQ / readsO))
    print('Read depth per bp:\t%d\t%d\t%.4f' % (readsPerBasePairO, readsPerBasePairQ, readsPerBasePairQ / readsPerBasePairO))
    print('Mean coverage:\t\t%.2f\t%.2f\t%.4f' % (coverageO, coverageQ, fraction))


    if fraction > 1:
        print('Subsampling not required for %s' % codeO)
    else:
        fraction = np.round(fraction, 5)
        # -s INT.FRAC decimal part is the fraction to keep; INT part is the seed, here set to 104
        cmd = 'samtools view -s %3.5f -b %s -o %s' % (104. + fraction, oncomineBamPath + '/Oncomine-files-added-020621/' + bamfile, subsampledBamPath + codeO + '.bam')
        returnValue = subprocess.run(cmd, capture_output=True, shell=True)

        assert returnValue.returncode == 0

        cmd = "bedtools coverage -a ../../gene-panels/Oncomine_Breast_cfDNA_v2.08212017.Designed.bed -b %s/%s.bam -d | awk '{ total += $10; count++ } END { print total/count }'" % (subsampledBamPath, codeO)

        returnValue = subprocess.run([cmd], capture_output=True, shell=True)
        coverageO2 = float(returnValue.stdout)

        samtools_view_cmd = SamtoolsViewCommandline(input_file=subsampledBamPath + codeO + '.bam')
        samtools_view_cmd.c = True

        returnValue = subprocess.run([str(samtools_view_cmd)], capture_output=True, shell=True)
        readsO2 = int(returnValue.stdout)

        readsPerBasePairO2 = readsO2 / basesOncomine


        print('\nAfter subsampling:')
        print('Total read depth:\t%d' % readsO2)
        print('Read depth per bp:\t%d' % readsPerBasePairO2)
        print('Mean coverage:\t\t%.2f' % coverageO2)
