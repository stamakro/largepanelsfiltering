import sys
import pandas as pd
import numpy as np

try:
	ext = sys.argv[1]
except:
	ext = ''


chromosomes = [str(i) for i in range(1,23)] + ['X', 'Y', 'MT']

N = 30

chromosome = []
position = []
refA = []
refC = []
refG = []
refT = []
isRepetitive = []
numAs = []
numCs = []
numGs = []
numTs = []
numNs = []
gcPercent = []
repetitivePercent = []



for chrom in chromosomes:

    with open('../../mutation-contexts' + ext + '/sequences_chr' + chrom + '.txt') as f:
        for line in f:
            region, sequence = line.split()

            L = len(sequence)
            assert L == 2 * N + 1

            sequenceUppercase = sequence.upper()

            chr, bases = region.split(':')
            base5, base3 = bases.split('-')
            base5 = int(base5)
            base3 = int(base3)

            assert base3 - base5 == 2*N+1

            pos = base3 - N

            # if pos == 912359:
            #     print('1111')
            #     sys.exit(0)

            refBase = sequenceUppercase[N]

            rep_position = int(sequence[N].islower())

            rep_percent = sum(map(str.islower, sequence)) / L

            nA = 0
            nC = 0
            nG = 0
            nT = 0
            nN = 0

            for s in sequenceUppercase:
                if s == 'A':
                    nA += 1
                elif s == 'C':
                    nC += 1
                elif s =='G':
                    nG += 1
                elif s =='T':
                    nT += 1
                else:
                    assert s == 'N'
                    nN += 1

            gc = (nG + nC) / L

            chromosome.append(chr)
            position.append(pos)

            if refBase == 'A':
                refA.append(1)
                refC.append(0)
                refG.append(0)
                refT.append(0)
            elif refBase == 'C':
                refA.append(0)
                refC.append(1)
                refG.append(0)
                refT.append(0)
            elif refBase == 'G':
                refA.append(0)
                refC.append(0)
                refG.append(1)
                refT.append(0)
            else:
                assert refBase == 'T'
                refA.append(0)
                refC.append(0)
                refG.append(0)
                refT.append(1)


            isRepetitive.append(rep_position)
            numAs.append(nA)
            numCs.append(nC)
            numGs.append(nG)
            numTs.append(nT)
            numNs.append(nN)
            gcPercent.append(gc)
            repetitivePercent.append(rep_percent)


mydict = {'chr': chromosome,
'pos': position,
'ref_is_A': refA,
'ref_is_C': refC,
'ref_is_G': refG,
'ref_is_T': refT,
'isRepetitive': isRepetitive,
'#As': numAs,
'#Cs': numCs,
'#Gs': numGs,
'#Ts': numTs,
'#Ns': numNs,
'GC%': gcPercent,
'repeatmsk%': repetitivePercent
}
data = pd.DataFrame(data=mydict)

data.to_csv('../../gene-panels/qiaseq-panel-details/genomic-context-muts' + ext + '.csv')
