import sys


filename = sys.argv[1]

if len(sys.argv) > 2:
	outfilename = sys.argv[2]
else:
	outfilename = '../..' + filename.split('.')[4] + '_filtered.vcf'

allowedChromosomes = set([str(i) for i in range(23)] + ['X'])

seenGTs = set()

with open(outfilename, 'w') as fw:
    with open(filename) as fr:
        for line in fr:
            if line[0] == '#':
                fw.write(line)

            fields = line.split('\t')

            if fields[0] not in allowedChromosomes:
                continue

            gtConsensus = fields[-2].split(':')[0]
            if gtConsensus == '1/1':
                # homozygous
                continue

            if 'rs' in fields[2]:
                continue

            gtNormal = fields[-3].split(':')[0]
            gtConsensus5 = fields[-1].split(':')[0]

            if '2' in gtNormal or '2' in gtConsensus or '2' in gtConsensus5:
                continue

            if '3' in gtNormal or '3' in gtConsensus or '3' in gtConsensus5:
                continue

            alt = fields[4]
            if ',' in alt:
                continue

            assert len(alt) == 1



            fw.write(line)
