import sys

try:
    ext = sys.argv[1]
except:
    ext = ''

N = 30

chromosomes = [str(i) for i in range(1,23)] + ['X', 'Y', 'MT']
writeFiles = [open('../../mutation-contexts' + ext + '/locations_chr' + i + '.bed', 'w') for i in chromosomes]


weird = set()
with open('../../all-qiaseq-calls' + ext + '.txt', 'r') as fr:
    for line in fr:
        chr, pos = line.split()

        try:
            ind = chromosomes.index(chr)
        except:
            weird.add(chr)
            continue

        pos = int(pos)

        if chr == 'MT':
            chr = 'M'

        writeFiles[ind].write('chr%s\t%d\t%s\n' % (chr, pos-N-1, pos+N))


# print(weird)
