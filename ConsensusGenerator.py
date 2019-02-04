import exrex
import os

newList = list(exrex.generate('[ACGTRYN]{6}'))

print(len(newList))

with open(os.path.expanduser('~/Desktop/consensus2.txt'), 'w+') as file:
    for x in newList:
        file.write('%s,' % x)