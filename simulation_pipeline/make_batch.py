
import sys
import os

directory = sys.argv[1]
batch_dir = sys.argv[2]
batchfile = batch_dir + '/batchfile.txt'
samples = range(1,7)

f = open(batchfile, 'w')
for sample in samples:
	sample = directory + '/' + str(sample) + '/sim_' + str(sample)
	sample1 = sample + '_1.fq.gz'
	sample2 = sample + '_2.fq.gz'
	f.write(str(sample) + '\t' + sample1 + '\t' + sample2 + '\n')
f.close()


