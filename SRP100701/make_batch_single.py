
import sys
import os
import glob
import string

directory = sys.argv[1]
batch_dir = sys.argv[2]
batchfile = batch_dir + '/batchfile.txt'
paths = glob.glob(directory + '/*.fq')
paths.sort()
samples = [os.path.basename(x) for x in paths]
samples = [x.replace('.trim', '') for x in samples]
samples.sort()
print(samples)

f = open(batchfile, 'w')
for i in range(len(samples)):
	f.write(str(samples[i]) + '\t' + str(paths[i]) + '\n')
f.close()


