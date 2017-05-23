import numpy as np

DELIM = "\t"

with open("/Users/lizkennedy/Dropbox/FROM_TARDIS/MESA_in/MESA_betas_NaN.csv", 'r') as infile:
	with open("avg_betarows.txt", 'w') as outfile:
		outfile.write("CpG\tavg_beta\n")
		for line in infile:
			cols = line.strip().split(DELIM)
			vals=np.array([cols[0:]]).astype(np.float)
			meanval = np.nanmean(vals)
			outfile.write(str(meanval))
			outfile.write('\n')
