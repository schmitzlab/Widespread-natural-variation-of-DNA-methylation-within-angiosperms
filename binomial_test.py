#!/usr/local/bin/python2.7
##Perform binomial test on each gene


import sys
import os
import getopt
from scipy import stats

def main(argv):
	input = argv[0] ##input methylation data. This should be output file generated from the "get_feature_methylation_data.py" script
	output = argv[1] ##output file, duh
	
	data = open(input)
	out = open(output, "w")
	
	out.write("Locus	total_CG	mCG	pval_CG	total_CHG	mCHG	pval_CHG	total_CHH	mCHH	pval_CHH" + "\n")
	
	CG = 0
	mCG = 0
	CHG = 0
	mCHG = 0
	CHH = 0
	mCHH = 0
	
	next(data)
	for line in data:
		per_row = line.split('\t')
		
		CG = int(per_row[1])
		mCG = int(per_row[2])
		CHG = int(per_row[6])
		mCHG = int(per_row[7])
		CHH = int(per_row[11])
		mCHH = int(per_row[12])
		
		pCG = stats.binom.sf(mCG-1,CG,0.236862727155) ##These values were determined from the coding sequence methylation levels of all 34 species in the study
		pCHG = stats.binom.sf(mCHG-1,CHG,0.0554156350293)
		pCHH = stats.binom.sf(mCHH-1,CHH,0.0115659429912)
		
		out.write(str(per_row[0]) + "\t")
		out.write(str(CG) + "\t")
		out.write(str(mCG) + "\t")
		out.write(str(pCG) + "\t")
		out.write(str(CHG) + "\t")
		out.write(str(mCHG) + "\t")	
		out.write(str(pCHG) + "\t")
		out.write(str(CHH) + "\t")
		out.write(str(mCHH) + "\t")	
		out.write(str(pCHH) + "\n")	
		
	data.close()
	out.close()
	
if __name__ == "__main__":
   main(sys.argv[1:])
		
		