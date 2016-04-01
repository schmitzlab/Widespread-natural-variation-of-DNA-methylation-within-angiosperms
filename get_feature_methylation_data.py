#!/usr/local/bin/python2.7
## get the methylation levels for features (genes, transposons, etc). Prior to use, Use Bedtools intersect and sort on allc file with appropriate bedfile to obtain input data:
## bedtools intersect -a <INPUT ALLC> -b <INPUT BED> -wb | cut -f4,5,6,7,11 | sort -s -k5,5 > x.tsv


import sys
import os
import getopt

def main(argv):
	input = argv[0] ##sorted intersecting data from allc 
	output = argv[1] ##output file
	
	data = open(input)
	out = open(output, "w")
	
	name = "none"
	
	tCG = 0
	tmCG = 0
	CG = 0
	mCG = 0
	tCHG = 0
	tmCHG = 0
	CHG = 0
	mCHG = 0
	tCHH = 0
	tmCHH = 0
	CHH = 0
	mCHH = 0
	tC = 0
	tmC = 0
	C = 0
	mC = 0
	
	for line in data:
		per_row = line.split('\t')
		
		if per_row[4] == name:
		
			tC = tC + float(str(per_row[2]))
			tmC = tmC + float(str(per_row[1]))
			C = C + 1
			if int(per_row[3]) == 1:
				mC = mC + 1
			
			if str(per_row[0]) == "CGA" or str(per_row[0]) == "CGT" or str(per_row[0]) == "CGG" or str(per_row[0]) == "CGC" or str(per_row[0]) == "CGN" or str(per_row[0]) == "CG":
 	 			tCG = tCG + float(str(per_row[2]))
 	 			tmCG = tmCG + float(str(per_row[1]))
 	 			CG = CG + 1
 	 			if int(per_row[3]) == 1:
 	 				mCG = mCG + 1
 	 		elif str(per_row[0]) == "CAG" or str(per_row[0]) == "CTG" or str(per_row[0]) == "CCG" or str(per_row[0]) == "CHG":
				tCHG = tCHG + float(str(per_row[2]))
				tmCHG = tmCHG + float(str(per_row[1]))
				CHG = CHG + 1
				if int(per_row[3]) == 1:
 	 				mCHG = mCHG + 1
 	 		elif str(per_row[0]) == "CAA" or str(per_row[0]) == "CAT" or str(per_row[0]) == "CAC" or str(per_row[0]) == "CTA" or str(per_row[0]) == "CTT" or str(per_row[0]) == "CTC" or str(per_row[0]) == "CCA" or str(per_row[0]) == "CCT" or str(per_row[0]) == "CCC" or str(per_row[0]) == "CHH":
				tCHH = tCHH + float(str(per_row[2]))
				tmCHH = tmCHH + float(str(per_row[1]))
				CHH = CHH + 1
				if int(per_row[3]) == 1:
 	 				mCHH = mCHH + 1
				
		else:
			
			if name == "none":
			
				name = per_row[4]
			
				out.write("Gene" + "\t" + "CG" + "\t" + "mCG" + "\t" + "total_CG" + "\t" + "methyl_CG" + "\t" + "weighted_CG" + "\t" + "CHG" + "\t" + "mCHG" + "\t" + "total_CHG" + "\t" + "methyl_CHG" + "\t" + "weighted_CHG" + "\t" + "CHH" + "\t" + "mCHH" + "\t" + "total_CHH" + "\t" + "methyl_CHH" + "\t" + "weighted_CHH" + "\t" + "C" + "\t" + "mC" + "\t" + "total_C" + "\t" + "methyl_C" + "\t" + "weghted_C" + "\n")
			
			else:
			
				if tCG == 0:
					wCG = "NA"
				else:
					wCG = float(tmCG/tCG)
				if tCHG == 0:
					wCHG = "NA"
				else:
					wCHG = float(tmCHG/tCHG)
				if tCHH == 0:
					wCHH = "NA"
				else:
					wCHH = float(tmCHH/tCHH)
				if tC == 0:
					wC = "NA"
				else:
					wC = float(tmC/tC)
			
				out.write( str(name.strip('\n')) + "\t" + str(CG) + "\t" + str(mCG) + "\t" + str(tCG)[:-2] + "\t" + str(tmCG)[:-2] + "\t" + str(wCG) + "\t" + str(CHG) + "\t" + str(mCHG) + "\t" + str(tCHG)[:-2] + "\t" + str(tmCHG)[:-2] + "\t" + str(wCHG) + "\t" + str(CHH) + "\t" + str(mCHH) + "\t" + str(tCHH)[:-2] + "\t" + str(tmCHH)[:-2] + "\t" + str(wCHH) + "\t" + str(C) + "\t" + str(mC) + "\t" + str(tC)[:-2] + "\t" + str(tmC)[:-2] + "\t" + str(wC) + "\n")
			
				name = per_row[4]
				
				tCG = 0
				tmCG = 0
				CG = 0
				mCG = 0
				tCHG = 0
				tmCHG = 0
				CHG = 0
				mCHG = 0
				tCHH = 0
				tmCHH = 0
				CHH = 0
				mCHH = 0
				tC = 0
				tmC = 0
				C = 0
				mC = 0
			
				tC = tC + float(str(per_row[2]))
				tmC = tmC + float(str(per_row[1])) 
				C = C + 1
				if int(per_row[3]) == 1:
					mC = mC + 1
			
				if str(per_row[0]) == "CGA" or str(per_row[0]) == "CGT" or str(per_row[0]) == "CGG" or str(per_row[0]) == "CGC" or str(per_row[0]) == "CGN" or str(per_row[0]) == "CG":
 	 				tCG = tCG + float(str(per_row[2]))
 	 				tmCG = tmCG + float(str(per_row[1]))
 	 				CG = CG + 1
 	 				if int(per_row[3]) == 1:
 	 					mCG = mCG + 1
 	 			elif str(per_row[0]) == "CAG" or str(per_row[0]) == "CTG" or str(per_row[0]) == "CCG" or str(per_row[0]) == "CHG":
					tCHG = tCHG + float(str(per_row[2]))
					tmCHG = tmCHG + float(str(per_row[1]))
					CHG = CHG + 1
					if int(per_row[3]) == 1:
 	 					mCHG = mCHG + 1
 	 			elif str(per_row[0]) == "CAA" or str(per_row[0]) == "CAT" or str(per_row[0]) == "CAC" or str(per_row[0]) == "CTA" or str(per_row[0]) == "CTT" or str(per_row[0]) == "CTC" or str(per_row[0]) == "CCA" or str(per_row[0]) == "CCT" or str(per_row[0]) == "CCC" or str(per_row[0]) == "CHH":
					tCHH = tCHH + float(str(per_row[2]))
					tmCHH = tmCHH + float(str(per_row[1]))
					CHH = CHH + 1
					if int(per_row[3]) == 1:
 	 					mCHH = mCHH + 1
	
	if tCG == 0:
		wCG = "NA"
	else:
		wCG = float(tmCG/tCG)
	if tCHG == 0:
		wCHG = "NA"
	else:
		wCHG = float(tmCHG/tCHG)
	if tCHH == 0:
		wCHH = "NA"
	else:
		wCHH = float(tmCHH/tCHH)
	if tC == 0:
		wC = "NA"
	else:
		wC = float(tmC/tC)
					
	out.write( str(name.strip('\n')) + "\t" + str(CG) + "\t" + str(mCG) + "\t" + str(tCG)[:-2] + "\t" + str(tmCG)[:-2] + "\t" + str(wCG) + "\t" + str(CHG) + "\t" + str(mCHG) + "\t" + str(tCHG)[:-2] + "\t" + str(tmCHG)[:-2] + "\t" + str(wCHG) + "\t" + str(CHH) + "\t" + str(mCHH) + "\t" + str(tCHH)[:-2] + "\t" + str(tmCHH)[:-2] + "\t" + str(wCHH) + "\t" + str(C) + "\t" + str(mC) + "\t" + str(tC)[:-2] + "\t" + str(tmC)[:-2] + "\t" + str(wC) + "\n")
	
	data.close()
	out.close()  
		
if __name__ == "__main__":
   main(sys.argv[1:])