# CoAct # - Hari Krishna Y & Basil Paul -last update 10/16/2019
# Usage CoAct.py -h #

import os, sys, glob, time, uuid,math, yaml, argparse, itertools
import pandas as pd
import uuid
import concurrent.futures as cf

def multiProcessPattern(passarray):
	try:
		out=passarray[0]+"."+passarray[3]
		patternGenes=extractP(passarray[1],passarray[2],passarray[5])
		overlap=overlapP(passarray[4],patternGenes)
		df = pd.read_csv(overlap, sep='\t', index_col=None)
		os.system("rm "+overlap)
		matches=list(df['PatternID'])
		frequency = {x:matches.count(x) for x in matches}
		df['PatternClass'] = df.apply (lambda row: label_pattern(row,frequency), axis=1)
		df.drop_duplicates(keep=False,inplace=True) 
		df.to_csv(out, sep='\t', index=False)
		return(1)
	except:
		print("Error in processing ",passarray[0])
		return(0)


def label_pattern (row,frequency):
	try:
		return(frequency[row['PatternID']])
	except:
		retunr("NA")

def overlapP(a,b):
	unique_filename = str(uuid.uuid4())
	os.system('bedtools sort -i ' + a + ' > ' + unique_filename+".a.bed")
	os.system( 'bedtools intersect -a ' + b + ' -b ' + unique_filename+".a.bed" + ' -u > '+unique_filename+'.bed')
	tempdf = pd.read_csv(unique_filename+".bed", sep='\t', index_col=None, header=None)
	#tempdf = tempdf.iloc[:, [0, 1, 2, 3, 4,5, 6,7,8,9,10,11]]
	#tempdf.columns = ['Chr1', 'Start1', 'End1', 'Peak1', 'Gene','Strand', 'Distance2TSS','PatternID','Chr2', 'Start2', 'End2', 'Peak2']
	tempdf = tempdf.iloc[:, [0, 1, 2, 3, 4, 5, 6, 7]]
	tempdf.columns = ['Chr1', 'Start1', 'End1', 'Peak1', 'Gene','Strand', 'Distance2TSS','PatternID']
	tempdf.to_csv(unique_filename+'3Peaks_Overlap.bed', sep='\t', index=False)
	os.system("rm "+unique_filename+".bed "+  unique_filename+".a.bed"+" "+b)
	return(unique_filename+'3Peaks_Overlap.bed')

def extractP(df,psd,m):
	unique_filename = str(uuid.uuid4())
	genes=list(set(list(df['Gene'])))
	pdf=pd.DataFrame()
	for gene in genes:
		id=1
		sdf=df[df['Gene']==gene]
		npeaks=sdf.shape[0]
		if npeaks >=m:
			peakst=list(sdf['Start'])
			peaked=list(sdf['End'])
			combinations=itertools.permutations(peakst,m)
			for c in combinations:
				c=(list(c))
				good=1
				for nc in range(0,len(c)-1):
					if (c[nc] < c[nc+1]):
						pass
					else:
						good = 0
				if ((good ==1) and abs(int(c[0])-int(c[-1])) <= psd):
					t=sdf[sdf['Start'].isin(c)].copy()
					t['PatternID']=gene+"_"+str(id)
					pdf=pdf.append(t,ignore_index=True)
					id=id+1
			''' Sequential scan 
			for i in range (0,npeaks-3,1):
				if abs(int(peakst[i+2])-int(peaked[i])) <=psd:
					t=sdf.iloc[i:i+3,:].copy()
					t['PatternID']=gene+"_"+str(id)
					pdf=pdf.append(t,ignore_index=True)
					id=id+1
			'''
	pdf.to_csv(unique_filename+'.bed', sep='\t', index=False, header=None)
	os.system('bedtools sort -i '+unique_filename+'.bed > '+unique_filename+'.bed_T')
	os.system('mv '+unique_filename+'.bed_T '+unique_filename+'.bed')
	return(unique_filename+'.bed')

def add_gene_name(a, b):
	os.system('bedtools sort -i ' + a + ' > ' + a.replace('.bed', '.sorted.bed'))
	cmd = 'bedtools closest -a ' + a.replace('.bed', '.sorted.bed') + ' -b ' + b + ' -id -D b -t first -k 1 > Peak2Gene.bed'
	os.system(cmd)
	tempdf = pd.read_csv('Peak2Gene.bed', sep='\t', index_col=None, header=None)
	tempdf = tempdf.iloc[:, [0, 1, 2, 3, 13, 15,-1]]
	tempdf.columns = ['Chr', 'Start', 'End', 'Peak', 'Gene','Strand', 'Distance']
	os.system('rm Peak2Gene.bed ' + a.replace('.bed', '.sorted.bed'))
	return (tempdf)

def getTSSbed(a):
	fw=open(a.replace(".bed",".TSS.bed"),"w")
	lines=open(a,"r").readlines()
	for line in lines:
		data=line.strip().split("\t")
		if data[5] == "+":
			data[2] = str(int(data[1])+2000)
		if data[5] == "-":
			data[1] = str(int(data[2])-2000)
		fw.write("\t".join(data)+"\n")
	fw.close()
	os.system("bedtools sort -i "+a.replace(".bed",".TSS.bed")+" > "+a.replace(".bed",".TSS.sorted.bed"))
	os.system("rm "+a.replace(".bed",".TSS.bed"))


def main ():
	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	# parse arguments #
	parser = argparse.ArgumentParser(description='CoAct:  -Inferring nuclear receptors from DNA binding patterns. \n')
	parser.add_argument('-p1',help='Peaks file 1 ',required='True',type=str)
	parser.add_argument('-p2',help='Peaks file 2',required='True',type=str)
	parser.add_argument('-r',help='Reference gene annotations in bed format', required='True',type=str)
	parser.add_argument('-d',help='Peak pattern span distance',type=int,default=10000)
	parser.add_argument('-m',help='No of peaks per pattern',type=int,default=3)
	parser.add_argument('-g',help='DEG gene list',type=str)
	parser.add_argument('-t',help='No of threads',required='True',type=int,default=1)
	parser.add_argument('-o',help='Output file',type=str,default="CoAct_Out.txt")
	args=parser.parse_args()
	
	checkfiles=[args.p1]+[args.p2]+[args.r]
	for checkfile in checkfiles:
		if os.path.exists(checkfile):
			pass
		else:
			print("\nError cannot locate "+checkfile+"\n")
			exit()

	if os.path.exists(args.o):
		os.system("rm "+args.o)

	localdate = time.strftime('%a %m/%d/%Y')
	localtime = time.strftime('%H:%M:%S')
	print('# Arguments checked : '+localdate+' at: ' + localtime+' \n')

	a=args.p1; b=args.r; c=args.p2 
	getTSSbed(b)
	df =  add_gene_name(a, b.replace(".bed",".TSS.sorted.bed"))
	os.system("rm "+b.replace(".bed",".TSS.sorted.bed"))
	
	chrlist=["chrX","chrY"]
	passlist=[]
	for i in range(1,23):
		chrlist.append("chr"+str(i))

	for ch in chrlist:
		passlist.append([ch,df[df['Chr']==ch],args.d,args.o,c,args.m])

	result=[]
	with cf.ProcessPoolExecutor(max_workers=int(args.t)) as (executor):
		result = list(executor.map(multiProcessPattern, passlist))

	files=glob.glob("*."+args.o)
	for file in files:
		os.system("cat "+file+" >> "+args.o)
		os.system("rm "+file)

	# Summary #
	fw=open("Summary."+args.o,"w")
	fw.write("CoAct Summary: \n")
	tempdf = pd.read_csv(args.o, sep='\t', index_col=None)
	for n in range(1,args.m+1):
		sdf=tempdf[tempdf['PatternClass']==str(n)]
		genes=list(sdf['Gene'])
		fw.write('No of genes in class '+str(n)+' '+str(len(list(set(genes))))+"\n")
		fw.write(", ".join(list(set(genes)))+'\n')

	if len(args.g) > 0:
		try:
			deg_genes=[]
			lines=open(args.g,"r").readlines()
			for line in lines:
				deg_genes.append(line.strip())
		except:
			print("Error reading DEG list ",args.g)
			exit()
		fw.write("\nCoAct DEG overlap: \n")
		tempdf = pd.read_csv(args.o, sep='\t', index_col=None)
		for n in range(1,args.m+1):
			sdf=tempdf[tempdf['PatternClass']==str(n)]
			genes=list(sdf['Gene'])
			pgenes=list(set(genes).intersection(set(deg_genes)))
			fw.write('No of DEG genes in class '+str(n)+' '+str(len(pgenes))+'\n')
			fw.write(", ".join(pgenes)+"\n")
		fw.close()
		fr=open(args.o,"r").readlines()
		fw=open(args.o,"w")
		for line in fr:
			if "PatternClass" in line:
				fw.write(line.strip()+"\tDEG\n")
			else:
				data=line.split("\t")
				try:
					if data[4] in deg_genes:
						fw.write(line.strip()+"\tYES\n")
					else:
						fw.write(line.strip()+"\tNO\n")
				except:
					pass			






if __name__ == "__main__":
	main()
