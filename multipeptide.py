import sys, os

def read_and_count_pep(n, faafile):# n indicates the number of continuous amino acids    
	res = {}
	seq = ''#the sequence of each protein
	total = 0#the total number of the multi-peptides
	while 1:
		line = faafile.readline()
		if line == '' or line[0] == '>':
			if seq == '':
				continue
			else:
				#construct hash table for the multi-peptides
				for i in range(len(seq)-n+1):
					total += 1
					k = seq[i:i+n]
					if k not in res:
						res[k] = 1
					else:
						res[k] += 1
				seq = ''
				if line == '':
					break
		else:
			l = line.rstrip('\n')
			seq += l
	return res, total

#add_num_to_table(table, myspecies, pep, total)
def	add_num_to_table(table, myspecies, pep, total, specieslist):
	for nmer in pep:
		if nmer not in table:
			table[nmer] = {}
			for bac in specieslist:
				table[nmer][bac] = 0
		table[nmer][myspecies] = pep[nmer] * 1000000.0 / total

def get_all_species(indir):
	tmp_list = []
	for subdir in os.listdir(indir): 
		filepath = os.path.join(indir, subdir)
		if not os.path.isdir(filepath):
			continue
		tmp_list.append(subdir)
	return ( sorted(tmp_list))

def main():
	# single letter variable indicates filehandle
	
	my_py_name = sys.argv[0]

	usage = """ 
			usage:     python3 """ + my_py_name + """ indir AA_num(3,4 or 5) 
			
			"""

	if(len(sys.argv) != 3):
		print(usage)
		sys.exit(0)

	indir = sys.argv[1]
	#numin is the number of amino acid, numout is the most or fewest n peptides in one bacterium
	numin = int(sys.argv[2])
	
	table = {}#the huge table that stores all
	exist = {}#hold all the amino acids combinations
	specieslist = []

	resname = "Kai_Xingang_" + str(numin) + '_peptide_freq_per_million.txt'
	if os.path.exists(resname):
		print (resname + " exists")
		sys.exit(0)

	specieslist = get_all_species(indir)

	for direct in os.listdir(indir):#all the sequences files store here
		filepath = os.path.join(indir, direct)
		if not os.path.isdir(filepath):
			continue
		
		myspecies = direct
		for myfile in os.listdir(filepath):
			if myfile.endswith("faa"):#ignore the other files, only one faa file per directory (the largest one when downloading from NCBI ftp)
				f = open(os.path.join(filepath, myfile))
				#print ('open file...')
				pep, total = read_and_count_pep(numin, f)# adjust the number to do the quad-peptides or penta-peptides
				f.close()
				print ( myspecies + "\t" + str(total) )
				add_num_to_table(table, myspecies, pep, total, specieslist)

	resfile = open( resname,'w')

	resfile.write('name\t')
	resfile.write("\t".join(specieslist))
	resfile.write("\n")

	myspecies_num = len(specieslist)
	for nmer in iter(table):
		resfile.write (nmer + "\t")
		for i in range(0,myspecies_num):
			tmp_num = table[nmer][specieslist[i]] 
			resfile.write( '%.1f' % tmp_num)
			if (i < myspecies_num - 1):
				resfile.write("\t")
			else:
				resfile.write("\n")
	resfile.close()

if __name__ == '__main__':
	main()        