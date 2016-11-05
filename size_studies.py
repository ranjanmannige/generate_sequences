import numpy,copy,random,os
from scipy.optimize import curve_fit

# Rough code for the study of DNA sequences at different L and Q

NNenergy={}

similaritymatrix={}
for s1 in ["A","T","G","C"]:
	similaritymatrix[s1] = {}
	for s2 in ["A","T","G","C"]:
		if s1 == s2:
			similaritymatrix[s1][s2] = 1
		else:
			similaritymatrix[s1][s2] = 0

#                   Y
#GX/CY     X        A               C               G                T	
NNenergy["GA"] = {"CA":0.17,	"CC":0.81,	"CG":-0.25,	"CT":-1.30}
NNenergy["GC"] = {"CA":0.47,	"CC":0.79,	"CG":-2.24,	"CT":0.62}
NNenergy["GG"] = {"CA":-0.52,	"CC":-1.84,	"CG":-1.11,	"CT":0.08}
NNenergy["GT"] = {"CA":-1.44,	"CC":0.98,	"CG":-0.59,	"CT":0.45}
#CX/GY
NNenergy["CA"] = {"GA":0.43,	"GC":0.75,	"GG":0.03,	"GT":-1.45}
NNenergy["CC"] = {"GA":0.79,	"GC":0.70,	"GG":-1.84, 	"GT":0.62}
NNenergy["CG"] = {"GA":0.11,	"GC":-2.17,	"GG":-0.11,	"GT":-0.47}
NNenergy["CT"] = {"GA":-1.28,	"GC":0.40,	"GG":-0.32,	"GT":-0.12}
#AX/TY
NNenergy["AA"] = {"TA":0.61,	"TC":0.88,	"TG":0.14,	"TT":-1.00}
NNenergy["AC"] = {"TA":0.77,	"TC":1.33,	"TG":-1.44,	"TT":0.64}
NNenergy["AG"] = {"TA":0.02,	"TC":-1.28,	"TG":-0.13,	"TT":0.71}
NNenergy["AT"] = {"TA":-0.88,	"TC":0.73,	"TG":0.07,	"TT":0.69}
#TX/AY
NNenergy["TA"] = {"AA":0.69,	"AC":0.92,	"AG":0.42,	"AT":-0.58}
NNenergy["TC"] = {"AA":1.33,	"AC":1.05,	"AG":-1.30,	"AT":0.97}
NNenergy["TG"] = {"AA":0.74,	"AC":-1.45,	"AG":0.44,	"AT":0.43}
NNenergy["TT"] = {"AA":-1.00,	"AC":0.75,	"AG":0.34,	"AT":0.68}

for seq1 in NNenergy.keys():
	seq1r = seq1[::-1]
	for seq2 in NNenergy[seq1].keys():
		seq2r = seq2[::-1]
		if not seq2r in NNenergy:
			NNenergy[seq2r] = {}
		NNenergy[seq2r][seq1r] = NNenergy[seq1][seq2]

pairwise = {}
pairwise["AT"] = NNenergy["AA"]["TT"]/2.0
pairwise["TA"] = NNenergy["AA"]["TT"]/2.0
pairwise["GC"] = NNenergy["GG"]["CC"]/2.0
pairwise["CG"] = NNenergy["GG"]["CC"]/2.0

# Uncomment if you want to set all native interaction energies to an arbitrary (average) value
#pairwise["AT"] = (NNenergy["AA"]["TT"]+NNenergy["GG"]["CC"])/4.0
#pairwise["TA"] = (NNenergy["AA"]["TT"]+NNenergy["GG"]["CC"])/4.0
#pairwise["GC"] = (NNenergy["AA"]["TT"]+NNenergy["GG"]["CC"])/4.0
#pairwise["CG"] = (NNenergy["AA"]["TT"]+NNenergy["GG"]["CC"])/4.0

singlelinear = {}
singlelinear["A"] = float(NNenergy["AA"]["TT"])/2.0
singlelinear["T"] = float(NNenergy["AA"]["TT"])/2.0
singlelinear["G"] = float(NNenergy["GG"]["CC"])/2.0
singlelinear["C"] = float(NNenergy["GG"]["CC"])/2.0

doublelinear = {}
doublelinear["GA"] = -1.30
doublelinear["GC"] = -2.24
doublelinear["GG"] = -1.84
doublelinear["GT"] = -1.44
doublelinear["CA"] = -1.45
doublelinear["CC"] = -1.84
doublelinear["CG"] = -2.17
doublelinear["CT"] = -1.28
doublelinear["AA"] = -1.00
doublelinear["AC"] = -1.44
doublelinear["AG"] = -1.28
doublelinear["AT"] = -0.88
doublelinear["TA"] = -0.58
doublelinear["TC"] = -1.30
doublelinear["TG"] = -1.45
doublelinear["TT"] = -1.00

Initiation_penalty = 1.96
Terminal_AT_penalty = 0.05
Symmetry_correction = 0.43

def calculate_free_energy(sequence1,sequence2):
	free_energy = 0.0
	for i in range(len(sequence1)-1):
		free_energy += NNenergy[sequence1[i:i+2]][sequence2[i:i+2]]
	#if sequence1[-1] == A
	# remember, no terminal AT punishment is utilized as of now
	free_energy += Initiation_penalty
	if sequence1 == sequence2:
		free_energy += Symmetry_correction
	return free_energy

def calculate_hybrid_free_energy_double(sequence):
	free_energy = 0.0
	for i in range(len(sequence)-1):
		free_energy += doublelinear[sequence[i:i+2]]
	free_energy += Initiation_penalty
	if sequence == sequence[::-1]:
		free_energy += Symmetry_correction
	return free_energy
#

def calculate_hybrid_free_energy_single(sequence):
	free_energy = 0.0
	for i in range(len(sequence)):
		free_energy += singlelinear[sequence[i]]
	free_energy += Initiation_penalty
	if sequence == sequence[::-1]:
		free_energy += Symmetry_correction
	return free_energy
#

def transcribe(s):
	is_hybrid = 1
	transcribed = ""
	for i in range(len(s)):
		if s[i] == "A":
			transcribed += "T"
		elif s[i] == "T":
			transcribed += "A"
		elif s[i] == "G":
			transcribed += "C"
		elif s[i] == "C":
			transcribed += "G"
	return transcribed

native_pairs = [("AA","TT"), ("TT","AA"), ("AT","TA"), ("AT","TA"), ("TA","AT"), ("TA","AT"), ("CA","GT"), ("TG","AC"), ("GT","CA"), ("AC","TG"), 
("CT","GA"), ("AG","TC"), ("GA","CT"), ("TC","AG"), ("CG","GC"), ("CG","GC"), ("GC","CG"), ("GC","CG"), ("GG","CC"), ("CC","GG")]

nonnative_pairs = []
for p1 in sorted(NNenergy.keys()):
	for p2 in sorted(NNenergy[p1].keys()):
		if not (p1,p2) in native_pairs:
			nonnative_pairs.append((p1,p2))


# From Hedges et al., we know:
# 1) Whatever the size of the assembly (Q), the distribution of energies of randomly chosen 
# structures is the same. This distribution is only DNA *length* dependent. So, we can study 
# the distribution of any Q (of large enough size to get decent statistics to arrive at the 
# same results about expected distributions.
if 1: 
	trials = 100
	for Q in [1000]:#,10000,50000]:#100000,1000000]:#,1000000000]:
		print "Q:",Q
		twoQ = 2*Q
		base_directory = "size_studies"
		if not os.path.isdir(base_directory):
			os.makedirs(base_directory)
		n_vs_average = {}
		n_vs_std = {}
		random_sample_average = {}
		random_sample_std = {}
		
		fn = base_directory+"/normalized_l_vs_sd.dat"
		print "----------------------------"
		print "Writing to:",fn
		print "----------------------------"
		l_vs_sd_normalized =  open(fn,"w")
		
		fn = base_directory+"/unnormalized_l_vs_sd.dat"
		print "----------------------------"
		print "Writing to:",fn
		print "----------------------------"
		l_vs_sd_unnormalized =  open(fn,"w")
		
		fn = base_directory+"/normalized_e_minus_2e.dat"
		print "----------------------------"
		print "Writing to:",fn
		print "----------------------------"
		e_minus_twoe_normalized  =  open(fn,"w")
		
		fn = base_directory+"/unnormalized_e_minus_2e.dat"
		print "----------------------------"
		print "Writing to:",fn
		print "----------------------------"
		e_minus_twoe_unnormalized  =  open(fn,"w")
		
		fn = base_directory+"/normalized_e_minus_0.25e.dat"
		print "----------------------------"
		print "Writing to:",fn
		print "----------------------------"
		e_minus_fourthe_normalized  =  open(fn,"w")
		
		fn = base_directory+"/unnormalized_e_minus_0.25e.dat"
		print "----------------------------"
		print "Writing to:",fn
		print "----------------------------"
		e_minus_fourthe_unnormalized  =  open(fn,"w")
		
		fn = base_directory+"/l_vs_energy.dat"
		print "----------------------------"
		print "Writing to:",fn
		print "----------------------------"
		l_vs_energy =  open(fn,"w")
		
		fn = base_directory+"/l_vs_q_limit.dat"
		print "----------------------------"
		print "Writing to:",fn
		print "----------------------------"
		l_vs_q_limit =  open(fn,"w")
		
		for L in range(3,19,1):
			print "\tL:",L
			choices = ["A","T","G","C"]
			
			sequence_sets = []
			energy_sets   = []
			for trial in range(trials):
				sequences = []
				energies  = []
				for instance in range(Q):
					seq = []
					for seqlen in range(L):
						seq.append(random.choice(choices))
					seq = "".join(seq)
					sequences.append(seq)
					energy = calculate_hybrid_free_energy_double(seq)/0.593
					energies.append(energy)
				sequence_sets.append(sequences)
				energy_sets.append(energies)
			
			all_energies_raw = []
			for es in copy.deepcopy(energy_sets):
				all_energies_raw+=es
			average_energy = numpy.average(copy.deepcopy(all_energies_raw))
			std_energy = numpy.std(copy.deepcopy(all_energies_raw))
			
			l_vs_energy.write(str(L)+"\t"+str(average_energy)+"\t"+str(std_energy)+"\n")
			l_vs_energy.flush()
			
			l_vs_q_limit.write(str(L)+"\t"+str(numpy.e**(3.0*abs(average_energy)/2.0))+"\n")
			l_vs_q_limit.flush()
			
			calculate_and_print_params = [[0,2.0,"xmgrace_parameters/histogram2.par"],
				                      [1,2.0,"xmgrace_parameters/histogram.par"],
				                      [0,0.25,"xmgrace_parameters/histogram4.par"],
				                      [1,0.25,"xmgrace_parameters/histogram3.par"]]
			create_agr_files_for = [6,8,13,18]
			for normalize,multiplier,parfn in calculate_and_print_params:
				all_energies = copy.deepcopy(all_energies_raw)
				all_energies_doubled = []
				for i in range(len(all_energies)):
					if normalize == 1:
						all_energies[i] = all_energies[i]/average_energy
					all_energies_doubled.append(copy.deepcopy(all_energies[i])*multiplier)
					
				
				if normalize:
					l_vs_sd_normalized.write(str(L)+"\t"+str(numpy.std(all_energies))+"\n")
					l_vs_sd_normalized.flush()
					
					if multiplier == 2.0:
						e_minus_twoe_normalized.write(str(L)+"\t"+str(numpy.average(all_energies_doubled)-numpy.average(all_energies))+"\n")
						e_minus_twoe_normalized.flush()
					else:
						e_minus_fourthe_normalized.write(str(L)+"\t"+str(numpy.average(all_energies_doubled)-numpy.average(all_energies))+"\n")
						e_minus_fourthe_normalized.flush()
				else:
					l_vs_sd_unnormalized.write(str(L)+"\t"+str(numpy.std(all_energies))+"\n")
					l_vs_sd_unnormalized.flush()
				
					if multiplier == 2.0:
						e_minus_twoe_unnormalized.write(str(L)+"\t"+str(numpy.average(all_energies_doubled)-numpy.average(all_energies))+"\n")
						e_minus_twoe_unnormalized.flush()
					else:
						e_minus_fourthe_unnormalized.write(str(L)+"\t"+str(numpy.average(all_energies_doubled)-numpy.average(all_energies))+"\n")
						e_minus_fourthe_unnormalized.flush()
						
					
				
				postpend = ""
				if float(multiplier) != 2.0:
					postpend = "_one_fourthed"
				
				prepend = "unnormalized"
				if normalize:
					prepend = "normalized"
				
				fn = base_directory+"/l"+str(L)+"_"+prepend+"_energy_histogram"+postpend+".dat"
				print "Writing to:",fn
				f = open(fn,"w")
				
				for values in [all_energies,all_energies_doubled]:
					# Draw a histogram of a ALL energies combined (not just separate sets)
					energymax = numpy.max(values)
					energymin = numpy.min(values)
					a,b = numpy.histogram(values,bins=numpy.arange(energymin,energymax+energymax*0.001,float(energymax-energymin)/100.0))
					asum = numpy.sum(a)
					x = []
					y = []
					for i in range(len(a)):
						v1 = float(b[i]+b[i+1])/2.0
						v2 = float(a[i])/float(asum)
						if v2:
							x.append(v1)
							y.append(v2)
							f.write(str(v1)+" "+str(v2)+"\n")
					f.write("\n\n")
					asum = numpy.sum(a)
					
					x = numpy.array(x)
					y = numpy.array(y)
					n = len(x)                          #the number of data
					mean = sum(x*y)                   #note this correction
					sigma = sum(y*(x-mean)**2)
					
					def gaus(x,a,x0,sigma):
						return a*numpy.exp(-(x-x0)**2/(2*sigma**2))
					
					popt,pcov = curve_fit(gaus,x,y,p0=[0.1,mean,sigma])
					f.write("#L\t=\t"+str(L)+"\n")
					f.write("#A\t=\t"+str(popt[0])+"\n")
					f.write("#X0\t=\t"+str(popt[1])+"\n")
					f.write("#SD\t=\t"+str(popt[2])+"\n")
					f.write("\n\n")
					
					#print sorted(x)
					#print numpy.arange(sorted(x)[0],sorted(x)[-1]+sorted(x)[-1]*0.0001,float(sorted(x)[-1]-sorted(x)[0])/10)
					#exit()
					
					for i in numpy.arange(sorted(x)[0],sorted(x)[-1]+sorted(x)[-1]*0.0001,float(sorted(x)[-1]-sorted(x)[0])/60):
						f.write(str(i)+" "+str(gaus(i,popt[0],popt[1],popt[2]))+"\n")
					f.write("\n\n")
				f.close()
				
				if create_agr_files_for.count(L) and float(multiplier) == 2.0: #!
					agrfn = fn[:-len(".dat")]+".agr"
					psfn = fn[:-len(".dat")]+".ps"
					command = "xmgrace "+fn+" -par "+parfn+" -saveall "+agrfn+" -hardcopy"
					os.system(command)
					command = "xmgrace "+agrfn+" -printfile "+psfn
					os.system(command)
					raw_input("(press enter)")
				
				"""
				for energies in energy_sets:
					energymax = numpy.max(energies)
					energymin = numpy.min(energies)
					a,b = numpy.histogram(energies,bins=numpy.arange(energymin,energymax+energymax*0.001,float(energymax-energymin)/20.0))
					asum = numpy.sum(a)
					for i in range(len(a)):
						v1 = float(b[i]+b[i+1])/2.0
						v2 = float(a[i])/float(asum)
						if v2:
							f.write(str(v1)+" "+str(v2)+"\n")
					f.write("\n\n")
				"""

			
		if 0:
			for aaindex in range(int(N)-1):
				new_sequences = []
				for aastring in sequences:
					for aa in ["A","T","G","C"]:
						new_sequences.append(aastring+aa)
				sequences = copy.deepcopy(new_sequences)
			
			f = open(base_directory+"/n"+str(N)+"_histogram1.dat","w")
			a,b = numpy.histogram(energies,bins=numpy.arange(numpy.min(energies),numpy.max(energies)+0.2,0.5))
			asum = numpy.sum(a)
			x = []
			y = []
			for i in range(len(a)):
				v1 = float(b[i]+b[i+1])/2.0
				v2 = float(a[i])/float(asum)
				if v2:
					x.append(v1)
					y.append(v2)
					f.write(str(v1)+" "+str(v2)+"\n")
			f.write("\n\n")
			asum = numpy.sum(a)
			
			x = numpy.array(x)
			y = numpy.array(y)
			n = len(x)                          #the number of data
			mean = sum(x*y)                   #note this correction
			sigma = sum(y*(x-mean)**2)
			
			def gaus(x,a,x0,sigma):
				return a*numpy.exp(-(x-x0)**2/(2*sigma**2))
			
			popt,pcov = curve_fit(gaus,x,y,p0=[0.1,mean,sigma])
			f.write("#N\t=\t"+str(N)+"\n")
			f.write("#A\t=\t"+str(popt[0])+"\n")
			f.write("#X0\t=\t"+str(popt[1])+"\n")
			f.write("#SD\t=\t"+str(popt[2])+"\n")
			f.write("\n\n")
			for i in numpy.arange(sorted(x)[0],sorted(x)[-1]+0.01,0.1):
				f.write(str(i)+" "+str(gaus(i,popt[0],popt[1],popt[2]))+"\n")
			f.write("\n\n")
			f.close()
			
			energies = energies-average_energy
			
			f = open(base_directory+"/n"+str(N)+"_random_sample.dat","w")
			tmp_random_sample_average = []
			tmp_random_sample_std = []
			if twoQ < len(energies):
				current_number_of_trials = 100
				for i in range(100):
					samples = random.sample(energies,twoQ)
					sample_average = numpy.average(samples)
					sample_std = numpy.std(samples)
					f.write(str(sample_average)+" "+str(sample_std)+"\n")
					tmp_random_sample_average.append(sample_average+average_energy)
					tmp_random_sample_std.append(sample_std)
			else:
				f.write(str(numpy.average(energies+average_energy))+" "+str(numpy.std(energies+average_energy))+"\n\n")
			f.close()
			
			if len(tmp_random_sample_average):
				random_sample_average[N] = [numpy.average(tmp_random_sample_average),numpy.std(tmp_random_sample_average)]
				random_sample_std[N]     = [numpy.average(tmp_random_sample_std)    ,numpy.std(tmp_random_sample_std)    ]
			else:
				random_sample_average[N] = [numpy.average(energies+average_energy),0]
				random_sample_std[N]     = [numpy.std(energies+average_energy)    ,0]
			
			f = open(base_directory+"/n"+str(N)+"_av_vs_std.dat","w")
			for i in range(0,len(energies)-twoQ+1):
				current = energies[i:i+twoQ]
				f.write(str(numpy.average(current))+" "+str(numpy.std(current))+"\n")
			f.write("\n")
			f.close()
			
			f = open(base_directory+"/n"+str(N)+"_histogram2.dat","w")
			a,b = numpy.histogram(energies,bins=numpy.arange(numpy.min(energies),numpy.max(energies)+0.2,0.5))
			asum = numpy.sum(a)
			x = []
			y = []
			for i in range(len(a)):
				v1 = float(b[i]+b[i+1])/2.0
				v2 = float(a[i])/float(asum)
				if v2:
					x.append(v1)
					y.append(v2)
					f.write(str(v1)+" "+str(v2)+"\n")
			f.write("\n\n")
			asum = numpy.sum(a)
			
			x = numpy.array(x)
			y = numpy.array(y)
			n = len(x)                          #the number of data
			mean = sum(x*y)                   #note this correction
			sigma = sum(y*(x-mean)**2)
			
			def gaus(x,a,x0,sigma):
				return a*numpy.exp(-(x-x0)**2/(2*sigma**2))
			
			popt,pcov = curve_fit(gaus,x,y,p0=[0.1,mean,sigma])
			f.write("#N\t=\t"+str(N)+"\n")
			f.write("#A\t=\t"+str(popt[0])+"\n")
			f.write("#X0\t=\t"+str(popt[1])+"\n")
			f.write("#SD\t=\t"+str(popt[2])+"\n")
			f.write("\n\n")
			for i in numpy.arange(sorted(x)[0],sorted(x)[-1]+0.01,0.01):
				f.write(str(i)+" "+str(gaus(i,popt[0],popt[1],popt[2]))+"\n")
			f.write("\n\n")
			f.close()
				
			f = open(base_directory+"/n_vs_average.dat","w")
			for n in n_vs_average.keys():
				f.write(str(n)+" "+str(n_vs_average[n])+"\n")
			f.close()
			
			f = open(base_directory+"/n_vs_std.dat","w")
			for n in n_vs_std.keys():
				f.write(str(n)+" "+str(n_vs_std[n])+"\n")
			f.close()
			
			f = open(base_directory+"/n_vs_average_random.dat","w")
			for n in random_sample_average.keys():
				f.write(str(n)+" "+str(random_sample_average[n][0])+" "+str(random_sample_average[n][1])+"\n")
			f.close()
			
			f = open(base_directory+"/n_vs_std_random.dat","w")
			for n in random_sample_std.keys():
				f.write(str(n)+" "+str(random_sample_std[n][0])+" "+str(random_sample_std[n][1])+"\n")
			f.close()
	
