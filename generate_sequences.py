#!/usr/bin/env python 
__author__ = 'Ranjan V. Mannige'

# To do:
#python generate_sequences.py -L 23 -Q 10000000 -k 4 -d 1.0 -o 0
#python generate_sequences.py -L 23 -Q 100000000 -k 4 -d 1.0 -o 0
#python generate_sequences.py -L 23 -Q 1000000000 -k 4 -d 1.0 -o 0

verbose_helpfile = """
 _____                           _         _____                                           
|  __ \                         | |       /  ___|                                          
| |  \/ ___ _ __   ___ _ __ __ _| |_ ___  \ `--.  ___  __ _ _   _  ___ _ __   ___ ___  ___ 
| | __ / _ \ '_ \ / _ \ '__/ _` | __/ _ \  `--. \/ _ \/ _` | | | |/ _ \ '_ \ / __/ _ \/ __|
| |_\ \  __/ | | |  __/ | | (_| | ||  __/ /\__/ /  __/ (_| | |_| |  __/ | | | (_|  __/\__ \
 \____/\___|_| |_|\___|_|  \__,_|\__\___| \____/ \___|\__, |\__,_|\___|_| |_|\___\___||___/
                                                         | |                               
                                                         |_| v 1.0.beta
# --------------
# QUICK SYNOPSIS
# --------------
# SETUP: 
# Assume that you want Q number of distinct DNA bricks to that are 
# expected to form a specific assembly. Assume that on average each DNA
# brick interacts with k other neighbors. Assume that the "sticky" regions 
# of the DNA are of length L.
# --------------
# ALL RESULTS ARE OUTPUTTED IN THE DIRECTORY ./sequences/L/
# For any successful run invluving L=<L>,Q=<Q>,and k=<k>, the following files 
# are created:
# 1) sequences[_tripletsallowed][_tight]_L<L>_Q<Q>_k<k>.seq.gz
# 2) sequences[_tripletsallowed][_tight]_L<L>_Q<Q>_k<k>.ene.gz
# 3) sequences[_tripletsallowed][_tight]_L<L>_Q<Q>_k<k>.his 
# "_tight" is used in the filename if -d or --distribution_thickness is provided. 
# "_tripletsallowed" is used in the filename if -a or --allow_triplets is used.
# The .seq and .ene files respectively list 2<Q> number of acceptable 
# sequences and  hybridization energies. The .his file displays the distribution 
# of energies for a random set of sequences, and for a set of sequences indicated
# by the user options. E.g., the second histogram would look exactly like the 
# first if there is not use of -d of --distribution_thickness. But the second 
# curve would range between -1/2*std and 1/2*std if -d 1 is provided (see below).
# --------------
# WHAT THIS CODE DOES:
# This script generates a "master set" of sequences to be used as sticky 
# patches within a dna brick assembly of size Q. As the faithful formation of
# these structures are strongly dependent on the extent to which the interacion
# energies are distributed (Hedges, Mannige, Whitelam, Soft Matter, 2014), We 
# have introduced a user input "distribution_thickness", which controls how tightly 
# distributed the energies of the set of sequences are (if distribution_thickness>0).
# ----------------------
# MORE DETAILED SYNOPSIS:
# -
# "Q" is the number of DNA bricks/subunits.
# 
# "k" is the number of sticky patches per block (coordination).
# 
# "L" is the length (in nucleotides) of DNA to be used as sticky patches.
# 
# "disallow_triplets", when 1, will not use any sequences with "AAA", "TTT", "GGG", "CCC" in them.
# 
# "distribution_thickness" (in units of KbT or energy) is the allowed range of sequence energy.
# if set to 0 or negative, then sequences from the entire sequence distribution will 
# be selected; then "offset_from_average" (below) will effectively be 0.
# 
# "offset_from_average" (in standard deviation) is a little more involved:
# This feature allows you to choose the average energies that are offset (in units of standard 
#  deviations) from the average energy. So, "offset_from_average=-1" will mean that we expect to
#  select sequences close to A - S<offset_from_average> in KbT, where A and S are the average and 
#  standard deviations of the energies in KbT (the allowed energies will then be 
#  A - S<offset_from_average>  +/- "distribution_thickness"/2; see above. A warning will be 
#  listed if there are not enough sequences available within that range in a combinatorial
#  sense).
# -------------------
# TESTED IN PYTHON 2.7. Please email ranjanmannige@gmail.com for any bugs, or suggestions as to 
# how to port to PYTHON >3.
# -------------------
"""

import numpy # for a number of math-related functions (min, max, average, std, histogram)
import copy,random,os,sys
import argparse # For parsing command-line arguments
from scipy.optimize import curve_fit # for fitting the energy histogram to a gaussian
import gzip # Write the *.seq and *.ene files in compressed formats (*.seq.gz, *.ene.gz).
            # This reduces the size of each file by ~1/10th.

#====================================================================
# SETTING DEFAULTS AND PARSING USER INPUTS
#====================================================================
# DEFAULT VALUES
# --------------
Q  = 100000 # number of distinct subunits in an assembly
k  = 4      # number of interactions that each subunit makes (excepting edges)
L  = 13     # The length of the sticky portion of each interacting component within the brick
disallow_triplets      = 1 # AAA, TTT, GGG, CCC are not allowed within a sticky sequence
distribution_thickness = 1.0 # (in KT) if 0, then no limits on the distribution 
                             # of interaction energies will be made.
# More control:
offset_from_average    = 0   # (in standard deviation)

# -------------------
# PARSING USER INPUTS
# -------------------
# OPENING A PARSER
parser = argparse.ArgumentParser(description='Program that produces DNA brick sequences of specific properties.', 
usage='\n--------\nDEFAULT: python generate_sequences.py -L '+str(L)+' -Q '+str(Q)+' -k '+str(k)+' -d '+str(distribution_thickness)+' -o '+str(offset_from_average)+'\n--------')

# ADD ARGUMENTS
parser.add_argument('--verbose_help',dest='verbose_help',
		    action='store_const', const=1,required=False,default=0,metavar='',
		    help='More details about the script and input variables')
parser.add_argument('-L','--length', dest='L',type=int,required=False,default=L,metavar='',
		    help='Length of each DNA sticky patch (default: %(default)s)')
parser.add_argument('-Q','--subunits',dest='Q',type=int,required=False,default=Q,metavar='',
		    help='Number of DNA bricks in assembly (default: %(default)s)')
parser.add_argument('-k','--coordination', dest='k',type=int,required=False,default=k,metavar='',
		    help='Number of sticky patches per block (default: %(default)s)')
parser.add_argument('-a','--allow_triplets',dest='disallow_triplets',
		    action='store_const',const=0,default=disallow_triplets,required=False,metavar='',
		    help='Exclude triplets (e.g., "AAA","TTT",...) in sequences (default: %(default)s)')
parser.add_argument('-d','--distribution_thickness',dest='distribution_thickness',
		    type=float,required=False,default=distribution_thickness,metavar='',
		    help='If >0, then tighten sequence distribution (in KbT; default: %(default)s)')
parser.add_argument('-o','--offset_from_average',dest='offset_from_average',
		    type=float,required=False,default=offset_from_average,metavar='',
		    help='Useful if distribution_thickness=0; offset from the average energy (in st.d. in KbT; default: %(default)s)')
#parser.add_argument('-b','--base_directory',dest='base_dir',type=str,required=False,default="./sequences/",metavar='',
#		    help='Base output directory for saving output files (default: ./sequences/<L>/)')
parser.add_argument('--run_with_defaults',dest='defaults', action='store_const',const=1,default=0,required=False,metavar='',
		    help='Run with defaults')

# ARRAY FOR ALL ARGUMENTS PASSED TO SCRIPT
args = parser.parse_args()

# ASSIGN ARGS TO VARIABLES
L                      = args.L
Q                      = args.Q
k                      = args.k
distribution_thickness = args.distribution_thickness
disallow_triplets      = args.disallow_triplets
offset_from_average    = args.offset_from_average
base_dir               = "./sequences/"+str(L)+"/"

# IF NO PARAMETERS WERE ASSIGNED, EVOKE HELP
if len(sys.argv) < 2:
	print
	parser.print_help()
	print
	exit()

# INVOKING VERBOSE HELP IF NEEDED
if args.verbose_help:
	print verbose_helpfile
	exit()

# CHECKING IF THE BASE DIRECTORY EXISTS
if not os.path.isdir(base_dir):
	if os.path.exists(base_dir):
		print "WARNING: '"+base_dir+"' does not appear to be a directory. Exiting."
		exit()
	else:
		print "'"+base_dir+"' does not exist. Creating."
		# MAKING IT IF IT DOES NOT EXIST
		print os.makedirs(base_dir)

#====================================================================
# SETTING UP THE MACHINERY TO CALCULATE HYBRIDIZATION FREE ENERGIES
# 
# Nucleotide-nucleotide hybridization free energies obtained from:
# John SantaLucia Jr.; Donald Hicks (2004). "The thermodynamics of DNA structural 
# motifs". Annual Review of Biophysics and Biomolecular Structure 33:415--440.
# 
#====================================================================

# BUILDING THE NUCLEOTIDE PAIR-PAIR ENERGY MATRIX AS A DICTIONARY
choices = ["A","T","G","C"] # THE DNA ALPHABET (DON'T CHANGE UNLESS YOU HAVE NEW NUCLEOTIDES 
                            # AND DETAILED INFORMATION ABOUT INTERACTION ENERGIES!)
NNenergy={} # THE ENERGIES OF INTERACTION BETWEEN PAIRS OF PAIRS (SantaLucia & Hicks, 2004).
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

# FILLING IN THE REDUNDANT VALUES
for seq1 in NNenergy.keys():
	seq1r = seq1[::-1]
	for seq2 in NNenergy[seq1].keys():
		seq2r = seq2[::-1]
		if not seq2r in NNenergy:
			NNenergy[seq2r] = {}
		NNenergy[seq2r][seq1r] = NNenergy[seq1][seq2]

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

Initiation_penalty  = 1.96
Terminal_AT_penalty = 0.05
Symmetry_correction = 0.43

range_l_minus_one = range(L-1)
def calculate_hybrid_free_energy_double(sequence):
	free_energy = 0.0
	for i in range_l_minus_one:
		free_energy += doublelinear[sequence[i:i+2]]
	free_energy += Initiation_penalty
	if sequence == sequence[::-1]:
		free_energy += Symmetry_correction
	return free_energy

#====================================================================
# FIRST WE TRY TO SAMPLE SEQUENCES ASSOCIATED WITH LENGTH L TO SEE 
# WHAT THE ENERGY DISTRIBUTION IS.
#====================================================================

print "# Obtaining a starting collection of sequences and energies"
all_sequences_to_energies_raw = {} # This will contain a dictionary of a sequence to its energy
all_sequences = [] # a list of all relevant sequences encountered in the search
all_energies  = [] # a corresponding list of interaction energies to each sequence's complement
initial_sequence_count = 10000
count = 0
non_triplet = 0.0
triplet     = 0.0
#alphapetchoices = choices #["A","T","G","C"]

triplet_frequenies = []
largest_count = copy.deepcopy(count)
for i in range(10):
	triplet = 0.0
	non_triplet = 0.0
	count = 0
	sequences = []
	energies = []
	while count < initial_sequence_count:
		seq = ""
		for seqlen in range(L):
			seq+=random.choice(choices)
			#Alternative method (not sure which is more efficient):
			#seq+=alphapetchoices[random.randrange(4)]
		
		# EXCLUDING REPEATS
		use_sequence = 1
		if disallow_triplets:
			if "AAA" in seq or "TTT" in seq or "GGG" in seq or "CCC" in seq:
				use_sequence = 0
				triplet += 1
			else:
				non_triplet +=1
		
		if use_sequence:
			# To show progress, uncomment the following:
			"""
			if count % 1000 == 0:
				print "#\t",count, "of",initial_sequence_count,"("+str(round(100.0*float(count)/initial_sequence_count,1))+"%)"
			"""
			energy = calculate_hybrid_free_energy_double(seq)/0.593 # KbT ~ 0.593 at standard temperature
			all_sequences_to_energies_raw[seq] = energy
			
			sequences.append(seq)
			energies.append(energy)
			
			count=len(sequences)
	all_sequences += copy.deepcopy(sequences)
	all_energies  += copy.deepcopy(energies)
	if disallow_triplets:
		triplet_frequenies.append(float(triplet)/(triplet+non_triplet))
#
triplet_frequency = 0.0
if disallow_triplets:
	triplet_frequency = numpy.average(triplet_frequenies)

print "# Frequency of triplets =",triplet_frequency

# NOW WE KNOW THE FOLLOWING PIECES OF INFORMATION:
# 1) triplet frequency (0, if disallow_triplets==1),
# 2) total number of sequences (4**L)
# So, the total number of ALLOWED sequences is:

allowed_sequences = (4**L)*(1-triplet_frequency)


print "# The number of possible sequences = (4^L)*(1-triplet_frequency) =",allowed_sequences

ave = numpy.average(all_energies)
std = numpy.std(all_energies)

# FIRST WE GET A RELATIONSHIP BETEWEN ENERGIES AND PROBABILITY OF 
# RANDOMLY FINDING SUCH AN ENERGY IN A SEQUENCE
energymax = numpy.average(all_energies) + 1.05*(numpy.max(all_energies) - numpy.min(all_energies))/2.0
energymin = numpy.average(all_energies) - 1.05*(numpy.max(all_energies) - numpy.min(all_energies))/2.0

# Calculating the histogram
energystep = 1.0
if distribution_thickness > 0: 
	# then we are interested in a smaller range than the whole (natural) energy range of 
	# squences of length L (sans triplet containing sequences)
	energystep = distribution_thickness
	a,b = numpy.histogram(all_energies,bins=numpy.arange(energymin,energymax+energystep/2.0,energystep))
else:
	a,b = numpy.histogram(all_energies,bins=numpy.arange(energymin,energymax+float(energymax-energymin)/100.0,float(energymax-energymin)/50.0))

asum = numpy.sum(a) # for normalization purposes
unconstrained_energy_distribution_x = []
unconstrained_energy_distribution_y = []
ymax = -1.0
for i in range(len(a)):
	x = float(b[i]+b[i+1])/2.0
	y = float(a[i])/float(asum)
	if y > ymax:
		ymax = y
	if y:
		unconstrained_energy_distribution_x.append(x)
		unconstrained_energy_distribution_y.append(y)

# ================================================
# NOW WE CREATE A GAUSSIAN FIT TO THE DISTRIBUTION
def gaus(x,a,x0,sigma):
	return a*numpy.exp(-(x-x0)**2/(2*sigma**2))
#
x = numpy.array(unconstrained_energy_distribution_x)
y = numpy.array(unconstrained_energy_distribution_y)
n = len(x)
mean = sum(x*y)
sigma = sum(y*(x-mean)**2)
popt,pcov = curve_fit(gaus,x,y,p0=[0.1,mean,sigma])
# ================================================

ave = numpy.average(all_energies)
std = numpy.std(all_energies)

target_energy = ( ave + offset_from_average * std ) #- float(energystep)/2.0

target_energy_frequency = gaus(target_energy,popt[0],popt[1],popt[2])

expected_abundance = allowed_sequences
if distribution_thickness<=0:
	expected_abundance = int(allowed_sequences*target_energy_frequency)

print "# Selected energy = ",target_energy 
print "# Expected number of sequences at selected energy = ",int(expected_abundance),"sequences"
print "# (distribution_thickness = "+str(distribution_thickness)

# We require the following number of sequences (abundance) to create the
required_abundance = int(float(Q)*float(k)/2.0)
if required_abundance > expected_abundance:
	print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	print "  WARNING: NOT ENOUGH SEQUENCES FOR THE LISTED PARAMETERS!"
	print "  required_abundance ("+str(required_abundance)+") > expected_abundance ("+str(expected_abundance)+")"
	print "# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	exit(0)

#====================================================================

final_sequences_to_energies = {}

emin = target_energy - float(distribution_thickness)/2
emax = target_energy + float(distribution_thickness)/2
if distribution_thickness <= 0.0:
	emin = -1000000000000000.0 # lazy person's method of ensuring that emin and emax
	emax = 1000000000000000.0  #  will not be the pre-set values.

print "# Pre-polulating the sequence list..."
for seq,ene in zip(all_sequences,all_energies):
	if emin <= ene and ene <= emax:
		final_sequences_to_energies[seq] = ene
		# Just printing it out at regular intervals (10000)
		if len(final_sequences_to_energies) % 10000 == 0:
			print "#\t",len(final_sequences_to_energies), "of",required_abundance,"("+str(round(100.0*float(len(final_sequences_to_energies))/required_abundance,2))+"%)"
#

print "# Further populating the sequence collection..."
count = len(final_sequences_to_energies)
time_count = 0
old_count = 0
inefficient_warning_count = 0
last_time_the_warning_was_raised = 0 # 
while count < required_abundance:
	time_count+=1
	seq = ""
	for seqlen in range(L):
		seq+=random.choice(choices)
	
	#if not seq in final_sequences:
	# EXCLUDING REPEATS
	use_sequence = 1
	if "AAA" in seq or "TTT" in seq or "GGG" in seq or "CCC" in seq:
		if disallow_triplets:
			use_sequence = 0
		#triplet += 1
		#pass
	
	if use_sequence:
		ene = calculate_hybrid_free_energy_double(seq)/0.593
		if emin <= ene and ene <= emax:
			#
			final_sequences_to_energies[seq] = ene
			
			count = len(final_sequences_to_energies)
			if count % 10000 == 0:
				if old_count == count:
					print "WARNING: not an efficient search (Q might be too high for L)."
					last_time_the_warning_was_raised += 1
					if last_time_the_warning_was_raised == 1:
						inefficient_warning_count += 1
						last_time_the_warning_was_raised = 0
					else:
						inefficient_warning_count = 0
					'''
					if inefficient_warning_count==5:
						print "OK, it seems like you might be asking too much from a poor DNA sequence of length L. I quit!"
						exit()
					'''
				print "#\t",count, "of",required_abundance,"("+str(round(100.0*float(count)/required_abundance,4))+"%)"
			old_count = count
	
print "#\t ..."
final_sequences = []
final_energies   = []
for seq,ene in final_sequences_to_energies.items():
	final_sequences.append(seq)
	final_energies.append(ene)
print "#\t...done"


#====================================================================
# FINALLY WRITING TO FILES!
#
#====================================================================

triplet_indicator = "_tripletsallowed"
if disallow_triplets:
	triplet_indicator = ""
tightness_indicator = "_tight"
if distribution_thickness <= 0:
	tightness_indicator = ""
#
base_fn = base_dir+"/sequences"+tightness_indicator+triplet_indicator+"_L"+str(L)+"_Q"+str(Q)+"_k"+str(k)+"_ave"+str(round(numpy.average(final_energies),2))+"_std"+str(round(numpy.std(final_energies),2))
if 1: # simple name
	base_fn = base_dir+"/sequences"+tightness_indicator+"_L"+str(L)+"_Q"+str(Q)+"_k"+str(k)

if 1:
	# =====================================================================
	# Cmpressing the seqfn and enefn files, as they are quite big otherwise
	# =====================================================================
	
	seqfn = base_fn+".seq.gz"
	print "Writing to:",seqfn
	f = gzip.open(seqfn,"w")
	f.write("# ==============================================================================\n")
	f.write("# SEQUENCE SET GENERATED BY GENERATE_SEQUENCE.PY USING THE FOLLOWING PARAMETERS:\n")
	f.write("# Q = "+str(Q)+"\n")
	f.write("# k = "+str(k)+"\n")
	f.write("# L = "+str(L)+"\n")
	f.write("# disallow_triplets      = "+str(disallow_triplets)+"\n")
	f.write("# distribution_thickness = "+str(distribution_thickness)+" (in KbT)\n")
	f.write("# offset_from_average    = "+str(offset_from_average)+"\t(in standard deviation of energy/KbT )\n")
	f.write("# ------------------\n")
	f.write("# ACTUAL PROPERTIES OF THE SEQUENCE SET:\n")
	f.write("# AVERAGE ENERGY:     "+str(numpy.average(final_energies))+"/KbT\n")
	f.write("# STANDARD DEVIATION: "+str(numpy.std(final_energies))+"/KbT\n")
	f.write("# ------------------\n")
	f.write("# PROPERTIES OF A RANDOMLY SELECTED SEQUENCE SET:\n")
	f.write("# AVERAGE ENERGY:     "+str(numpy.average(all_energies))+"/KbT\n")
	f.write("# STANDARD DEVIATION: "+str(numpy.std(all_energies))+"/KbT\n")
	f.write("# ------------------\n")
	
	for s in final_sequences:
		f.write(s+"\n")
	f.close()
	
	enefn = base_fn+".ene.gz"
	print "Writing to:",enefn
	f = gzip.open(enefn,"w")
	f.write("# ==============================================================================\n")
	f.write("# SEQUENCE SET GENERATED BY GENERATE_SEQUENCE.PY USING THE FOLLOWING PARAMETERS:\n")
	f.write("# Q = "+str(Q)+"\n")
	f.write("# k = "+str(k)+"\n")
	f.write("# L = "+str(L)+"\n")
	f.write("# disallow_triplets      = "+str(disallow_triplets)+"\n")
	f.write("# distribution_thickness = "+str(distribution_thickness)+" (in KbT)\n")
	f.write("# offset_from_average    = "+str(offset_from_average)+"\t(in standard deviation of energy/KbT )\n")
	f.write("# ------------------\n")
	f.write("# ACTUAL PROPERTIES OF THE SEQUENCE SET:\n")
	f.write("# AVERAGE ENERGY:     "+str(numpy.average(final_energies))+"/KbT\n")
	f.write("# STANDARD DEVIATION: "+str(numpy.std(final_energies))+"/KbT\n")
	f.write("# ------------------\n")
	f.write("# PROPERTIES OF A RANDOMLY SELECTED SEQUENCE SET:\n")
	f.write("# AVERAGE ENERGY:     "+str(numpy.average(all_energies))+"/KbT\n")
	f.write("# STANDARD DEVIATION: "+str(numpy.std(all_energies))+"/KbT\n")
	f.write("# ------------------\n")
	
	for e in final_energies:
		f.write(str(e)+"\n")
	f.close()
	
	hisfn = base_fn+".his"
	print "Writing to:",hisfn
	f = open(hisfn,"w")
	
	f.write("# THIS FILE CONTAINS TWO DISTRIBUTIONS:\n")
	f.write("# 2. THE GAUSSIAN DISTRIBUTION OF ENERGIES DISPLAYED BY SEQUENCES\n")
	f.write("#    OF LENGTH '"+str(L)+"' AVAILABLE IN '"+seqfn+"'\n")
	f.write("# 1. A GUASSIAN DISTRIBUTION OF UNCONSTRAINED ENERGIES OF SEQUENCE LENGTH '"+str(L)+"'\n")
	f.write("\n")
	f.write("\n")
	
	f.write("# =======================================================================\n")
	f.write("# DISTRIBUTION #1:\n")
	f.write("# ----------------\n")
	f.write("# ENERGY DISTRIBUTION OF SEQUENCES IN '"+seqfn+"'\n")
	f.write("# GENERATED BY 'generate_sequence.py' USING THE FOLLOWING PARAMETERS:\n")	
	f.write("# Q = "+str(Q)+"\n")
	f.write("# k = "+str(k)+"\n")
	f.write("# L = "+str(L)+"\n")
	f.write("# disallow_triplets      = "+str(disallow_triplets)+"\n")
	f.write("# distribution_thickness = "+str(distribution_thickness)+" (in KbT)\n")
	f.write("# offset_from_average    = "+str(offset_from_average)+"\t(in standard deviation of energy/KbT )\n")
	f.write("# ------------------\n")
	f.write("# ACTUAL PROPERTIES OF THE SEQUENCE SET:\n")
	f.write("# AVERAGE ENERGY:     "+str(numpy.average(final_energies))+"/KbT\n")
	f.write("# STANDARD DEVIATION: "+str(numpy.std(final_energies))+"/KbT\n")
	f.write("# ------------------\n")
	f.write("# PROPERTIES OF A RANDOMLY SELECTED SEQUENCE SET:\n")
	f.write("# AVERAGE ENERGY:     "+str(numpy.average(all_energies))+"/KbT\n")
	f.write("# STANDARD DEVIATION: "+str(numpy.std(all_energies))+"/KbT\n")
	f.write("# ------------------\n")
	
	values = final_energies
	energymax = numpy.max(values)
	energymin = numpy.min(values)
	a,b = numpy.histogram(values,bins=numpy.arange(energymin,energymax+float(energymax-energymin)/100.0,float(energymax-energymin)/50.0))
	asum = numpy.sum(a)
	x = []
	y = []
	ymax = -1.0
	for i in range(len(a)):
		v1 = float(b[i]+b[i+1])/2.0
		v2 = float(a[i])/float(asum)
		if v2 > ymax:
			ymax = v2
		if v2:
			x.append(v1)
			y.append(v2)
			f.write(str(v1)+" "+str(v2)+"\n")
	
	f.write("\n\n# =======================================================================\n")
	f.write("# DISTRIBUTION #2:\n")
	f.write("# ----------------\n")
	f.write("# DISTRIBUTION OF A RANDOMLY SELECTED SET OF SEQUENCES FIT TO A GAUSSIAN\n")
	f.write("# THE GAUSSIAN STATISTICS ARE:\n")
	f.write("#L\t=\t"+str(L)+"\n")
	f.write("#A\t=\t"+str(popt[0])+"\n")
	f.write("#X0\t=\t"+str(popt[1])+"\n")
	f.write("#SD\t=\t"+str(popt[2])+"\n")
	f.write("# =======================================================================\n")
	sx = sorted(unconstrained_energy_distribution_x)
	for i in numpy.arange(sx[0],sx[-1]+float(sx[-1]-sx[0])/200,float(sx[-1]-sx[0])/100):
		f.write(str(i)+" "+str(gaus(i,popt[0],popt[1],popt[2]))+"\n")
	
	f.close()

print "DONE!"
print
