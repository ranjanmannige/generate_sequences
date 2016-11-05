<pre>
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
</pre>
