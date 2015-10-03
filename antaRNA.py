import numpy
import sys
import random
import subprocess
import re
import decimal
import math
import os
import shutil
import time
import types
import uuid
import argparse
from argparse import RawTextHelpFormatter


##############################
# ARGPARSE TYPES AND FUNCTIONS
##############################

def parseGC(string):
	"""
		Argparse Type Function
		parses a string of "FLOAT:INT-INT", where ":INT-INT" is set to be optional.
		Returns a single value or a triple (FLAOT, INT, INT).
	"""
	m = re.match('^\s*(\d+\.\d+):?(\d+)?\s*?-?\s*?(\d+)?\s*?$', string)
	if m is None:
		print "'" + string + "' is not a valid input for tGC Expected forms like '0-5'."
		exit(1)
	tGC = float(m.group(1))
	if m.group(2) and m.group(3):
		start = int(m.group(2))
		end = int(m.group(3))
		return (tGC, start, end)
	else:
		return (tGC)
    
def convert_arg_line_to_args(arg_line):
    for arg in arg_line.split():
        if not arg.strip():
            continue
        yield arg
        
#####################
# SUPPORT FUNCTIONS
#####################

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False
    
    
def print2file(f, i, m):
	"""
		print content i to file f in mode m
	"""
	line = str(i)
	if m == "a":
		call = "echo \"" +  line + "\" >> " + f
	elif m == "w":
		call = "echo \"" +  line + "\" > " + f
	os.system(call)

def getUsedTime(start_time):
	"""
		Return the used time between -start time- and now.
	"""
	end_time = time.time()
	return end_time - start_time

def substr(x, string, subst):
	"""
		Classical substring function
	"""
	s1 = string[:x-1]
  
	s2 = string[x-1:x]
	s3 = string[x:]
  #s2 = s[x+len(string)-x-1:]
  
	return s1 + subst + s3
	

####################################################
# STRUCTURE AND SEQUENCE INTEGRITY CHECK FUNCTIONS
####################################################

def checkSequenceConstraint(SC):
	"""
		Checks the Sequence constraint for illegal nucleotide characters
	"""
	for c in SC:
		if c not in "ACGURYSWKMBDHVNacgu": 
			print "\tIllegal Character in the constraint sequence!"
			print "\tPlease use the IUPAC nomenclature for defining nucleotides in the constraint sequence!"
			print "\tA   	Adenine"
			print "\tC   	Cytosine"
			print "\tG   	Guanine"
			print "\tT/U 	Thymine/Uracil"
			print "\tR 	A or G"
			print "\tY 	C or T/U"
			print "\tS 	G or C"
			print "\tW 	A or T/U"
			print "\tK 	G or T/U"
			print "\tM 	A or C"
			print "\tB 	C or G or T/U"
			print "\tD 	A or G or T/U"
			print "\tH 	A or C or T/U"
			print "\tV 	A or C or G"
			print "\tN	any base"
			print "\tOr their lowercase soft constraint variants of ACGU!"
			exit(1)
  
def checkSimilarLength(s, SC):
	"""
		Compares sequence and structure constraint length
	"""
	if len(s) != len(SC):
		print "Sequence and Structure-Constraint provide two different lengths..."
		print "Structure Constraint length : " + str(len(s))
		print "Sequence Constraint length  : " + str(len(SC))
		exit(1)
  
def isStructure(s):
	"""
		Checks if the structure constraint only contains "(", ")", and "." and legal fuzzy structure constraint characters.
	"""
	returnvalue = 1
	for a in range(0,len(s)):
		if s[a] not in  ".()[]{}<>":
			if s[a] not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz":
				returnvalue = 0
	return returnvalue   

def isBalanced(s):
	"""
		Check if the structure s is of a balanced nature
	"""
	
	balance = 1
	for bracket in ["()", "[]", "{}", "<>"]:
		counter = 0
		for a in xrange(len(s)):
			if s[a] in bracket[0]:
				counter += 1
			elif s[a] in bracket[1]:
				counter -= 1
		if counter != 0:
			balance = 0
	return balance

def fulfillsHairpinRule(s):
	"""
		CHECKING FOR THE 3 nt LOOP INTERSPACE
			for all kind of basepairs, even wihtin the pdeudoknots 
	"""
	
	fulfillsRules = 1
	for bracket in ["()", "[]", "{}", "<>"]:
		last_opening_char = 0
		check = 0
		for a in xrange(len(s)):
			if s[a] == bracket[0]:
				last_opening_char = a
				check = 1
			elif s[a] == bracket[1] and check == 1:
				check = 0
				if a - last_opening_char < 4:
					return 0
	return 1
	
def isValidStructure(s):
	"""
		Checks, if the structure s is a valid structure
	"""
	Structure = isStructure(s)
	Balanced = isBalanced(s)
	HairpinRule = fulfillsHairpinRule(s)
	if not Structure == 1 or not Balanced == 1 or not HairpinRule == 1:
		print "Structure Integrity Issue!"
		exit(1)
		
def isStructureCompatible(lp1, lp2 ,bp): 
	"""
		Checks, if the region within lp1 and lp2 is structurally balanced
	"""
	x = lp1 + 1
	while (x < lp2):
		if (bp[x] <= lp1 or bp[x] > lp2):
			return False
		if x == bp[x]:
			x += 1
		else:
			x = bp[x] + 1
	return x == lp2
		
def checkConstaintCompatibility(args):
	"""
		Checks if the constraints are compatible to each other
	"""
	for id1 in args.BPstack:  # key = (constraint , (pos, constraint)))
		constr1 = args.BPstack[id1][0]
		id2 = args.BPstack[id1][1][0]
		constr2 = args.BPstack[id1][1][1]    
		if id1 != id2 and not isCompatible(constr1, constr2, args.IUPAC_compatibles):
			print "Contraint Compatibility Issue:"
			print "Nucleotide constraint " + str(constr1) + " at position " + str(id1) + " is not compatible with nucleotide constraint " + str(constr2) + " at position " + str(id2) + "\n"
			exit(1)

def transform(seq):
	"""
		Transforms "U" to "T" for the processing is done on DNA alphabet
	"""
	S = ""
	for s in seq:
		if s == "T":
			S += "U"
		elif s == "t":
			S += "u"
		else:
			S += s
	return S

def complementBase(c):
	"""
		Returns the complement RNA character of c (without GU base pairs)
	"""
	retChar = ""
	if c == "A" :
		retChar = "U"
	elif c == "U":
		retChar = "A"
	elif c == "C":
		retChar = "G"
	elif c == "G":
		retChar = "C"
	return retChar  


############################################
# IUPAC LOADINGS AND NUCLEOTIDE MANAGEMENT
############################################

def loadIUPACcompatibilities(args):
	"""
		Generating a hash containing all compatibilities of all IUPAC RNA NUCLEOTIDES
	"""
	compatible = {}
	for nuc1 in args.IUPAC: # ITERATING OVER THE DIFFERENT GROUPS OF IUPAC CODE
		sn1 = list(args.IUPAC[nuc1])
		for nuc2 in args.IUPAC: # ITERATING OVER THE DIFFERENT GROUPS OF IUPAC CODE
			sn2 = list(args.IUPAC[nuc2])
			compatib = 0
			for c1 in sn1: # ITERATING OVER THE SINGLE NUCLEOTIDES WITHIN THE RESPECTIVE IUPAC CODE:
				for c2 in sn2: # ITERATING OVER THE SINGLE NUCLEOTIDES WITHIN THE RESPECTIVE IUPAC CODE:
					# CHECKING THEIR COMPATIBILITY
					if args.noGUBasePair == False:
						if (c1 == "A" and c2 == "U") or (c1 == "U" and c2 == "A") or (c1 == "C" and c2 == "G") or (c1 == "G" and c2 == "C") or (c1 == "G" and c2 == "U") or (c1 == "U" and c2 == "G"): 
							compatib = 1
					else:
						if (c1 == "A" and c2 == "U") or (c1 == "U" and c2 == "A") or (c1 == "C" and c2 == "G") or (c1 == "G" and c2 == "C"): 
							compatib = 1
			compatible[nuc1 + "_" + nuc2] = compatib # SAVING THE RESPECTIVE GROUP COMPATIBILITY, REVERSE SAVING IS NOT REQUIRED, SINCE ITERATING OVER ALL AGAINST ALL
	return compatible
	
def isCompatibleToSet(c1, c2, IUPAC_compatibles):
	"""
		Checks compatibility of c1 wihtin c2
	"""
	compatible = True
	for setmember in c2:
		if isCompatible(c1, setmember, IUPAC_compatibles) == False: 
			return False
	return compatible
	
	
def isCompatible(c1, c2, IUPAC_compatibles):
	"""
		Checks compatibility between character c1 and c2
	"""
	if IUPAC_compatibles[c1.upper() + "_" + c2.upper()] == 1:
		return True
	else:
		return False


#######################################
# STRUCTURE AND DICTIONARY MANAGEMENT
#######################################

def getLP(BPSTACK):  
	"""
	Retreives valid lonley base pairs from a base pair stack
	"""
	#20 ('N', (>BLOCK<, 'N'))

	# getting single base pairs
	stack = {}
	LP = {}
	if type(BPSTACK[random.choice(BPSTACK.keys())]) == types.TupleType:
		for i in BPSTACK.keys():
			#if str(BPSTACK[i][1][0]) not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
			stack[i] = int(BPSTACK[i][1][0])
		#print i , BPSTACK[i][1][0]
	else: 
		for i in BPSTACK.keys():
			#if str(BPSTACK[i]) not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
			stack[i] = BPSTACK[i]

	# removing redundant base pair indices 
	for i in stack.keys():
		if i >= stack[i]:
			del stack[i]

	# actual checking for single lonley base pairs
	for i in stack.keys():
		if not (i-1 in stack and stack[i-1] == stack[i] + 1) and not (i+1 in stack and stack[i+1] == stack[i] - 1): 
			LP[i] = stack[i]

	##actual removal of 2er lonley base pairs
	for i in stack.keys():
		if not (i-1 in stack and stack[i-1] == stack[i] + 1) and  (i+1 in stack and stack[i+1] == stack[i] - 1) and not (i+2 in stack and stack[i+2] == stack[i] - 2): 
			LP[i] = stack[i]
			LP[i+1] = stack[i+1]
	return LP
  
def getBPStack(args):
	"""
		Returns a dictionary of the corresponding basepairs of the structure s and the sequence constraint seq.
	"""
	tmp_stack = {"()":[], "{}":[], "[]":[], "<>":[]}
	BPstack = {}
	for i in xrange(len(args.Cstr)):
		
    # REGULAR SECONDARY STRUCTURE DETECTION
		if args.Cstr[i] in "(){}[]<>":

			no = 0
			### opening
			if args.Cstr[i] in "([{<":
				if args.Cstr[i] == "(":
					tmp_stack["()"].append((i, args.Cseq[i]))
				elif args.Cstr[i] == "[":
					tmp_stack["[]"].append((i, args.Cseq[i]))
				elif args.Cstr[i] == "{":
					tmp_stack["{}"].append((i, args.Cseq[i]))
				elif args.Cstr[i] == "<":
					tmp_stack["<>"].append((i, args.Cseq[i]))

			#closing
			elif args.Cstr[i] in ")]}>":
				if args.Cstr[i]  == ")":
					no, constr = tmp_stack["()"].pop() 
				elif args.Cstr[i]  == "]":
					no, constr = tmp_stack["[]"].pop() 
				elif args.Cstr[i]  == "}":
					no, constr = tmp_stack["{}"].pop() 
				elif args.Cstr[i]  == ">":
					no, constr = tmp_stack["<>"].pop() 
				BPstack[no] = (constr, (i, args.Cseq[i])) 
				BPstack[i] = (args.Cseq[i] ,(no, constr)) 

		elif args.Cstr[i] == ".": 
			BPstack[i] = (args.Cseq[i], (i, args.Cseq[i])) 
		elif args.Cstr[i] in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz":
			BPstack[i] = (args.Cseq[i], (i, args.Cseq[i])) 
	 
	return (BPstack, getLP(BPstack))
	

def getbpStack(structure): 
	"""
		Returns a dictionary of the corresponding basepairs of the structure s and the sequence constraint seq.
	"""
	tmp_stack = {"()":[], "{}":[], "[]":[], "<>":[]}
	bpstack = {}

	for i in xrange(len(structure)):
		if structure[i] in "(){}[]<>":

			no = 0
			### opening
			if structure[i] in "([{<":
				if structure[i] == "(":
					tmp_stack["()"].append(i)
				elif structure[i] == "[":
					tmp_stack["[]"].append(i)
				elif structure[i] == "{":
					tmp_stack["{}"].append(i)
				elif structure[i] == "<":
					tmp_stack["<>"].append(i)

			#closing
			elif structure[i] in ")]}>":
				if structure[i] == ")":
					no = tmp_stack["()"].pop() 
				elif structure[i] == "]":
					no = tmp_stack["[]"].pop() 
				elif structure[i] == "}": 
					no = tmp_stack["{}"].pop() 
				elif structure[i] == ">": 
					no = tmp_stack["<>"].pop() 
				bpstack[no] = i # save basepair in the format {opening base id (opening seq constr,(closing base id, closing seq constr))}
				bpstack[i] = no # save basepair in the format {closing base id (closing seq constr,(opening base id, opening seq constr))}

		elif structure[i] == ".": # no structural constaint given: produce entry, which references itself as a base pair partner....
			bpstack[i] = i
	
		elif structure[i] in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz":
			bpstack[i] = i
    
	return (bpstack, getLP(bpstack))

############################################
# TERRAIN GRAPH MANAGEMENT
############################################

def maprange( a, b, s):
  """
   Mapping function
  """
  (a1, a2), (b1, b2) = a, b
  return  b1 + ((s - a1) * (b2 - b1) / (a2 - a1))

def getConstraint(TE, args):
	"""
	Dependend on the situation in the constraint an the respective path section, setting wether a specific constraint can be given or not (for that path section)
	"""
	# TE :: transition element / path section under dispute
	# id1 :: id of the position of the caharacter to which the transition is leading to
	# id2 :: id of the position of the character, which is listed in the BPinformation, it can be id1 as well, when no bp is present
	# val :: BPstack information of the specific position
	# constr1 :: constraining character of pos id1
	# constr2 :: constraining character of pos id2

	id1 = int(TE.split(".")[0])
	val = args.BPstack[id1] # check out the value of the destination character in the basepair/constraint stack
	constr1 = val[0] # getting the constraint character of position id1
	id2 = int(val[1][0]) # getting position id2
	constr2 = val[1][1] # getting the sequence constraint for position id2
	targetNucleotide = TE.split(".")[1][-1:] # where the edge is leading to
	
	lowerCase = constr1.islower()
	
	if lowerCase:
		return 1
	else:
		c1 = set(args.IUPAC[constr1.upper()]) # getting all explicit symbols of c1
		c2 = set(args.IUPAC_reverseComplements[constr2.upper()]) # getting the reverse complement explicit symbols of c2
		
		if targetNucleotide in c1:
			if id1 == id2:
				return 1
			else:
				if targetNucleotide in c2:
					return 1
				else:
					return 0
		else:
			return 0

def initTerrain(args): 
	"""
		Initialization of the terrain graph  in equidistant mode.
		Vertices are modeled implicitly.
	"""
	nt = ["A","C","G","U"] 
	nt2 = ["AA","AC","AG","AU","CA","CC","CG","CU","GA","GC","GG","GU","UA","UC","UG","UU"] # Allowed dinucleotides
	args.terrain = {}
	pathlength = 1
	pheromone = 1
	for p in xrange(len(args.Cstr)):
		if p == 0:
			for i in nt:
				args.terrain["%s.%s"%(p,i)] = (pheromone, pathlength)
		elif p > 0:
			for n in nt2:
				args.terrain["%s.%s"%(p,n)] = (pheromone, pathlength)

	
def applyTerrainModification(args):
	"""
		Dependent on the input, this function modifies the terrain accordingly. 
		It reflects the pheromone and path length adjustment.
		Directly affected edges and their implicit neighbor edges are removed 
		from the graph.
	"""
	actGC = {}
	for i in args.GC:
		v, s1, s2 = i
		for j in range(s1-1,s2):
			actGC[j] = v
			
	dels = []
	for terrainelement in sorted(args.terrain):
		pheromone, pathlength = args.terrain[terrainelement]
		pheromone = getConstraint(terrainelement, args)
		pathlength = getConstraint(terrainelement, args)
		pathlength = applyGCcontributionPathAdjustment(pathlength, actGC, terrainelement)
		if pheromone * pathlength == 0: dels.append(terrainelement)
		args.terrain[terrainelement] = (pheromone, pathlength,[])
	
	further_dels = {}
	for terrainelement in sorted(dels):
		pos, nucs = terrainelement.split(".")
		if int(pos) < len(args.Cstr)-1:
			to_nt = nucs[-1:]
			successor_pos = int(pos) + 1
			for i in ["A", "C", "G", "U"]:
				del_element = str(successor_pos) + "." + to_nt + i
				further_dels[del_element] = 1
		further_dels[terrainelement] = 1
	# deleting the inbound and outbound edges, which are forbidden
	for terrainelement in further_dels:
		del args.terrain[terrainelement]
	# allocate the appropriate children of edges 
	for terrainelement in args.terrain:
		pheromone, pathlength, children = args.terrain[terrainelement]
		pos, nucs = terrainelement.split(".")
		if int(pos) < len(args.Cstr):
			to_nt = nucs[-1:]
			successor_pos = int(pos) + 1
			for i in ["A", "C", "G", "U"]:
				if str(successor_pos) + "." + to_nt + i in args.terrain:
					children.append(str(successor_pos) + "." + to_nt + i)
		args.terrain[terrainelement] = (pheromone, pathlength,children)
	# ADDING THE START EDGES
	starts = []
	for i in ["A", "C", "G", "U"]:
		if str(0) + "." + i in args.terrain:
			starts.append(str(0) + "." + i)
	args.terrain["00.XY"] = (1, 1, starts)

	
	
def applyGCcontributionPathAdjustment(pathlength, actGC, terrainelement):
	"""
		GC path length contribution calculation.
	"""	
	GCadjustment = 1.5
	minimum = 0.5
	upper = GCadjustment
	lower = minimum
	position , nts = terrainelement.split(".")
	nt = nts[-1:]
	
	if nt == "A" or nt == "U":
		pathlength = pathlength * maprange( (0, 1) , (lower, upper), actGC[int(position)])
	elif nt == "G" or nt == "C":
		pathlength = pathlength * maprange( (1, 0) , (lower, upper), actGC[int(position)])
	return pathlength
  
def printTerrain(terrain):
	"""
		Prints a given terrain
	"""
	#print sorted(terrain.keys())
	tmp_i = "0"
	tmp_c = 0
	terrain = terrain[0]
	
	for a, i in enumerate(sorted(terrain.keys())):
		#print a
		if i.split(".")[0] != tmp_i:
			print "\nElements:", tmp_c,"\n#########################\n", i, terrain[i]
			
			tmp_c = 1
			tmp_i = i.split(".")[0]
		else:
			print i, terrain[i]
			tmp_c += 1
			
	print "\nElements:", tmp_c
	print "#########################"
	print len(terrain)
	
########################################
# SEQUENCE ASSEMBLY a.k.a. The ANTWALK
########################################

def pickStep(tmp_steps, summe):
	"""
		Selects a step within the terrain
	"""
	if len(tmp_steps) == 1:
		return tmp_steps[0][1] # returning the nucleotide of the only present step
	else:
		rand = random.random() # draw random number
		mainval = 0
		for choice in xrange(len(tmp_steps)):
			val, label = tmp_steps[choice]
			mainval += val/float(summe)
			if mainval > rand: # as soon, as the mainval gets larger than the random value the assignment is done
				return label

def getPath(args):
	"""
		Performs a walk through the terrain and assembles a sequence, while respecting the structure constraint and IUPAC base complementarity
		of the base pairs GU, GC and AT
	"""
	nt = ["A","C","G","U"]
	prev_edge = "00.XY"
	sequence = ""
	while len(sequence) <  len(args.Cstr):
		coming_from = sequence[-1:]
		summe = 0
		steps = []
		i = len(sequence)
		allowed_nt = "ACGU"
		# base pair closing case check, with subsequent delivery of a reduced allowed nt set
		
		if i > args.BPstack[i][1][0]:
			jump =  args.BPstack[i][1][0]
			nuc_at_jump = sequence[jump]
			allowed_nt = args.IUPAC_reverseComplements[nuc_at_jump]

		# Checking for every possible nt if it is suitable for the selection procedure
		for edge in args.terrain[prev_edge][-1]:
			if edge[-1:] in allowed_nt:
				pheromone, PL , children = args.terrain[edge]
				value = ((float(pheromone * args.alpha)) + ((1/float(PL)) * args.beta))
				summe += value
				steps.append((value, edge))
		prev_edge = pickStep(steps, summe)
		sequence += prev_edge[-1:]
		
	return sequence
	
def getPathFromSelection(args, RNAfold, RNAfold_pattern):
	"""
		Returns the winning path from a selection of pathes...
	"""
	win_path = 0
	for i in xrange(args.ants_per_selection):
		# Generate Sequence
		sequence = getPath(args)
		# Measure sequence features and transform them into singular distances
		distance_structural = float(getStructuralDistance(args, sequence, RNAfold, RNAfold_pattern)) 
		distance_GC = float(getGCDistance(args.GC, sequence))
		distance_seq = float(getSequenceEditDistance(args.Cseq, sequence))
		# Calculate Distance Score
		D = distance_structural + distance_GC + distance_seq
      
		# SELECT THE BEST-OUT-OF-k-SOLUTIONS according to distance score
		if i == 0:
			win_path = (sequence, D, distance_structural, distance_GC, distance_seq)
		else:
			if D < win_path[1]:
				win_path = (sequence, D, distance_structural, distance_GC, distance_seq)
	return win_path
	
################################
# STRUCTURE PREDICTION METHODS
################################

def getPKStructure(sequence, args):
	"""
		Initialization pKiss mfe pseudoknot prediction
	"""
	p2p = "pKiss_mfe"
	p2p = "/usr/local/pkiss/2014-03-17/bin/pKiss_mfe"
	p = subprocess.Popen( ([p2p, '-T', str(args.temperature), '-s', args.strategy, sequence]),
				#shell = True,
				stdin = subprocess.PIPE,
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE,
				close_fds = True)

	# get linewise stdout output 
	pKiss_output = p.stdout.read().split("\n");

	# init with empty sequence
	structure = "." * len(sequence)
	if (len(pKiss_output) > 1) :
		structure = "".join(pKiss_output[1].split(" ")[3])
	else :
		print "DEBUG getPKStructure() : call stdout ="
		print pKiss_output
		print "DEBUG getPKStructure() : call stderr ="
                print p.stderr.read()
		print "DEBUG END (empty structure returned)"
		# close subprocess handle
		p.communicate()
		exit(-1)
	# close subprocess handle
	p.communicate()
	return structure

def getIPKnotStructure(sequence):
    """
        Initialization IPKnot pseudoknot prediction
    """
    f = str(uuid.uuid4())
    with open(f, 'w') as FASTA:
        FASTA.write(">tmp1\n")
        FASTA.write(sequence)
    p = subprocess.Popen( (['ipknot', f ]),
                #shell = True,
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                close_fds = False)
    IPKnot_output = p.stdout.read().split("\n");
    # init with empty sequence
    structure = "." * len(sequence)
    if (len(IPKnot_output) > 1) :
        #print IPKnot_output
        structure = IPKnot_output[2]
    else :
        print "DEBUG getHKStructure() : call stdout ="
        print IPKnot_output
        print "DEBUG getHKStructure() : call stderr ="
        print p.stderr.read()
        print "DEBUG END (empty structure returned)"
        # close subprocess handle
        p.communicate()
        exit(-1)
    p.communicate()
    os.remove(f)
    return structure

def getHKStructure(sequence, args):
    """
        Initialization HotKnots pseudoknot prediction
    """
    curr_dir = os.getcwd()
    os.chdir(args.HotKnots_PATH)
    args = 'HotKnots -m CC -s ' + sequence
    p = subprocess.Popen( (args),
                shell = True,
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                close_fds = False)

    HotKnots_output = p.stdout.read().split("\n");

    # init with empty sequence
    structure = "." * len(sequence)
    if (len(HotKnots_output) > 1) :
        structure = "".join(HotKnots_output[2].replace(" ", "#").replace("##","#").split("#")[1].split("\t")[0])
    else :
        print "DEBUG getHKStructure() : call stdout ="
        print HotKnots_output
        print "DEBUG getHKStructure() : call stderr ="
        print p.stderr.read()
        print "DEBUG END (empty structure returned)"
        # close subprocess handle
        p.communicate()
        exit(-1)
    p.communicate()	
    os.chdir(curr_dir)
    return structure


def init_RNAfold(version, temperature, paramFile = ""):
	"""
		Initialization RNAfold listener
	"""
	p2p = ""
	t = "-T " + str(temperature)
	P = ""
	if paramFile != "":
		P = "-P " + paramFile
	p2p = "RNAfold"
	p = subprocess.Popen( ([p2p, '--noPS', '-d 2', t, P]),
				#shell = True,
				stdin = subprocess.PIPE,
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE,
				close_fds = True)
	return p
		
def consult_RNAfold(seq, p):
	"""
		Consults RNAfold listener
	"""
	p.stdin.write(seq+'\n')
	out = ""
	for i in xrange(2):
		out += p.stdout.readline()
	return out


def getRNAfoldStructure(struct2, process1):
	"""
		Retrieves folded structure of a RNAfold call
	"""
  
	RNAfold_pattern = re.compile('.+\n([.()]+)\s.+')
	RNAfold_match = RNAfold_pattern.match(consult_RNAfold(struct2, process1))
	current_structure = ""
	return RNAfold_match.group(1)
  
  

###########################################
# DISTANCE MEASURES AND RELATED FUNCTIONS
###########################################
  
def getInducingSequencePositions(args):
	"""
		Delimiting the degree of structure inducement by the supplied sequence constraint.
		0 : no sequence induced structure constraint
		1 : "ACGT" induce structure (explicit nucleotide structure inducement level)
		2 : "MWKSYR" and "ACGT" (explicit and double instances)
		3 : "BDHV" , "MWKSYR" and "ACGT" (explicit, double, and triple instances)
	"""
	setOfNucleotides = "" # resembling the "0"-case
	if args.level == 1:
		setOfNucleotides = "ACGU"
	elif args.level == 2:
		setOfNucleotides = "ACGUMWKSYR"
	elif args.level == 3:
		setOfNucleotides = "ACGUMWKSYRBDHV"
		
	tmpSeq = ""
	listset = setOfNucleotides
	for pos in args.Cseq:
		if pos not in listset:
			tmpSeq += "N"
		else:
			tmpSeq += pos
	
	return setOfNucleotides, tmpSeq

def getBPDifferenceDistance(stack1, stack2):
	"""
		Based on the not identical amount of base pairs within both structure stacks
	"""
	d = 0
	for i in stack1.keys():
		# check base pairs in stack 1
		if i < stack1[i] and stack1[i] != stack2[i]:
			d += 1
		# check base pairs in stack 2
	for i in stack2.keys():
		if i < stack2[i] and stack1[i] != stack2[i]:
			d += 1
	return d

def getStructuralDistance(args, sequence, RNAfold, RNAfold_pattern):
	"""
		Calculator for Structural Distance
	"""
	# fold the current solution's sequence to obtain the structure
	current_structure = ""
	### Selection of specific folding mechanism
	if args.pseudoknots:
		if args.pkprogram == "pKiss":
			current_structure = getPKStructure(sequence, args)
		elif args.pkprogram == "HotKnots":
			current_structure = getHKStructure(sequence, args)
		elif args.pkprogram == "IPKnot":
			current_structure = getIPKnotStructure(sequence)
	else:
		RNAfold_match = RNAfold_pattern.match(consult_RNAfold(sequence, RNAfold))
		current_structure = RNAfold_match.group(1)

	# generate the current structure's base-pair stack
	bp = getbpStack(current_structure)[0]
	
	# deriving a working copy of the target structure# add case-dependend structural constraints in case of lonley basepairs formation
	tmp_target_structure_bp = getbpStack(args.Cstr)[0]
	
	### LONELY BASE PAIR MANAGEMENT
	# add case-dependend structural constraints in case of lonley basepairs formation
	if args.noLBPmanagement:
		for lp in args.LP:
			if bp[lp] == args.LP[lp]: # if the base pair is within the current solution structure, re-add the basepair into the constraint structure.
				tmp_target_structure_bp[lp] = args.LP[lp]
				tmp_target_structure_bp[args.LP[lp]] = lp
	###
	### "ABCDEFGHIJKLMNOPQRSTUVWXYZ" -> HARD CONSTRAINT
	### "abcdefghijklmnopqrstuvwxyz" -> SOFT CONSTRAINT
	###

	# FUZZY STRUCTURE CONSTRAINT MANAGEMENT - HARD CONSTRAINT
	# check for all allowed hard implicit constraint block declarators
	dsoftCstr = 0
	for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz":
		occurances = []
		for m in re.finditer(c, args.Cstr): # search for a declarator in the requested structure
			occurances.append(m.start()) # save the corresponding index
		# transform declarator into single stranded request
		for i in occurances: 
			tmp_target_structure_bp[i] = i
		# infer a base pair within the block declarated positions, if the current structure provides it.
		fuzzy_block_penalty = 1
		for i in occurances:
			for j in occurances:
				if i < j:
					if bp[i] == j:
						fuzzy_block_penalty = 0 # if a base pair is present within a current fuzzy_block, no penalty is risen
						tmp_target_structure_bp[i] = bp[i]
						tmp_target_structure_bp[bp[i]] = i
		if len(occurances) > 0:
			if c.isupper():
				dsoftCstr += fuzzy_block_penalty
			
	# INDUCES STRUCTURE MANAGEMENT
	# Checking Cseq influence and it's induced basepairs
	IUPACinducers, tmp_Cseq = getInducingSequencePositions(args)
	if len(args.Cseq.strip("N")) > 0:
    #print "Processing Cseq influence"
		# Iterate over all positions within the Base Pair stack
		for i in args.BPstack: # Check for each base index i 
			if i < bp[i]: # if the current index is samller that the affiliated in the basepair stack of the current solution
				bp_j = bp[i] # Actual j index of the current solution
				BP_j = args.BPstack[i][1][0] # j index of the requested structure
				if (i != bp_j and i == BP_j and args.BPstack[i][0] in IUPACinducers ): # if i pairs with some other base in the current structure, and i is requested single stranded and the Sequence constraint is allowed to induce...
					if (args.BPstack[bp_j][1][0] == bp_j and args.BPstack[bp_j][0] in IUPACinducers):# If position j is requested singlestranded and position j nucleotide can induce base pairs
						#if isCompatible(bp[i][0], bp[i][1][1], IUPAC_compatibles): # If both nucleotides, i and j are actually compatible
						tmp_target_structure_bp[i] = bp[i]
						tmp_target_structure_bp[bp_j] = i
						
	dsreg = getBPDifferenceDistance(tmp_target_structure_bp, bp)
		
	# CHECK FOR ALL DETERMINED LONELY BASE PAIRS (i<j), if they are formed
	failLP = 0
	if args.noLBPmanagement:
		for lp in args.LP: 
			if bp[lp] != args.LP[lp]: 
				isComp = isCompatible(sequence[lp],sequence[args.LP[lp]], args.IUPAC_compatibles)
				isStru = isStructureCompatible(lp, args.LP[lp] ,bp)
				if not ( isStru and isStru ): # check if the bases at the specific positions are compatible and check if the 
					# basepair can be formed according to pseudoknot free restriction. If one fails, a penalty distance is raised for that base pair
					failLP += 1

	dsLP = float(failLP)
	return (dsreg + dsLP + dsoftCstr) /float(len(tmp_target_structure_bp))
  
def getGC(sequence):
	"""
		Calculate GC content of a sequence
	"""
	GC = 0
	for nt in sequence:
		if nt == "G" or nt == "C":
			GC = GC + 1
	GC = GC/float(len(sequence))
	return GC
	
def getGCDistance(GC, sequence):
	"""
	Calculate the pseudo GC content distance 
	"""
	D = 0
	for tGC in GC:
		v, s1, s2 = tGC # (tGC, start, stop)
		tmp_seq = sequence[s1 : s2 + 1]
		L = len(tmp_seq)
		gc = getGC(tmp_seq)
		nt_coeff = L * v
		pc_nt = (1/float(L))*100

		d = gc - v
		d = d * 100
		f = math.floor(nt_coeff)
		c = math.ceil(nt_coeff)

		if d < 0: 
			d = d + (abs(nt_coeff - f)) * pc_nt
		elif d > 0: 
			d = d - abs(nt_coeff - c) * pc_nt
		elif d == 0:
			pass
		
		d = abs(round(d, 7))
		D += d
	return D
	

def getSequenceEditDistance(SC, path):
	"""
	Calculate sequence edit distance of a solution to the constraint
	"""
	IUPAC = {"A":"A", "C":"C", "G":"G", "U":"U", "R":"AG", "Y":"CU", "S":"GC", "W":"AU","K":"GU", "M":"AC", "B":"CGU", "D":"AGU", "H":"ACU", "V":"ACG", "N":"ACGU"}         
	edit = 0
	for i in xrange(len(SC)):
		if path[i] not in IUPAC[SC[i].upper()]:
			edit += 1
	return edit/float(len(path))


#########################
# TERRAIN GRAPH EDITING
#########################

def evaporate(args): 
	"""
	Evaporate the terrain's pheromone trails
	"""
	c = 1
	for key in args.terrain:
		p,l,c = args.terrain[key]
		p *= (1-args.ER)
		args.terrain[key] = (p, l, c)
		
def trailBlaze(sequence, current_structure, ds, dgc, dseq, args):
	"""
		Pheromone Update function accorinding to the quality of the solution
	"""
	bpstack, LP = getbpStack(current_structure)
	bs = updateValue(ds, args.Cstrweight, args.omega)
	bGC = updateValue(dgc, args.Cgcweight, args.omega)
	bSeq = updateValue(dseq, args.Cseqweight, args.omega)
	d = bs + bGC + bSeq

	transitions = getTransitions(sequence)

	for trans in xrange(len(transitions)): # for each transition in the path
		id1 = int(transitions[trans].split(".")[0])
		tar_id2 = int(args.BPstack[id1][1][0]) # getting requested  position id2
		curr_id2 = int(bpstack[id1]) # getting the current situation
		multiplicator = 0
		if tar_id2 == curr_id2 and id1 != tar_id2 and id1 != curr_id2: # case of a base pair, having both brackets on the correct position
			multiplicator = 1
		elif tar_id2 == curr_id2 and id1 == tar_id2 and id1 == curr_id2: # case of a single stranded base in both structures
			multiplicator = 1
		p, l, c = args.terrain[transitions[trans]] # getting the pheromone and the length value of the single path transition
		p +=  d * multiplicator
		args.terrain[transitions[trans]] = (p, l, c) # updating the values wihtin the terrain's


def getTransitions(p):
	"""
		Retreive transitions of a specific path/sequence
	"""
	transitions = []
	for pos in xrange(len(p)):
		if pos == 0:
			transitions.append(str(pos) + "." + p[pos])

		else:
			insert = p[pos-1] + p[pos]
			transitions.append(str(pos) + "." + insert)

	return transitions

def updateValue(distance, correction_term, omega):
	"""
	Retrieves a distance dependend pheromone value
	"""
	if correction_term == 0:
		return 0
	else:
		if distance == 0:
			return omega * correction_term
		else:
			return (1/float(distance)) * correction_term
      
def updateTerrain(sequence, current_structure, ds, dgc, dseq, args):
	"""
		General updating function
	"""
	evaporate(args)
	trailBlaze(sequence, current_structure, ds, dgc, dseq, args)

	
	
#######################
# CONVERGENCE MEASURE 
#######################

def inConvergenceCorridor(d_struct, d_gc, d_seq, BS_d_struct, BS_d_gc, BS_d_seq):
	"""
		Check if a solutions qualities are within the convergence corridor
	"""
	struct_var = (BS_d_struct + BS_d_struct/float(100) * 5)
	gc_var = (BS_d_gc + BS_d_gc/float(100) * 5)
	seq_var = (BS_d_seq + BS_d_seq/float(100) * 5)
	if d_struct <= struct_var and d_gc <= gc_var and d_seq <= seq_var:
		return True
	else:
		return False



##############################
# TARGET GC VALUE MANAGEMENT
##############################

def getGCSamplingValue(GC, tGCmax, tGCvar):
	"""
	Returns a suitable GC value, dependend on the user input: Either returning the single GC value,
	which the user entered, or a smpled GC value
	from a designated distribution in it's interavals
	"""
	returnval = 0
	if tGCmax == -1.0 and tGCvar == -1.0: # regular plain tGC value as requested 
		return GC
	elif tGCmax != -1.0 and tGCvar == -1.0: # uniform distribution tGC value sampling
		if GC < tGCmax:
			tmp_GC = tGCmax
			tGCmax = GC
			GC = tmp_GC
		while returnval <= 0:
			returnval = float(numpy.random.uniform(low=GC, high=tGCmax, size=1))
		return returnval
	elif tGCmax == -1.0 and tGCvar != -1.0: # normal distribution tGC value sampling
		while returnval <= 0:
			returnval = float(numpy.random.normal(GC, tGCvar, 1))
		return returnval
	

  
  






			



	
def exe():
	"""
	MAIN EXECUTABLE WHICH PARSES THE INPUT COMMAND LINE
	"""

	argument_parser = argparse.ArgumentParser(
	fromfile_prefix_chars='@',
	description = """
    
	#########################################################################
	#       antaRNA - ant assembled RNA                                     #
	#       -> Ant Colony Optimized RNA Sequence Design                     #
	#       ------------------------------------------------------------    #
	#       Robert Kleinkauf (c) 2015                                       #
	#       Bioinformatics, Albert-Ludwigs University Freiburg, Germany     #
	#########################################################################
  
	- For antaRNA only the VIENNNA RNA Package must be installed on your linux system.
	  antaRNA will only check, if the executables of RNAfold  of the ViennaRNA package can be found. If those programs are 
	  not installed correctly, no output will be generated, an also no warning will be prompted.
	  So the binary path of the Vienna Tools must be set up correctly in your system's PATH variable in order to run antaRNA correctly!
	- If you want to use the pseudoknot functionality, pKiss_mfe of the RNAshapes studio OR HotKnots OR IPKnot must be installed and callable as standalone in order to execute antaRNA.
    
    - antaRNA was only tested under Linux.
    
    - For questions and remarks please feel free to contact us at http://www.bioinf.uni-freiburg.de/
	
	""",
	
	epilog = """   
	Example calls:
		python antaRNA_vXY.py -Cstr "...(((...)))..." -tGC 0.5 -n 2 -v
		python antaRNA_vXY.py -Cstr ".........aaa(((...)))aaa........." -tGC 0.5 -n 10 --output_file /path/to/antaRNA_TESTRUN -v
		python antaRNA_vXY.py -Cstr "BBBBB....AAA(((...)))AAA....BBBBB" -Cseq "NNNNANNNNNCNNNNNNNNNNNGNNNNNNUNNN" --tGC 0.5 -n 10

	#########################################################################
	#       --- Hail to the Queen!!! All power to the swarm!!! ---          #
	#########################################################################
		""",
	formatter_class=RawTextHelpFormatter
	)
	
	argument_parser.convert_arg_line_to_args = convert_arg_line_to_args
	
	constraints = argument_parser.add_argument_group('Constraint Variables', 'Use to define an RNA constraint system.')
	constraints.add_argument("-Cstr", "--Cstr", 
								help="Structure constraint using RNA dotbracket notation with fuzzy block constraint. \n(TYPE: %(type)s)\n\n", 
								type=str, 
								required=True)
								
	constraints.add_argument("-Cseq", "--Cseq", 
								help="Sequence constraint using RNA nucleotide alphabet {A,C,G,U} and wild-card \"N\". \n(TYPE: %(type)s)\n\n", 
								type=str, 
								default = "") 
								
	constraints.add_argument("-l", "--level", 
								help="Sets the level of allowed influence of sequence constraint on the structure constraint [0:no influence; 3:extensive influence].\n(TYPE: %(type)s)\n\n", 
								type=int, 
								default = 1)

	constraints.add_argument("-tGC", "--tGC", 
								help="Objective target GC content in [0,1].\n(TYPE: %(type)s)\n\n", 
								type=parseGC, 
								required=True,
								action = 'append')
								
	constraints.add_argument("-tGCmax", "--tGCmax", 
								help = "Provides a maximum tGC value [0,1] for the case of uniform distribution sampling. The regular tGC value serves as minimum value.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=float, 
								default=-1.0)
	
	constraints.add_argument("-tGCvar", "--tGCvar", 
								help = "Provides a tGC variance (sigma square) for the case of normal distribution sampling. The regular tGC value serves as expectation value (mu).\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=float, 
								default=-1.0)

	constraints.add_argument("-T", "--temperature", 
								help = "Provides a temperature for the folding algorithms.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=float, 
								default=37.0)
	
	constraints.add_argument("-P", "--paramFile", 
								help = "Changes the energy parameterfile of RNAfold. If using this explicitly, please provide a suitable energy file delivered by RNAfold. \n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=str, 
								default="")
								
	constraints.add_argument("-noGU", "--noGUBasePair", 
								help="Forbid GU base pairs. \n\n", 
								action="store_false")
								
	constraints.add_argument("-noLBP", "--noLBPmanagement", 
								help="Disallowing antaRNA lonely base pair-management. \n\n", 
								action="store_false", 
								default =True)

	pk = argument_parser.add_argument_group('Pseudoknot Variables', 'Use in order to enable pseudoknot calculation. pKiss_mfe needs to be installed.')
	pk.add_argument("-p", "--pseudoknots", 
								help = "Switch to pseudoknot based prediction using pKiss. Check the pseudoknot parameter usage!!!\n\n", 
								action="store_true")
	
	pk.add_argument("-pkP", "--pkprogram",
								help = "Select a pseudoknot prediction program.\nIf HotKnots is used, please specify the bin folder of Hotknots with absolute path using HK_PATH argument.\n(DEFAULT: %(default)s, TYPE: %(type)s, Choice: [pKiss|HotKnots|IPKnot])\n\n",
								type=str,
								default="pKiss")
	pk.add_argument("-pkPar", "--pkparameter", 
								help = "Enable optimized parameters for the usage of pseudo knots (Further parameter input ignored).\n\n", 
								action="store_true")	
	
	pk.add_argument("-HKPATH", "--HotKnots_PATH",
								help = "Set HotKnots absolute path, like /path/to/HotKnots/bin.\nIf HotKnots is used, please specify the bin folder of Hotknots with absolute path using HK_PATH argument.\n(DEFAULT: %(default)s, TYPE: %(type)s\n\n",
								type=str,
								default="")

	pk.add_argument("--strategy", 
								help = "Defining the pKiss folding strategy.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=str, 
								default="A")


	output = argument_parser.add_argument_group('Output Variables', 'Tweak form and verbosity of output.')
	output.add_argument("-n", "--noOfColonies", 
								help="Number of sequences which shall be produced. \n(TYPE: %(type)s)\n\n", 
								type=int,  
								default=1)

	output.add_argument("-of","--output_file", 
								help="Provide a path and an output file, e.g. \"/path/to/the/target_file\". \n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=str, 
								default="STDOUT")
	
	output.add_argument("-rPY","--py", 
								help="Switch on PYTHON compatible behavior. \n(DEFAULT: %(default)s)\n\n", 
								action="store_true") 
	
	output.add_argument("--name", 
								help="Defines a name which is used in the sequence output. \n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=str, 
								default="antaRNA")
								
	output.add_argument("-v", "--verbose", 
								help="Displayes intermediate output.\n\n", 
								action="store_true") 
	
	output.add_argument("-ov", "--output_verbose", 
								help="Prints additional features and stats to the headers of the produced sequences. Also adds the structure of the sequence.\n\n", 
								action="store_true")



	aco = argument_parser.add_argument_group('Ant Colony Variables', 'Alter the behavior of the ant colony optimization.')
	aco.add_argument("-s", "--seed", 
								help = "Provides a seed value for the used pseudo random number generator.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=str, 
								default="none")
								
	aco.add_argument("-ip", "--improve_procedure", 
								help = "Select the improving method.  h=hierarchical, s=score_based.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=str, 
								default="s")  
								
	aco.add_argument("-r", "--Resets", 
								help = "Amount of maximal terrain resets, until the best solution is retuned as solution.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=int, 
								default=5)
	
	aco.add_argument("-aps", "--ants_per_selection", 
								help = "best out of k ants.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=int, 
								default=10)

	aco.add_argument("-CC", "--ConvergenceCount", 
								help = "Delimits the convergence count criterion for a reset.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=int, 
								default=130)  
	
	aco.add_argument("-aTC", "--antsTerConv", 
								help = "Delimits the amount of internal ants for termination convergence criterion for a reset.\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=int, 
								default=50)
	aco.add_argument("-a", "--alpha", 
								help="Sets alpha, probability weight for terrain pheromone influence. [0,1] \n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=float, 
								default=1.0)
	
	aco.add_argument("-b", "--beta", 
								help="Sets beta, probability weight for terrain path influence. [0,1]\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=float, 
								default=1.0)
	
	aco.add_argument("-er", "--ER", 
								help="Pheromone evaporation rate. \n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=float, 
								default=0.2)
	
	aco.add_argument("-Cstrw", "--Cstrweight", 
								help="Structure constraint quality weighting factor. [0,1]\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=float, 
								default=0.5)
	
	aco.add_argument("-Cgcw", "--Cgcweight", 
								help="GC content constraint quality weighting factor. [0,1]\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n", 
								type=float, 
								default=5.0)
	
	aco.add_argument("-Cseqw", "--Cseqweight", 
								help="Sequence constraint quality weighting factor. [0,1]\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n\n", 
								type=float, 
								default=1.0)
	
	aco.add_argument("-o", "--omega", 
								help="Scoring parameter\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n\n", 
								type=float, 
								default=2.23)

	aco.add_argument("-t", "--time", 
								help="Limiting runtime [seconds]\n(DEFAULT: %(default)s, TYPE: %(type)s)\n\n\n", 
								type=int, 
								default=600)		
	
	argparse_arguments = argument_parser.parse_args()

	
	hill = AntHill()
	hill.params.readArgParseArguments(argparse_arguments)
	hill.params.py = False
	hill.params.varCheck()
	hill.findSequence()
	
	
	#antaRNA_variables = AntaRNAVariables()
	#antaRNA_variables.readArgParseArguments(argparse_arguments)
	#antaRNA_variables.varCheck()
	#findSequence(antaRNA_variables)

	

#########################
### CLASSES
#########################
class AntHill:
	"""
		antaRNA AntHill. Can perform varoius actions! :)
	"""
	
	def __init__(self):
	
		self.params = AntaRNAVariables()
		self.tmp_sequence = ""
		self.tmp_structure = ""
		self.tmp_stats = ""
		self.tmp_result = ""
		self.result = []
		
		
	###################################################  
	#  PRINCIPAL EXECUTION ROUTINES OF THE ACO METHOD
	###################################################

	def swarm(self):
		"""
			Execution function of a single ant colony finding one solution sequence
		"""
		
		# Constraint Checks prior to execution
		checkSimilarLength(self.params.Cstr, self.params.Cseq)
		isValidStructure(self.params.Cstr)
		checkSequenceConstraint(self.params.Cseq)
		
		self.params.BPstack , self.params.LP = getBPStack(self.params)
		checkConstaintCompatibility(self.params)

		retString = ""
		retString2 = []

		start_time = time.time()

		# INITIALIZATION OF Vienna RNAfold
		RNAfold = init_RNAfold(213, self.params.temperature, self.params.paramFile)
		RNAfold_pattern = re.compile('.+\n([.()]+)\s.+')
		
		initTerrain(self.params) 
		applyTerrainModification(self.params)
		global_ant_count = 0
		global_best_ants = 0
		criterion = False
		met = True  
		ant_no = 1
		prev_res = 0
		seq = ""

		counter = 0
		
		dstruct_log = []
		dGC_log = []
			
		convergence_counter = 0
		
		resets = 0
		path = ""
		curr_structure = ""

		Dscore = 100000
		ds = 10000
		dGC = 10000
		dseq = 10000
		best_solution = (path, curr_structure, Dscore, ds, dGC, dseq)
		best_solution_local = (path, curr_structure, Dscore, ds, dGC, dseq)
		
		best_solution_since = 0
		
		
		
		# IN CASE OF LP-MANAGEMENT
		if self.params.noLBPmanagement:
			if len(self.params.LP) > 0 :
				for lp in self.params.LP:
					self.params.Cstr = substr(lp + 1, self.params.Cstr, ".")
					self.params.Cstr = substr(self.params.LP[lp] + 1, self.params.Cstr, ".")

		init = 1
		used_time = getUsedTime(start_time)
		while criterion != met and used_time < self.params.time:
			iteration_start = time.time()
			global_ant_count += 1
			global_best_ants += 1

			path_info = getPathFromSelection(self.params, RNAfold, RNAfold_pattern)

			distance_structural_prev = ds
			distance_GC_prev = dGC
			distance_seq_prev = dseq

			sequence, Dscore , ds, dGC, dseq = path_info
			curr_structure = ""
			if self.params.pseudoknots:
				if self.params.pkprogram == "pKiss":
					curr_structure = getPKStructure(sequence, self.params)
				elif self.params.pkprogram == "HotKnots":
					curr_structure = getHKStructure(sequence, self.params)
				elif self.params.pkprogram == "IPKnot":
					curr_structure = getIPKnotStructure(sequence)
			else:
				curr_structure = getRNAfoldStructure(sequence, RNAfold)
				
			curr_solution = (sequence, curr_structure, Dscore, ds, dGC, dseq)
			# BEST SOLUTION PICKING
			if self.params.improve_procedure == "h": # hierarchical check
				# for the global best solution
				if ds < best_solution[3] or (ds == best_solution[3] and dGC < best_solution[4]):
					best_solution = curr_solution
					ant_no = 1
				# for the local (reset) best solution
				if ds < best_solution_local[3] or (ds == best_solution_local[3] and dGC < best_solution_local[4]):
					best_solution_local = curr_solution
					
			elif self.params.improve_procedure == "s": #score based check
				# store best global solution
				if Dscore < best_solution[2]:
					best_solution = curr_solution
					ant_no = 1
				# store best local solution for this reset
				if Dscore < best_solution_local[2]:
					best_solution_local = curr_solution

			if self.params.verbose:
				print "SCORE " + str(Dscore) + " Resets " + str(resets) + " #Ant " + str(global_ant_count) + " out of " + str(self.params.ants_per_selection)  + " cc " + str(convergence_counter)

				print self.params.Cstr, " <- target struct" 
				print best_solution[0] , " <- BS since ", str(best_solution_since), "Size of Terrrain:", len(self.params.terrain)
				print best_solution[1] , " <- BS Dscore " + str(best_solution[2]) + " ds " + str(best_solution[3]) + " dGC " + str(best_solution[4]) + " dseq " + str(best_solution[5])+ " LP " + str(len(self.params.LP)) + " <- best solution stats"
				print curr_structure, " <- CS"
				print sequence,
				print " <- CS", "Dscore", str(Dscore), "ds", ds, "dGC", dGC, "GC", getGC(sequence)*100, "Dseq", dseq

			#### UPDATING THE TERRAIN ACCORDING TO THE QUALITY OF THE CURRENT BESTO-OUT-OF-k SOLUTION
			updateTerrain(sequence, curr_structure, ds, dGC, dseq, self.params) 

			if self.params.verbose: print "Used time for one iteration", time.time() - iteration_start
				
				
			# CONVERGENCE AND TERMINATION CRITERION MANAGEMENT
			if inConvergenceCorridor(curr_solution[3], curr_solution[4], curr_solution[5], best_solution_local[3], best_solution_local[4], best_solution_local[5]):
				convergence_counter += 1
			if distance_structural_prev == ds and distance_GC_prev == dGC and distance_seq_prev == dseq:
				convergence_counter += 1

			if best_solution[3] == self.params.objective_to_target_distance:
				if best_solution[4] == 0.0:
					if best_solution[5] == 0.0:
						break
				ant_no = ant_no + 1
				convergence_counter -= 1
			else:
				ant_no = 1

			if ant_no == self.params.antsTerConv or resets >= self.params.Resets or global_ant_count >= 100000 or best_solution_since == 5:
				break

			# RESET
			if ant_no < self.params.antsTerConv and convergence_counter >= self.params.ConvergenceCount:

				initTerrain(self.params)
				applyTerrainModification(self.params)
				criterion = False
				met = True  
				ant_no = 1
				prev_res = 0
				sequence = ""
				curr_structure = ""
				counter = 0
				Dscore = 100000
				ds = 10000
				dGC = 10000
				dseq = 10000
				best_solution_local = (sequence, curr_structure, Dscore, ds, dGC, dseq)

				convergence_counter = 0

				if resets == 0:
					sentinel_solution = best_solution
					best_solution_since += 1
				else:
					if best_solution[2] < sentinel_solution[2]:
						sentinel_solution = best_solution
						best_solution_since = 0
					else:
						best_solution_since += 1

				resets += 1
			
			used_time = getUsedTime(start_time)
			
		duration  = used_time

		self.result_stats += "|Ants:" + str(global_ant_count)
		self.result_stats += "|Resets:" + str(resets) + "/" + str(self.params.Resets)
		self.result_stats += "|AntsTC:" + str(self.params.antsTerConv) 
		self.result_stats += "|CC:" + str(self.params.ConvergenceCount) 
		self.result_stats += "|IP:" + str(self.params.improve_procedure) 
		self.result_stats += "|BSS:" + str(best_solution_since)
		self.result_stats += "|LP:" + str(len(self.params.LP))
		self.result_stats += "|ds:" + str(getStructuralDistance(self.params, sequence, RNAfold, RNAfold_pattern))
		self.result_stats += "|dGC:" + str(best_solution[4])
		self.result_stats += "|GC:" + str(getGC(sequence)*100)
		self.result_stats += "|dseq:" + str(getSequenceEditDistance(self.params.Cseq, sequence))
		self.result_stats += "|L:" + str(len(sequence))
		self.result_stats += "|Time:" + str(duration)

		self.result_sequence = best_solution[0]
		self.result_structure = best_solution[1]		

		# CLOSING THE PIPES TO THE PROGRAMS
		if (RNAfold is not None) :
			RNAfold.communicate()


	
	def findSequence(self):
		"""
			MAIN antaRNA - ant assembled RNA
		"""

		if self.params.seed != "none":
			random.seed(self.params.seed)
		
		print_to_STDOUT = (self.params.output_file == "STDOUT")

		if self.params.py == False:
			if print_to_STDOUT == False:
				outfolder = '/'.join(self.params.output_file.strip().split("/")[:-1])
				curr_dir = os.getcwd()
				if not os.path.exists(outfolder):
					os.makedirs(outfolder)
				os.chdir(outfolder)  
			
		self.params.Cseq = transform(self.params.Cseq)
	  
		# Allowed deviation from the structural target:
		self.params.objective_to_target_distance = 0.0

		# Loading the IUPAC copatibilities of nuleotides and their abstract representing symbols
		self.params.IUPAC = {"A":"A", "C":"C", "G":"G", "U":"U", "R":"AG", "Y":"CU", "S":"GC", "W":"AU","K":"GU", "M":"AC", "B":"CGU", "D":"AGU", "H":"ACU", "V":"ACG", "N":"ACGU"}         
		self.params.IUPAC_compatibles = loadIUPACcompatibilities(self.params)

		
		if self.params.noGUBasePair == True: ## Without the GU basepair
			self.params.IUPAC_reverseComplements = {"A":"U", "C":"G", "G":"C", "U":"A", "R":"UC", "Y":"AG", "S":"GC", "W":"UA","K":"CA", "M":"UG", "B":"AGC", "D":"ACU", "H":"UGA", "V":"UGC", "N":"ACGU"}         
		else: ## allowing the GU basepair
			self.params.IUPAC_reverseComplements = {"A":"U", "C":"G", "G":"UC", "U":"AG", "R":"UC", "Y":"AG", "S":"UGC", "W":"UAG","K":"UCAG", "M":"UG", "B":"AGCU", "D":"AGCU", "H":"UGA", "V":"UGC", "N":"ACGU"}         
		
		for col in xrange(self.params.noOfColonies):
			# Checking the kind of taget GC value should be used
			self.params.GC = []
			if len(self.params.tGC) == 1:
				self.params.GC.append((getGCSamplingValue(self.params.tGC[0][0], self.params.tGCmax, self.params.tGCvar), self.params.tGC[0][1], self.params.tGC[0][2]))
			else:
				self.params.GC = self.params.tGC

			# Actual execution of a ant colony procesdure
			swarm()

			# Post-Processing the output of a ant colony procedure
			self.tmp_result = ">" + self.params.name + "#" + str(col)
			if self.params.output_verbose:
				
				GC_out = ""
				for i in self.params.GC:
					v, s1, s2 = i
					GC_out += str(s1) + "-" + str(s2) + ">" + str(v) + ";"
				GC_out = GC_out[:-1]
				
				self.tmp_result += "|Cstr:" + self.params.Cstr + "|Cseq:" + self.params.Cseq + "|Alpha:" + str(self.params.alpha) + "|Beta:" + str(self.params.beta) + "|tGC:" + str(GC_out) + "|ER:" + str(self.params.ER) + "|Struct_CT:" + str(self.params.Cstrweight) + "|GC_CT:" + str(self.params.Cgcweight) + "|Seq_CT:" + str(self.params.Cseqweight) + self.tmp_stats + "\n" + self.tmp_sequence "\n" + self.tmp_structure + "\n"
			else:
				self.tmp_result += "\n" + self.result[0]
			if self.params.py == False:
				if print_to_STDOUT:
					print self.tmp_result
				else:
					if col == 0:
						print2file(self.params.output_file, self.tmp_result, 'w')
					else:
						print2file(self.params.output_file, self.tmp_result, 'a')
			else:
				self.result.append(self.tmp_result)
		if print_to_STDOUT == False:    
			os.chdir(curr_dir)
  
	
class AntaRNAVariables:
	"""
		antaRNA Variables management.
	"""
	def __init__(self):
		self.Cstr = ""
		self.Cseq = ""
		self.tGC = []
		self.level = 1
		self.tGCmax = -1.0
		self.tGCvar = -1.0
		self.temperature = 37.0
		self.paramFile = ""
		self.noGUBasePair = False
		self.noLBPmanagement = True
		self.pseudoknots = False
		self.pkprogram = "pKiss"
		self.pkparameter = False
		self.HotKnots_PATH = ""
		self.strategy = "A"
		self.noOfColonies = 1
		self.output_file = "STDOUT"
		self.py = True
		self.name="antaRNA"
		self.verbose = False 
		self.output_verbose = False
		self.seed ="none"
		self.improve_procedure = "s"
		self.Resets = 5
		self.ants_per_selection = 10
		self.ConvergenceCount = 130
		self.antsTerConv = 50
		self.alpha = 1.0
		self.beta = 1.0
		self.ER = 0.2
		self.Cstrweight = 0.5
		self.Cgcweight = 5.0
		self.Cseqweight = 1.0
		self.omega = 2.23
		self.time = 600
		
	def readArgParseArguments(self, args):
		self.Cstr = args.Cstr
		self.Cseq = args.Cseq
		self.tGC = args.tGC
		self.level = args.level
		self.tGCmax = args.tGCmax
		self.tGCvar = args.tGCvar
		self.temperature = args.temperature
		self.paramFile = args.paramFile
		self.noGUBasePair = args.noGUBasePair
		self.noLBPmanagement = args.noLBPmanagement
		self.pseudoknots = args.pseudoknots
		self.pkprogram = args.pkprogram
		self.pkparameter = args.pkparameter
		self.HotKnots_PATH = args.HotKnots_PATH
		self.strategy = args.strategy
		self.noOfColonies = args.noOfColonies
		self.output_file = args.output_file
		self.py = args.py
		self.name = args.name
		self.verbose = args .verbose
		self.output_verbose = args.output_verbose
		self.seed = args.seed
		self.improve_procedure = args.improve_procedure
		self.Resets = args.Resets
		self.ants_per_selection = args.ants_per_selection
		self.ConvergenceCount = args.ConvergenceCount
		self.antsTerConv = args.antsTerConv
		self.alpha = args.alpha
		self.beta = args.beta
		self.ER = args.ER
		self.Cstrweight = args.Cstrweight
		self.Cgcweight = args.Cgcweight
		self.Cseqweight = args.Cseqweight
		self.omega = args.omega
		self.time = args.time

	def varCheck(self):
		"""
			CHECK THE COMMAND LINE STUFF
		"""

		if self.Cseq == "":
			self.Cseq = "N" * len(self.Cstr)

		self.parse_GC_management()

		self.checkForViennaTools()
		if self.pseudoknots:
			if self.pkprogram == "pKiss":
				self.checkForpKiss()
				if self.pkparameter == True:
					self.alpha = 1.0
					self.beta = 0.1
					self.ER = 0.2 
					self.Cstrweight = 0.1 
					self.Cgcweight = 1.0 
					self.Cseqweight = 0.5 
					self.Cseqweight = 50 
					self.ConvergenceCount = 100
				
			elif self.pkprogram == "HotKnots" and self.HotKnots_PATH != "":
				self.checkForHotKnots(args)
				if self.pkparameter == True:
					self.alpha = 1.0
					self.beta = 0.1
					self.ER = 0.2 
					self.Cstrweight = 0.1 
					self.Cgcweight = 1.0 
					self.Cseqweight = 0.5 
					self.Cseqweight = 50 
					self.ConvergenceCount = 100
				
			elif self.pkprogram == "IPKnot":
				self.checkForIPKnot()
				if self.pkparameter == True:
					self.alpha = 1.0
					self.beta = 0.1
					self.ER = 0.2 
					self.Cstrweight = 0.1 
					self.Cgcweight = 1.0 
					self.Cseqweight = 0.5 
					self.Cseqweight = 50 
					self.ConvergenceCount = 100	
			else:
				print " Please choose a suitable pseudoknot predictor: [pKiss|Hotknots|IPKnot]"
				exit(1)
	def reachableGC(self):
		"""
			Checks if a demanded GC target content is reachable in dependence with the given sequence constraint.
			For each explicit and ambiguous character definition within the constraint, the respective possibilities
			are elicited: "A" counts as "A", but for "B", also an "U" is possible, at the same time, "G" or "C" are possible
			as well. So two scenarios can be evaluated: A minimum GC content, which is possible and a maximum GC content.
			For this, for all not "N" letters their potential towards  their potentials is evaluated and counted in the
			respective counters for min and max GC content. Only those characters are taken into concideration which would enforce
			a categorial pyrimidine/purine decision. (ACGU, SW)
			
		"""
		nucleotide_contribution = 1/float(len(self.Cseq)) * 1
		
		minGC = 0.0
		maxGC = 1.0
		for i in self.Cseq:
			if i != "N":
				if i == "A" or i == "U":
					maxGC -= nucleotide_contribution
				elif i == "C" or i == "G":
					minGC += nucleotide_contribution
				elif i == "S":#(G or C)
					minGC += nucleotide_contribution
				elif i == "W":#(A or T/U):
					maxGC -= nucleotide_contribution
		return (minGC, maxGC)	
	
	
	def parse_GC_management(self):

		if len(self.tGC) == 1 and type(self.tGC[0]) is float: # CASE Only one tGC value is defined, which needs to account for the whole terrain
			tgc = self.tGC.pop()
			self.tGC.append((tgc, 1, len(self.Cstr)))
			
		for t in self.tGC:
			if len(t) != 3:
				print "Error :: Not enough tGC and affiliated areas declarations"
				exit(1)
				
		check_set = set(range(1,len(self.Cstr) + 1))
		curr_set = set()
		for i, area in enumerate(self.tGC): # CHECK if the areas are consistent and do not show disconnectivity.
			v, s1, s2 = area
			if i < 0 or i > 1:
				print "Error: Chosen tGC > %s < not in range [0,1]" % (i)
				exit(1)
			tmp_set = set(range(int(s1), int(s2 + 1)))
			if len(curr_set.intersection(tmp_set)) == 0:
				curr_set = curr_set.union(tmp_set)
			else: 
				print "Error: Double defined tGC declaration area sector detected. Nucleotide positions", ", ".join(str(e) for e in curr_set.intersection(tmp_set)), "show(s) redundant tGC declaration"
				exit(1)
		if len(curr_set.symmetric_difference(check_set)) != 0:
			print "Error: Undefined tGC area sectors detected. Nucleotide positions", ", ".join(str(e) for e in curr_set.symmetric_difference(check_set)), "is/are not covered."
			exit(1)
			
		for tgc in self.tGC: # CHECK if the specified GC values can be reached at all...
			v, start, stop = tgc
			tmp_sc = self.Cseq[start:stop + 1]
			minGC, maxGC = self.reachableGC()
			if v > maxGC or v < minGC:
				print >> sys.stderr, "WARNING: Chosen target GC %s content is not reachable. The selected sequence constraint contradicts the tGC constraint value." % (v) 
				print >> sys.stderr, "Sequence Constraint allows tGC only to be in [%s,%s]" % (minGC, maxGC) 
				exit (1)


	
	##########################
	# PROGRAM PRESENCE CHECK
	##########################



	def checkForViennaTools(self):
		"""
		Checking for the presence of the Vienna tools in the system by which'ing for RNAfold and RNAdistance
		"""
		RNAfold_output = subprocess.Popen(["which", "RNAfold"], stdout=subprocess.PIPE).communicate()[0].strip()
		if len(RNAfold_output) > 0 and RNAfold_output.find("found") == -1 and RNAfold_output.find(" no ") == -1:
			return True
		else:
			print "It seems the Vienna RNA Package is not installed on your machine. Please do so!"
			print "You can get it at http://www.tbi.univie.ac.at/"
			exit(0)

		
	def checkForpKiss(self):
		"""
			Checking for the presence of pKiss
		"""
		pKiss_output = subprocess.Popen(["which", "pKiss_mfe"], stdout=subprocess.PIPE).communicate()[0].strip()
		#pKiss_output = subprocess.Popen(["which", "/usr/local/pkiss/2014-03-17/bin/pKiss_mfe"], stdout=subprocess.PIPE).communicate()[0].strip()
		if len(pKiss_output) > 0 and pKiss_output.find("found") == -1 and pKiss_output.find(" no ") == -1:
			return True
		else:
			print "It seems that pKiss is not installed on your machine. Please do so!"
			print "You can get it at http://bibiserv2.cebitec.uni-bielefeld.de/pkiss"
			exit(0)

	def checkForIPKnot(self):
		"""
			Checking for the presence of IPKnot
		"""
		pKiss_output = subprocess.Popen(["which", "ipknot"], stdout=subprocess.PIPE).communicate()[0].strip()
		if len(pKiss_output) > 0 and pKiss_output.find("found") == -1 and pKiss_output.find(" no ") == -1:
			return True
		else:
			print "It seems that IPKnot is not installed on your machine. Please do so!"
			print "You can get it at http://rtips.dna.bio.keio.ac.jp/ipknot/"
			exit(0)

	def checkForHotKnots(self):
		"""
			Checking for the presence of HotKnots
		"""
		cmd = self.HotKnots_PATH + "/HotKnots"
		
		#HotKnots_output = subprocess.Popen([cmd], stdout=subprocess.PIPE).communicate()[0].strip()
		#if len(HotKnots_output) > 0 and HotKnots_output.find("found") == -1 and HotKnots_output.find(" no ") == -1:
			#return True
		
		if os.path.exists(cmd):
			return True
		else:
			print "It seems that HotKnots is not installed on your machine. Please do so!"
			print "You can get it at http://www.cs.ubc.ca/labs/beta/Software/HotKnots/"
			exit(0)

#########################
# ENTRY TO THE ANT HIVE
#########################

if __name__ == "__main__":

	exe()
    

  

 





  

 
