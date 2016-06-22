import subprocess
import os
import argparse
from argparse import RawTextHelpFormatter
import re
import math

def getCorrespondingPositions(positions, n):
    interactions = []
    for entry in positions:
        for i in xrange(1,entry):
            interactions.append((i, entry))
        for j in xrange(entry,n+1):
            interactions.append((entry, j))
    return interactions

def removeUnpairedFrom_bpstack(struct_stack):
    """
        Produce stack for the calculation of accuracy
    """
    tmp_struct_stack = {}
    for i in struct_stack.keys():
        if struct_stack[i] > i:
            tmp_struct_stack[i + 1] = struct_stack[i] + 1
    return tmp_struct_stack

def removePairedAndUndefined_From_bpstack(P, struct_stack):
    """
        Produce stack for the calculation of accessibility
    """
    tmp_struct_stack = {}
    for i in struct_stack.keys():
        if struct_stack[i] == i and P[i] == "x":
            tmp_struct_stack[i + 1] = struct_stack[i] + 1
    return tmp_struct_stack    

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
    
    return (bpstack)

def AccuracyFeature(string):
    """
        argparse type definition of accuracy feature
    """
    m = re.match('^\s*([.()]*)\s*(\w+)\s*(\d+\.?\d+)\s*$', string)
    if m is None:
        print "'" + string + "' is not a valid AccuracyFeature input of type STRUCTURE STRING STRING/FLOAT"
        exit(1)
    s1 = m.group(1)
    s2 = m.group(2)
    s3 = float(m.group(3))
    #print (s1, s2, s3)
    return ( (s1, s2, s3) )
    
def AccessibilityFeature(string):
    """
        argparse type definition of accessibility feature
    """
    m = re.match('^\s*([.x]*)\s*(\w+)\s*(\d+\.?\d+)\s*$', string)
    if m is None:
        print "'" + string + "' is not a valid AccessibilityFeature input of type STRUCTURE STRING STRING/FLOAT"
        exit(1)
    s1 = m.group(1)
    s2 = m.group(2)
    s3 = float(m.group(3))
    if len(s1.replace(".", "")) != len(s1.strip(".")):
        print "defined accessibility", s1, "is not within one stretch! Please reconsider."
        exit(1)

    return ( (s1, s2, s3) )
    
def DiffAccuracyFeature(string):
    """
        argparse type definition of differential accuracy feature
    """
    m = re.match('^\s*([.()]*)\s*(\w+)\s*(\d+\.?\d+)\s*(\w+)\s*(\d+\.?\d+)\s*$', string)
    if m is None:
        print "'" + string + "' is not a valid DiffAccuracyFeature input of type STRUCTURE STRING FLOAT STRING FLOAT"
        exit(1)
    s1 = m.group(1)
    s2 = m.group(2)
    s3 = float(m.group(3))
    s4 = m.group(4)
    s5 = float(m.group(5))
    return ( (s1, s2, s3, s4, s5) )
    
def DiffAccessibilityFeature(string):
    """
        argparse type definition of differential accessibility feature
    """
    m = re.match('^\s*([.x]*)\s*(\w+)\s*(\d+\.?\d+)\s*(\w+)\s*(\d+\.?\d+)\s*$', string)
    if m is None:
        print "'" + string + "' is not a valid DiffAccessibilityFeature input of type STRUCTURE STRING FLOAT STRING FLOAT"
        exit(1)
    s1 = m.group(1)
    s2 = m.group(2)
    s3 = float(m.group(3))
    s4 = m.group(4)
    s5 = float(m.group(5))
    
    if len(s1.replace(".", "")) != len(s1.strip(".")):
        print "defined diff_accessibility", s1, "is not within one stretch! Please reconsider."
        exit(1)
    return ( (s1, s2, s3, s4, s5) )
    
def RNAfold(sequence,Constraint = "", temperature= 37, paramFile = ""):
	
	p = "-p"

	t = "-T " + str(temperature)
	C = ""
	if Constraint != "":
		C = "-C"
	P = ""
	if paramFile != "":
		P = "-P " + paramFile
	p2p = "RNAfold"

	print [p2p, '-d 2', t, P, C, p]
	p = subprocess.Popen( ([p2p, '-d 2', t, P, C, p]),
				stdin = subprocess.PIPE,
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE,
				close_fds = True)
	p.stdin.write(sequence+'\n')
	p.stdin.write(Constraint+'\n\n')
	print p.communicate()
	
	
frame_lines = [
"/frame { %size x y box - draws box centered on x,y\n",
"   2 index 0.5 mul sub            % x -= 0.5\n",
"   exch 2 index 0.5 mul sub exch  % y -= 0.5\n",
"   3 -1 roll dup rectstroke\n",
"} bind def\n",
"\n",
"/uframe {\n",
"   logscale {\n",
"      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n",
"   } if\n",
"   3 1 roll\n",
"   exch len exch sub 1 add frame\n",
"}bind def\n",
"\n",
"/lframe {\n",
   "3 1 roll\n",
"   len exch sub 1 add frame\n",
"} bind def\n"
]



argument_parser = argparse.ArgumentParser(
description = """

	Plot dotplots to the sequence and the given structure. 
	One time as unconstraint UB.ps and the other time the B.ps where the bound constraint is
	plotted
	""",
	formatter_class=RawTextHelpFormatter
	)


argument_parser.add_argument("-Cseq", "--Cseq", 
								help="RNA Sequence\n(TYPE: %(type)s)\n\n", 
								type=str, 
								default = "") 

argument_parser.add_argument("-Cstr", "--Cstr", 
								help="RNA Secondary Constraint Structure.\n(TYPE: %(type)s)\n\n", 
								type=str, 
								default = "") 

argument_parser.add_argument("--accuracy", 
                                help="Define an accuracy evaluation block.\n\n", 
                                type=AccuracyFeature,
                                default=None,
                                action='append'
                                )
argument_parser.add_argument("--accessibility", 
                            help="Define an accessibility evaluation block.\n\n", 
                            type=AccessibilityFeature,
                            default=None,
                            action='append'
                            )

argument_parser.add_argument("--diff-accuracy", 
                            #help="Define an differential accuracy evaluation block.\n\n", 
                            type=DiffAccuracyFeature,
                            default=None,
                            action='append'
                            )

argument_parser.add_argument("--diff-accessibility", 
                            help="Define an differential accessibility evaluation block.\n\n", 
                            type=DiffAccessibilityFeature,
                            default=None,
                            action='append'
                            )

args = argument_parser.parse_args()



sequence  = args.Cseq
structure = args.Cstr




#Folding PART


# UNCONSTRAINED
RNAfold(sequence, "", 37, "")
file_name = "dot.ps"
data = [] 
with open(file_name, 'r') as SOURCE:
	data = SOURCE.readlines()
for i, line in enumerate(data):
	if line.find("dot.ps") != -1:
		data[i] = line.replace("dot.ps", "UNBOUND")

with open(file_name, 'w') as TARGET:
	for i, line  in enumerate(data):
		TARGET.write(line)
new_file_name = "UB.ps"
os.rename(file_name, new_file_name)

if structure :
	# CONSTRAINED
	RNAfold(sequence, structure, 37, "")
	file_name = "dot.ps"
	data = [] 
	with open(file_name, 'r') as SOURCE:
		data = SOURCE.readlines()
	for i, line in enumerate(data):
		if line.find("dot.ps") != -1:
			data[i] = line.replace("dot.ps", "BOUND")

	with open(file_name, 'w') as TARGET:
		for i, line  in enumerate(data):
			TARGET.write(line)
	new_file_name = "B.ps"
	os.rename(file_name, new_file_name)


os.remove("rna.ps")


# PARSING THE CONSTRAINTS
accur = {}
accur["UB"] = []
accur["B"] = []
access = {}
access["UB"] = []
access["B"] = []

if args.accuracy:
    for feature in args.accuracy:
        structure, system, value = feature
        accur[system].append(feature)

if args.accessibility:
    for feature in args.accessibility:
        structure, system, value = feature
        access[system].append(feature)

if args.diff_accuracy:
    print "Yes, indeed" 
    for feature in args.diff_accuracy:
        structure, system1, value1, system2, value2 = feature
        print feature
        accur[system1].append((structure, system1, value1))
        accur[system2].append((structure, system2, value2))

if args.diff_accessibility:

    for feature in args.diff_accessibility:
        structure, system1, value1, system2, value2 = feature 

        access[system1].append((structure, system1, value1))
        access[system2].append((structure, system2, value2))


# HIGHLIGHTING THE CONSTRAINTS
# CONSTRAINT CASE
coloring_B_accur = []
coloring_UB_accur = []
coloring_B_access = []
coloring_UB_access = []

if args.Cstr:
    accur["B"].append((args.Cstr, "B", 1.0))
    s = getbpStack(args.Cstr)
    s = removeUnpairedFrom_bpstack(s)
    for i in s:
        j = s[i]
        coloring_B_accur.append("%s %s %s %s %s %s\n" % ("1 0 0", "setrgbcolor", i, j, 1.0, "lframe") )

if len(accur["B"]) > 0:

    for feature in accur["B"]:
        structure, system, value = feature
        s = getbpStack(structure)
        s = removeUnpairedFrom_bpstack(s)
        print "accur, B", s
        for i in s:
            j = s[i]
            coloring_B_accur.append("%s %s %s %s %s %s\n" % ("0 1 0", "setrgbcolor", i, j, math.sqrt(value), "lframe") )



if len(access["B"]) > 0:
    for feature in access["B"]:
        structure, system, value = feature
        s = getbpStack(structure)
        s = removePairedAndUndefined_From_bpstack(structure, s)
        print "access B" , s
        s = getCorrespondingPositions(s, len(args.Cseq))
        for pair in s:
            i,j = pair
            coloring_B_access.append("%s %s %s %s %s %s\n" % ("0.95 0.95 0.95", "setrgbcolor", i, j, 0.8, "lbox") )

if len(accur["UB"]) > 0:
    for feature in accur["UB"]:
        structure, system, value = feature
        s = getbpStack(structure)
        s = removeUnpairedFrom_bpstack(s)
        print "accur, UB", s
        for i in s:
            j = s[i]
            coloring_UB_accur.append("%s %s %s %s %s %s\n" % ("0 1 0", "setrgbcolor", i, j, math.sqrt(value), "uframe") )

if len(access["UB"]) > 0:
    for feature in access["UB"]:
        structure, system, value = feature
        s = getbpStack(structure)
        s = removePairedAndUndefined_From_bpstack(structure, s)
        print "access UB" , s
        s = getCorrespondingPositions(s, len(args.Cseq))
        for pair in s:
            i,j = pair
            coloring_UB_access.append("%s %s %s %s %s %s\n" % ("0.95 0.95 0.95", "setrgbcolor", i, j, 0.8, "ubox") )
info_B= []
if len(access["B"]) > 0 or len(accur["B"]) > 0:
    b_content = []

    for c_b in coloring_B_access:
        info_B.append(c_b)
    with open("B.ps", 'r') as SOURCE:
        b_content = SOURCE.readlines()
    for i, line in enumerate(b_content):
        if line.strip("\n").endswith("ubox"):
            info_B.append("0 0 0 setrgbcolor " + line.replace("ubox", "lbox"))
    for c_b in coloring_B_accur:
        info_B.append(c_b)


d = []
tail = []
head = []
infoUB = []



with open("UB.ps", 'r') as SOURCE:
    d = SOURCE.readlines()

tail.append(d.pop())
tail.append(d.pop())
tail.append(d.pop())


while(not d[0].startswith("/ubox")):
    head.append(d.pop(0))

for f_line in frame_lines:
    head.append(f_line)
head.append("\n")

for i, line in enumerate(d):
    if not line.strip("\n").endswith("lbox"):
        infoUB.append(line)

dots = []

for i, line  in enumerate(infoUB):
    if line.strip("\n").endswith("ubox"):
        dots.append(line)
for i, line in enumerate(dots):
    del infoUB[infoUB.index(line)]

if info_B:
    for i, line in enumerate(infoUB):
        if line.strip("\n").endswith("(UNBOUND) show"):
                tmp_line = line.replace("UNBOUND", "BOUND")
                tmp_line_split = tmp_line.split(" ")
                tmp_line_split[0] = str(0)
                tmp_line_split[1] = str(int(tmp_line_split[1]) - 100)
                infoUB.insert(i+1, " ".join(tmp_line_split))

infoUB = head + infoUB


if coloring_UB_access:
    for c_ub in coloring_UB_access:
        infoUB.append(c_ub)


for i, line in enumerate(dots):
    infoUB.append("0 0 0 setrgbcolor " + line)

if coloring_UB_accur:
    for c_ub in coloring_UB_accur:
        infoUB.append(c_ub)


if info_B:
    for i, line in enumerate(info_B):
        infoUB.append(line)

while(len(tail) > 0):
    infoUB.append(tail.pop())


with open("UB.ps", 'w') as TARGET:
    for i, line in enumerate(infoUB):
        TARGET.write(line)





# UNCONSTRAINT CASE





# if structure :
# 	os.system("gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=Constraint_System.pdf B.ps UB.ps")
# 	os.remove("B.ps")
# else:
os.system("gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=Constraint_System.pdf UB.ps")
os.remove("UB.ps")


"""
    uc = fold_path + "/" + name + "|uc_dp.ps"
    c = fold_path + "/" + name + "|c_dp.ps"

    blanck_dotplot_file = c

    UC_bps = {}
    C_bps = {}

    sequence_length = 0
    BLANKFILE = "Constraint.ps"
    BLANKFILE_DATA = []


    ### CONSTRAINT and BLANKFILE_MANAGEMENT

    with open(c, 'r') as SOURCE:
        data = SOURCE.readlines()
        BLANKFILE_DATA = [i for i in data if not i.endswith("ubox\n") and not i.endswith("lbox\n")]
        deleto = 0
        for i in xrange(len(BLANKFILE_DATA)):
            if BLANKFILE_DATA[i].find("/Helvetica findfont 14 scalefont setfont") != -1:
                deleto = (i, i-1)
            elif BLANKFILE_DATA[i].startswith("/sequence"):
                sequence_length = len(BLANKFILE_DATA[i+1].strip("\n"))
        a,b = deleto
        del BLANKFILE_DATA[a]
        del BLANKFILE_DATA[b]

        for i in xrange(1, sequence_length + 1):
            for j in xrange(1, sequence_length + 1):
                if j > i:
                    UC_bps[str(i) + "_" + str(j)] = 0
                    C_bps[str(i) + "_" + str(j)] = 0

        datas = [i.strip("\n") for i in data if i.endswith("ubox\n") and not i.startswith("%")]
        for entry in datas:
            es = entry.split(" ")
            C_bps[es[0] + "_" + es[1]] = float(es[2])

        C_acc_sqstr = accuracy(sequestor_hp, C_bps)
        C_acc_aptmr = accuracy(apta_hp, C_bps)


    ### UNCONSTRAINT
    with open(uc, 'r') as SOURCE:
        data = SOURCE.readlines()
    data = [i.strip("\n") for i in data if i.endswith("ubox\n") and not i.startswith("%")]
    for entry in data:
        es = entry.split(" ")
        UC_bps[es[0] + "_" + es[1]] = float(es[2])
    UC_acc_sqstr = accuracy(sequestor_hp, UC_bps)
    UC_acc_aptmr = accuracy(apta_hp, UC_bps)

    acc_diff[C_acc_sqstr - UC_acc_sqstr] =  (name, C_acc_aptmr-UC_acc_aptmr) 
    # before it was just the UC_acc_aptmr which was used as indicative value for the aptamer
    
    DIFF_BP = {}

    for i in xrange(1, sequence_length + 1):
        for j in xrange(1, sequence_length + 1):
            if j > i:
                difference = C_bps[str(i) + "_" + str(j)] - UC_bps[str(i) + "_" + str(j)]
                color = ""
                #if difference > 0:
                    #color = "0 " + str(abs(difference)) + " 0"
                #elif difference < 0:
                    #color = str(abs(difference)) + " 0 0"
                #else:
                    #color = "0 0 " + str(abs(difference))
                if difference > 0:
                    color = "0 1 0"
                elif difference < 0:
                    color = "1 0 0"
                else:
                    color = "0 0 1"
                DIFF_BP[str(i) + "_" + str(j)] = (abs(difference), color)
    SD_bp = {}
    if SD != "":
        for i in xrange(len(SD)):
            if SD[i] == "x":
                for p_i in xrange(1,i + 1):
                    SD_bp[str(i + 1) + "_" + str(p_i)] = 1

                for p_j in xrange(i + 2,len(SD)+1):
                    SD_bp[str(p_j) + "_" + str(i + 1)] = 1
    #print SD_bp
    os.chdir(diff_path)
    with open(BLANKFILE, 'w') as TARGET:
        for i in xrange(len(BLANKFILE_DATA) - 3):
            TARGET.write(BLANKFILE_DATA[i])

        for sq_i in SD_bp.keys():
            i, j = sq_i.split("_")
            TARGET.write("%s %s %s %s %s %s\n" % ("0.9 0.9 0.9", "setrgbcolor", j, i, 1, "ubox") )
        for i in sorted(DIFF_BP.keys()):
            value, color = DIFF_BP[i]
            if value != 0:
                a,b = i.split("_")
                #TARGET.write("%s %s %s %s %s %s\n" % ("0 0 0", "setrgbcolor", a, b, value, "ubox") )
                TARGET.write("%s %s %s %s %s %s\n" % (color, "setrgbcolor", a, b, value, "ubox") )
        for sq_i in apta_hp.keys():
            if apta_hp[sq_i] > sq_i:
                TARGET.write("%s %s %s %s %s %s\n" % ("0 0 1", "setrgbcolor", sq_i + 1, apta_hp[sq_i] + 1 , 1, "lbox") )
        for sq_i in sequestor_hp.keys():
            if sequestor_hp[sq_i] > sq_i:
                TARGET.write("%s %s %s %s %s %s\n" % ("0 1 0", "setrgbcolor", sq_i + 1, sequestor_hp[sq_i] + 1 , 1, "lbox") )
        last = BLANKFILE_DATA.pop()
        prelast = BLANKFILE_DATA.pop()
        preprelast = BLANKFILE_DATA.pop()
        TARGET.write(preprelast)
        TARGET.write(prelast)
        TARGET.write(last)


    os.chdir(fold_path)

os.chdir(path)
"""