import subprocess
import os
import argparse
from argparse import RawTextHelpFormatter


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
args = argument_parser.parse_args()



sequence  = args.Cseq
structure = args.Cstr

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


if structure :
	os.system("gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=Constraint_System.pdf B.ps UB.ps")
	os.remove("B.ps")
else:
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