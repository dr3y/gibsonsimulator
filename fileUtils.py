import random
import glob
import os
def readPrimers(infile):
    """takes in an output from primer generator program"""
    fle = open(infile,'r')
    primers = []
    for lin in fle:
        slin = lin.split()
        primers += [(slin[1],slin[3])]
    fle.close()
    return primers
def readOligoSeqs(infile,acceptable=[]):
    """reads in oligo sequences from a fasta file"""
    fle = open(infile,'r')
    outseqs = []
    good = 1
    for lin in fle:
        lspl = lin.split("_")
        if(len(lspl)>1 and "C" in lspl and acceptable !=[]):
            conum = int(lspl[lspl.index("C")+1])
            good = conum in acceptable
        if(lin[0] == ">" and good):
            outseqs+=[(lin.strip("\r\n"),'')]
        elif(lin[0] in ["a","g","c","t","A","G","C","T"] and outseqs[-1][1]=='') :
            outseqs[-1] = (outseqs[-1][0]+"_"+str(infile.split("\\")[-1]),lin.strip("\r\n"))
        elif(lin[0] in ["a","g","c","t","A","G","C","T"] and outseqs[-1][1]!='') :
            outseqs[-1] = (outseqs[-1][0],outseqs[-1][1]+lin.strip("\r\n"))
    return outseqs

def writeList(outlist,path = ""):
    """writes a list as a CSV"""
    if(path==""):
        path = guttDNA+"temp\\test"+str(random.choice(range(10000)))+".csv"
    f1 = open(path,"wb")
    [f1.write("{},{}\n".format(a[0],a[1])) for a in outlist]
    f1.close()
    print "wrote file {}".format(path.split("\\")[-1])
def readOligoLists(folder,acceptable = []):
    """reads in a bunch of fasta files that are in a folder"""
    flist =[]
    for infile in glob.glob( os.path.join(folder, '*.fas') ):
        flist += [readOligoSeqs(infile,acceptable)]
    return flist