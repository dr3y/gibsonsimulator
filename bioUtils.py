from copy import deepcopy
import Levenshtein

def dc(inpt):
    """wrapper for deep copy"""
    return deepcopy(inpt)
def rc(seq):
    """reverse complement of the sequence!!"""
    revdict = {"A":"T",
               "T":"A",
               "G":"C",
               "C":"G",
               "a":"t",
               "t":"a",
               "g":"c",
               "c":"g"}
    try:
        outseq = ''.join([revdict[a] for a in seq][::-1])
    except KeyError:
        outseq = seq
    return outseq
def gcContent(seq):
    """counts number of C and G in the seqence"""
    seq = seq.upper()
    gc  = float(seq.count('G')+seq.count('C'))/len(seq)
    return gc
def allpermut(mutlist):
    """gives a list of lists which contains all permuatations of an input"""
    outlst = []
    if(len(mutlist)>1):
        for a in mutlist:
            modmutlist = dc(mutlist)
            modmutlist.remove(a)
            outlst+= [[a] + b for b in allpermut(modmutlist)]
    else:
        return [mutlist]
    return outlst

def allComb(partslist):
    '''recursively finds all possible paths through the partslist'''
    if(len(partslist)==1):
        return [[a] for a in partslist[0]]
    else:
        result = []
        for a in partslist[0]:
            result+=[[a]+b for b in allComb(partslist[1:])]
        return result
def differences(seq1,seq2,difct = 0):
    """counts the number of differences between seq1 and 2, if they are aligned at the
    left end. Does no internal alignments!! Stops after Difct"""
    diff = 0
    while len(seq1) >0:
        diff += seq1[0]!=seq2[0]
        if(diff >= difct and difct > 0):
            break;
        seq1 = seq1[1:]
        seq2 = seq2[1:]
    return diff
def levenWeight(kmer,match,maxdif=10):
    x = Levenshtein.distance(kmer,match)
    return 50/((x+1)**2)