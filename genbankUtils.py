from InitCommands import *

def header(acnum = "NR_001463", size = 19817, poly="DNA", shape = "linear"):
    """makes the header of the genbank file"""
    return "LOCUS       {}              {} bp    {}     {}   ROD 05-NOV-2013\r\n".format(acnum,size,poly,shape)

def features(featlist):
    """featlist contains [type,[args]]. So far definitions include exon, misc_feature, and primer_bind
    This is a total kludge, need to fix"""
    outtxt = "FEATURES             Location/Qualifiers\r\n"
    types = {"exon":"     exon            {}..{}\r\n                     /gene=\"{}\"\r\n                     /gene_synonym=\"{}\"\r\n                     /inference=\"{}\"\r\n", \
             "misc_feature":"     misc_feature    {}..{}\r\n                     /note=\"{}\"\r\n", \
             "primer_bind":"     primer_bind     {}..{}\r\n                     /note=\"{}\"\r\n                     /note=\"{}\"\r\n                     /note=\"{}\"\r\n                     /note=\"{}\""
    }
    if featlist != []:
        for feature in featlist:
            #print feature
            try:
                outtxt+=types[feature[0]].format(*feature[1])
            except KeyError:
                pass
    return outtxt

def locus(seq):
    """outputs LOCUS block for genbank file. This contains the sequence"""
    outtxt = "ORIGIN      \r\n"
    seqline = ""
    ind = 1
    col = 0
    while len(seq) > 0:
        seqline += str(ind).rjust(9)
        while (len(seq) >= 10) and (col < 6):
            seqline += " "+seq[:10]
            seq = seq[10:]
            ind += 10
            col += 1
        if  (col < 6) and (len(seq) < 10):
            seqline += " "+seq
            seq = ""
        seqline += "\r\n"
        outtxt+= dc(seqline)
        seqline = ""
        col = 0
    outtxt += "//"
    return outtxt

def gbout(gblist,filename = guttDNA+"\\temp\\temp.gb"):
    """writes a genbank file given the appropriate list. This list has
    [[acnum,size,poly,shape],featlist,seq]"""
    #seq,featlist=[],filename="temp.gb",acnum = "NR_001463", size = 19817, poly="DNA", shape = "linear"):
    [[acnum, size, poly, shape],featlist,seq] = gblist
    outfl = open(filename,'wb')
    outtxt = ""
    outtxt+= header(filename.split('\\')[-1][:-3],size,poly,shape)
    #print outtxt
    outtxt+= features(featlist)
    outtxt+= locus(seq)
    outfl.write(outtxt)
    outfl.close()
    print "wrote file {}".format(filename)
    
def gbin(filename):
    """reads a genbank file and outputs a list of the contents"""
    acnum = ""
    size = 0
    poly = ""
    shape = ""
    mode = 0
    featlist = []
    curfeat = ""
    curfeatlst = []
    seq = ""
    gbopen = open(filename)
    for line in gbopen:
        #print line
        line  = line.replace("Exported File", "Exported_File")
        lsplit = line.split()
        if(lsplit == []):
            lsplit = ["_"]
        if lsplit[0] == 'LOCUS' and mode == 0:
            acnum = lsplit[1]
            size = int(lsplit[2])
            poly = lsplit[4]
            shape = lsplit[5]
        if lsplit[0] == "FEATURES" and mode == 0:
            mode = 2
        #if mode == 2:
            #print lsplit
        if (".." in line) and (len(lsplit) == 2) and (mode == 2):
            featlist+=[[lsplit[0],[]]]
            featlist[-1][1]+=[int(lsplit[1].replace("complement(","").split("..")[0])]
            featlist[-1][1]+=[int(lsplit[1].replace(")","").split("..")[1])]
        if lsplit[0][0] == "/" and mode == 2:
            featlist[-1][1]+=[lsplit[0].split("=")[1].replace("\"","")]
        if lsplit[0] == "ORIGIN":
            mode = 1
        if lsplit[0][-1] == "1" and mode == 1:
            seq+= "".join(lsplit[1:])
    return [[acnum, size, poly, shape],featlist,seq]
def adjustFeatures(featurelocs,features):
    """given a list of features in order, and the list of locations they should be at, it adjusts
    the coordinates of the features so the sizes are correct and outputs the fixed feature list"""
    #print featurelocs, features
    fsizes = [a[1][1]-a[1][0] for a in features]
    adjustedpos = []
    offset = 0
    for feature in features:
        feature[1][0] = featurelocs[0]+offset
        feature[1][1] = featurelocs[0]+fsizes[0]+offset
        adjustedpos +=[feature]
        offset+=fsizes[0]
        featurelocs = featurelocs[1:]
        fsizes = fsizes[1:]
    #print adjustedpos
    return adjustedpos

def updateFlist(ofst, featLIST):
    """offset everything in the feature list with a set value"""
    return [(a[0]-ofst,a[1]-ofst) for a in featLIST]

def infeature(seqGBList,klen,softlen = 0):
    posdict = {}
    for feature in seqGBList[1]:
        start = int(feature[1][0])-klen/2+softlen
        end = int(feature[1][1])+(klen-klen/2)-softlen
        for a in range(start,end+1):
            posdict[a] = True
    return posdict
