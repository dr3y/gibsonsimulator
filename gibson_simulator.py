import Bio
import os
import numpy as np
import random as r
import time
import sys
#print sys.path
sys.path.append("C:\\Users\\Andrey\\SkyDrive\\Guttlab\\Python")
import pydna

import genbankUtils as gbU
import fileUtils as filU
import bioUtils as biU
import stochastickinetics as stoc

def makeOverhangs(seqList,length,digest='5prime',cutoff=0,concentration = 100):
    """makes a list of overhangs produced by gibson assembly
    from a list of dsDNA sequences. Cutoff is for 5' restriction
    overhangs!"""
    overarry = None
    endlist = [[] for a in range(len(seqList)*2)]
    polymerlist = {}
    X = []
    for sequenceNum in range(len(seqList)):
        sequence = seqList[sequenceNum][1]
        #print "sequence now is {}".format(sequence)
        if(cutoff > 0):
            sequence = sequence[cutoff:-cutoff]
        odelseq = sequence[length:-length]
        concentration = int(stoc.random.normalvariate(concentration,concentration*.05))
        
        #calculate a random concentration right here!!
        
        X+=[concentration,concentration] #both fragments at the same concentration!
        
        over1 = {'oseq':sequence[:length], \
                 'midseq':odelseq,\
                 'attheend':0,\
                 'sequencenumber':sequenceNum}
        
        over2 = {'oseq':sequence[-length:],\
                 'midseq':odelseq,\
                 'attheend':1,\
                 'sequencenumber':sequenceNum}
        
        polymername = "{:03}".format(sequenceNum)
        
        polymerlist[polymername]=(concentration,sequenceNum*2,sequenceNum*2+1)
        
        endlist[sequenceNum*2].append((polymername,0))
        endlist[sequenceNum*2+1].append((polymername,1))
        
        if(digest == '5prime'):
            over1['oseq'] = biU.rc(over1['oseq'])
        else:
            over2['oseq'] = biU.rc(over2['oseq'])
        #print "over1 is {}".format(over1['oseq'])
        #print "over2 is {}".format(over2['oseq'])
        if(overarry == None):
            #print "newarray"
            overarry = np.array([over1,over2])
        else:
            #print "append"
            overarry = np.append(overarry,[over1,over2])
        #print overarry
        
    enddatabase = open(gbU.guttDNA+'\\temp\\testdb.fas','wb')
    for end in overarry:
        enddatabase.write(">{}_{}\r\n{}\r\n".format(end['sequencenumber'],end['attheend'],end['oseq']))
    enddatabase.close()
    return overarry,endlist,polymerlist,X

def readBlast(infile,threshold = 25):
    blastfile = open(infile,'rb')
    matchlist = []
    for line in blastfile:
        #print "the line is {}".format(line)
        if( len(line) < 5):
            continue;
        lspl = line.split(',')
        query = lspl[0].split("_")
        subject = lspl[1].split("_")
        
        reactant1 = int(query[0])*2+int(query[1])#"{:03}".format(int(query[0])*2+int(query[1]))
        reactant2 = int(subject[0])*2+int(subject[1])#"{:03}".format(int(subject[0])*2+int(subject[1]))
         
        aligned = lspl[3]
        mismatched = lspl[4]
        matched = int(aligned)-int(mismatched)
        
        score = float(matched)/float(threshold)*0.5
        
        if(matched > threshold):
            matchlist += [(reactant1,reactant2,score)]
    return matchlist
        
def makeBLASTDatabase(dbfilename,makeblastdb = "",outfilename="outdb"):
    """creates a BLAST database using the makeblastdb function"""
    if(makeblastdb==""):
        makeblastdb = "makeblastdb -dbtype {} -out \"{}\" -title {} -in \"{}\""
    os.system(makeblastdb.format('nucl',gbU.guttDNA+"temp\\"+outfilename,outfilename,dbfilename))
    print "made db!"
    
def blastSearch(querypath,outpath=gbU.guttDNA+"temp\\searchout.txt",dbpath=gbU.guttDNA+"temp\\outdb",blastdbcommand = ""):
    """searches a given database with a given query file"""
    if(blastdbcommand == ""):
        blastdbcommand = "blastn -query \"{}\" -db \"{}\" -evalue {} -task {} -strand plus -outfmt 10 -out \"{}\""
    os.system(blastdbcommand.format(querypath,dbpath,0.01,'blastn-short',outpath))
    print "search done!"


def partD(X,R,P,endlist,time = 1.0,makeMovie = True,frames=300,movieout = "C:\\Users\\Andrey\\Desktop\\movie"):
    #time = 1.0 #1 second
    frametime = time/frames
    #X = [1,1,100]
    #X[2] = int(random.random()*20+90)
    
    accX = [biU.dc(X)]
    accT = [0]
    worktime = 0.0
    sizelist = [0]*(len(X))
    stoc.plt.ion()
    fig = stoc.plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    line, = ax.plot(sizelist)
    stoc.plt.xlim([1,len(X)/2])
    stoc.plt.ylim([0,X[0]/2])
    count = 0
    lastframe = 0
    while worktime < time:
        count+=1
        #print 'iteration!'
        t,newX,newendlist,newP = stoc.timestep(X,R,P,endlist)
        #print t
        #accX=newP
        ''
        klist = newP.keys()
        sizelist = [0]*(len(X))
        ''
        mult = 400
        
        if(int(worktime/frametime)>lastframe and makeMovie):
            nc = int(worktime/frametime)
            for polymer in klist:
                sizelist[len(polymer.split('_'))]+=newP[polymer][0]
            accX+=[biU.dc(sizelist)]
            #stoc.plt.clear()
            line.set_ydata(sizelist)
            stoc.plt.draw()
            fname = movieout+'\\_tmp%06d.png'%nc
            #print 'Saving frame', fname
            fig.savefig(fname)
            lastframe = int(worktime/frametime)
            #files.append(fname)

            #stoc.plt.show()
        #'''
        X = newX
        endlist = newendlist
        P = newP
        
        worktime+=t
        accT+=[worktime]
        print worktime/time
    
    #actual    
    #print accX
    #print levelslist
    '''
    x1 = [a[0] for a in accX]
    x2 = [a[1] for a in accX]
    x3 = [a[2] for a in accX]
    run = [accT,x1,x2,x3]
    runscore+=[run]
    stoc.plt.plot(accT,x1)
    stoc.plt.plot(accT,x2)
    stoc.plt.plot(accT,x3)
    stoc.plt.legend(['A','B','Z'])
    '''
    return accT,sizelist,newP


def run140716():
    
    desktopPath = "C:\\Users\\Andrey\\Desktop\\"
    filePath = gbU.guttDNA
    print 'reading seqs'
    inseqs = filU.readOligoSeqs(filePath+'Xist Assembly\\140207SynthXist\\WTassy.fas')#'Assemblies\\malat1\\malat1.fas')
    
    print 'making overhangs'
    overarry, endlist, polymerlist,X = makeOverhangs(inseqs,50,concentration=2000)
    print 'endlist {} '.format(endlist[3])
    print 'polymerlist {} '.format(polymerlist)
    
    endsquery = open(filePath+'temp\\testquery.fas','wb')
    for end in overarry:
        endsquery.write("\n>{}_{}\n{}".format(end['sequencenumber'],end['attheend'],biU.rc(end['oseq'])))
    endsquery.close()
    
    print "making blast database"
    #makeBLASTDatabase(filePath+"temp\\testdb.fas")
    print "performing a query"
    #blastSearch(filePath+"temp\\testquery.fas")
    
    print 'matching ends'
    matchlist = readBlast(filePath+"temp\\searchout.txt")
    
    times,answer,polymers = partD(X,matchlist,polymerlist,endlist,frames=300,time = .5,movieout=desktopPath+"movie5")
    #print polymers
    keypol = polymers.keys()
    pollist = []
    for key in keypol:
        pollist += [(polymers[key][0],key)]
    pollist = sorted(pollist)
    print pollist[-5:]
    #print len(answer)
    #print len(answer[0])
    #stoc.plt.plot(answer)#[times,answer])#imshow(answer[:200])
    #stoc.plt.show()
    
    #print matchlist

run140716()

'''
C:\Program Files\NCBI\blast-2.2.29+\bin>makeblastdb.exe -dbtype nucl -out "C:\Users\Andrey\SkyDrive\Guttlab\DNA\temp\newDB" -title mydatabase -in "C:\Users\Andrey\SkyDrive\Guttlab\DNA\temp\WTassy.fas"
'''