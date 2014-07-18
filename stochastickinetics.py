#from pylab import figure, plot, xlabel, grid, hold, legend, title, savefig
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import random
from scipy.misc import comb
import math
from copy import deepcopy as cp
from matplotlib.mlab import griddata
import numpy as np
import time


def timestep(X=[],R=[],P={},endlist=[]):
    '''time until the next reaction happens, and also tells me which reaction it was'''
    #X is a list of numbers of ends.
    #R is a list of which ends react with which
    #also it tells us what sequence of fragments the molecules have
    
    if(X==[]):
        X=[100,100,100,100]
        #sum up all the values in a row
    if(R==[]):
        R = [[1,2,.1],[1,4,.1],[1,5,.05],
            [2,1,.1],[2,4,.01],
            ]
        #this lists two ends which can react, and their rate
        #so when a reaction happens we need to remove the proper ends
    if(P=={}):
        P={'001_002r_003':(20,0,4),'000_003_002r':(10,5,3)}
        #each polymer has (amount, leftend, rightend)
    if(endlist ==[]):
        endlist = [[['001_002r_003',0],['000_003_002r',1]],[['001_002r_003',1],['000_003_002r',0]]]
        #[who it belongs to,whether its the right end]
    A=[]
    A_0 = 0
    i=0
    for reaction in R:
        #no fancy combinatorial calcs, just
        #multiply the amount of one end that reacts
        #by the amount of the other end that reacts
        h=X[reaction[0]]*X[reaction[1]]
            
        ApoA = h*reaction[2]
        A+=[ApoA]
        A_0 += ApoA
        i+=1
    
    r1 = random.random()#+.0001
    r2 = random.random()
    #if(A_0 <=0.001):
    #    A_0 == 0.001
    #mlog = math.log(1.0/r1)
    #if(mlog <= 0.001):
    #    mlog = 0.001
    Tau = 1.0/A_0*math.log(1.0/r1)
    reactionchoice = r2*A_0
    ApoMoo = 0.0
    Moo = 0
    for aval in A:
        ApoMoo += aval
        if(ApoMoo >= reactionchoice):
            break
        Moo += 1
    #print 'moo'
    #print Moo
    reaction = R[Moo]
     
    #now, remove both ends that reacted
    #somehow we also need to record that this fragment was made and
    #how many of them exist
    
    #.0...1...2...3...4.    
    #0 1 2 3 4 5 6 7 8 9
    
    #now we subtract from the ends that reacted
    #print P
    '''
    print "here is the reaction {}".format(reaction)
    print "X values {}".format([X[reaction[0]],X[reaction[1]]])
    print "pre-endlist"
    print endlist[reaction[0]]
    print "pre-endlist"
    print endlist[reaction[1]]
    print "pre-counts"
    c1 =[P[a[0]][0] for a in endlist[reaction[0]]]
    c2 = [P[a[0]][0] for a in endlist[reaction[1]]]
    print c1
    print c2
    
    if(sum(c1) != X[reaction[0]]):
        print "ERRRRROOOOORRRRRR"
        time.sleep(10)
    if(sum(c2) != X[reaction[1]]):
        print "ERRRROOOOOOORRRRRR"
        time.sleep(10)
    '''    
    pol1pick = random.randint(0,X[reaction[0]])
    pol2pick = random.randint(0,X[reaction[1]])
    pol1 = ''
    pol2 = ''
    #print X
    #print P
    #print "endlist reaction 1"
    #print endlist[reaction[0]]
    #print "endlist reaction 2"
    #print endlist[reaction[1]]
    
    
    for polymer in endlist[reaction[0]]:
        # print "pol1pick is {}".format(pol1pick)
        polamt = 0
        try:
            polamt = P[polymer[0]][0]
        except KeyError:
            continue
            #if the thing we are looking for does not exist,
            #remove that entry in endlist and also set polamt to zero
            #print endlist[reaction[0]]
            #badind = endlist[reaction[0]].index(polymer)
            #del endlist[reaction[0]][badind]
            
            
        #print "polamt is {}".format(polamt)
        pol1pick-= polamt
        if(pol1pick <= 0):
            pol1 = cp(polymer)
            break;
    
    for polymer in endlist[reaction[1]]:
        #print "pol2pick is {}".format(pol2pick)
        polamt = 0
        try:
            polamt = P[polymer[0]][0]
        except KeyError:
            continue
            #badind = endlist[reaction[0]].index(polymer)
            #del endlist[reaction[0]][badind]
        #print "polamt is {}".format(polamt)
        pol2pick-= polamt
        if(pol2pick <= 0):
            pol2 = cp(polymer)
            break;
    #print "pol1pick {}".format(pol1pick)
    #print "pol2pick {}".format(pol2pick)
   # print "pol1 is {}".format(pol1)
    #print "pol2 is {}".format(pol2)
    #decide which one you need to split
    flipp1 = not pol1[1]
    flipp2 = pol2[1]
    
    p1str = pol1[0]#.split('_')
    p2str = pol2[0]#.split('_')
    
    #construct the new polymer code
    leftend = P[pol1[0]][1]
    rightend = P[pol2[0]][2]
    if(flipp1 and flipp2):
        newpol = '{}_{}'.format(p2str,p1str)
        leftend = P[pol2[0]][1]
        rightend = P[pol1[0]][2]
    else:
        if(flipp1):
            leftend = P[pol1[0]][2]
            p1str = rcPolymer(p1str)
        if(flipp2):
            rightend = P[pol2[0]][1]
            p2str = rcPolymer(p2str)
        
            
        newpol = '{}_{}'.format(p1str,p2str)
    
    #reactants are gone
    
    P[pol1[0]] = (P[pol1[0]][0]-1,P[pol1[0]][1],P[pol1[0]][2])
    if(P[pol1[0]][0]==0): #if there are zero left, delete it
        del P[pol1[0]]
    if(pol2[0] != pol1[0]):
        P[pol2[0]] = (P[pol2[0]][0]-1,P[pol2[0]][1],P[pol2[0]][2])
        if(P[pol2[0]][0]==0):
            del P[pol2[0]]
    #ends are gone
    X[reaction[0]]-=1
    X[reaction[1]]-=1
    #print "reaction!"
    #new polymer is recorded
    reverse = False
    
    npolrc = rcPolymer(newpol)

      
    try:    
        P[newpol] =(P[newpol][0]+ 1,P[newpol][1],P[newpol][2])
        #print "found forwards"
    except KeyError:
        try:
            
            npolrc = rcPolymer(newpol)
            P[npolrc]=(P[npolrc][0] + 1,P[npolrc][1],P[npolrc][2])
            reverse=True
            #print "found RC"
        except KeyError:
            
            P[newpol] = (1,leftend,rightend)
           # print "found nothing"
    if(reverse):
        newpol = rcPolymer(newpol)
        lend = leftend
        rend = rightend
        leftend = rend
        rightend = lend
    
    #new ends are recorded
   # print "appended!!"
    lp = (newpol,0)
    rp = (newpol,1)
    if(not lp in endlist[leftend]):
        endlist[leftend].append(lp)

    if(not rp in endlist[rightend]):
        endlist[rightend].append(rp)
    
    
        
        
    #print X
    '''
    print "after reaction"
    print "X-values {}".format([X[reaction[0]],X[reaction[1]]])
    print "new polymer {} {}".format(newpol, P[newpol])
    print "endlists"
    print endlist[leftend]
    print endlist[rightend]
    print "post-counts"
    c1 =[P[a[0]][0] for a in endlist[reaction[0]]]
    c2 = [P[a[0]][0] for a in endlist[reaction[1]]]
    print c1
    print c2
    
    if(sum(c1) != X[reaction[0]]):
        print "ERRRRROOOOORRRRRR"
        time.sleep(10)
    if(sum(c2) != X[reaction[1]]):
        print "ERRRROOOOOOORRRRRR"
        time.sleep(10)
     '''   
    return Tau,X,endlist,P
#"""    
def rcPolymer(polinpt):
   # print "input {}".format(polinpt)
    polinpt = polinpt.split("_")
    polstr = polinpt[::-1]
    for e in range(len(polstr)):
        #print polstr
        if(polstr[e][-1]=='r'):
            polstr[e] = polstr[e][:-1]
        else:
            polstr[e] = polstr[e]+'r'
    polstr = "_".join(polstr)
    #print polstr
    return polstr

def avgRuns(runscore):
    '''compute the average plot for a lot of different plots'''
    binsize = 0.01
    time = 1.0
    tlist = []
    alist = []
    blist = []
    clist = []
    for run in runscore:
        i=0
        for reaction in range(len(run[0])):
            time = run[0][reaction]
            conca = run[1][reaction]
            concb = run[2][reaction]
            concc = run[3][reaction]
            if(len(tlist)<i+1):
                tlist+=[[time]]
            else:
                tlist[i]+=[time]
            if(len(alist)<i+1):
                alist+=[[conca]]
            else:
                alist[i]+=[conca]
            if(len(blist)<i+1):
                blist+=[[concb]]
            else:
                blist[i]+=[concb]
            if(len(clist)<i+1):
                clist+=[[concc]]
            else:
                clist[i]+=[concc]
            i+=1
    outlist = [[float(sum(a))/len(a) for a in tlist],
               [float(sum(a))/len(a) for a in alist],
               [float(sum(a))/len(a) for a in blist],
               [float(sum(a))/len(a) for a in clist]]
    return outlist

def partC():
    t=[stoptime * float(i) / (numpoints-1.0) for i in range(numpoints)]
    
    wsol = odeint(autocat,w0,t,args=(rateconstants,))#,atol=abserr,rtol=relerr)
    #plt.figure()
    #print wsol
    X=[]
    Y=[]
    Z=[]
    for el in wsol:
        X += [el[0]]
        Y+=[el[1]]
        Z+=[el[2]]
    print len(t)
    print len(X)
    plt.plot(t,X,'--')
    plt.plot(t,Y,'--')
    plt.plot(t,Z,'--')
    plt.legend(['A','B','Z'])
    #plt.show()

def partD(runscore = []):
    time = 1.0 #1 second
    
    X = [1,1,100]
    #X[2] = int(random.random()*20+90)
    
    accX = [cp(X)]
    accT = [0]
    worktime = 0.0
    while worktime < time:
        t,newX = timestep(X)
        accX+=[cp(newX)]
        X = newX
        worktime+=t
        accT+=[worktime]
        
    #print accX
    #print levelslist
    x1 = [a[0] for a in accX]
    x2 = [a[1] for a in accX]
    x3 = [a[2] for a in accX]
    run = [accT,x1,x2,x3]
    runscore+=[run]
    plt.plot(accT,x1)
    plt.plot(accT,x2)
    plt.plot(accT,x3)
    plt.legend(['A','B','Z'])
    return runscore

def partE():
    runscore = []
    for a in range(100):
        runscore = partD(runscore)
    print len(runscore)
    avruns = avgRuns(runscore)
    plt.axis([0.0,1.0,0,100])
    plt.show()
    
    plt.plot(avruns[0],avruns[1],'--')
    plt.plot(avruns[0],avruns[2],'--')
    plt.plot(avruns[0],avruns[3],'--')
    plt.legend(['A','B','Z'])
    return runscore

def partF(runscore):
    finalVals=[[],[]]
    for run in runscore:
        finalVals[0]+=[run[1][-1]]
        finalVals[1]+=[run[2][-1]]
    #z = griddata(
    z = plt.hist2d(finalVals[0],finalVals[1],bins=11)[0]
    plt.show()
    x = [a*10 for a in range(11)]#finalVals[0]
    y = [a*10 for a in range(11)]#finalVals[1]
    print x
    print y
    print z
    plt.contourf(x,y,z)
    plt.colorbar()
    #hist2d(finalVals[0],finalVals[1],bins=25)
if(__name__=="__main__"):
        
    w0 = [1,1,100] #initial conditions of variables
    rateconstants = [0.09,0.09] #rate constants
    
    abserr = 1.0e-8
    relerr = 1.0e-6
    stoptime = 1.0
    numpoints = 250
    
    #runscore = partE()
    
    partC()
    plt.xlabel("time")
    plt.ylabel("Number of molecules")
    
    plt.axis([0.0,1.0,0,100])
    
    plt.show()
    '''
    partF(runscore)
    
    plt.xlabel("A")
    plt.ylabel("B")
    
    plt.axis([0,100,0,100])
    plt.show()'''