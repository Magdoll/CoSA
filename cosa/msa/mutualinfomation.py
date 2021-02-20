import sys
from collections import Counter
import math

################################
def readmsa( filename ):
    
    msaid = []
    msaseq = []

    thisid = None
    thisseq = []

    for line in open(filename):
        if line[0] == ">":
            if thisid != None:
                msaid.append(thisid)
                toadd = "".join(thisseq)
                msaseq.append(toadd)
                print("msa",thisid,len(toadd))
            thisid = line.strip()
            thisseq = []
        else:
            thisseq.append(line.strip()) # .toupper()) mafft is all lowercase

    msaid.append(thisid)
    toadd = "".join(thisseq)
    msaseq.append(toadd)
    print("msa",thisid,len(toadd),file=sys.stderr)
    return( (msaid, msaseq) )
          
################################
def pullcol( msa, col ):
    coldat = []
    for row in range(len(msa)):
        coldat.append(msa[row][col])
    return(coldat)

################################
def printmatrix( mm ):

    numkey=14
    key2ind = {"a":0, "c":1, "g":2, "t":3, "n":4, "-":5, "b":6, "d":7, "s":8, "w":9, "y":10, "h":11, "r":12, "v":13}
    ind2key = {v: k for k,v in key2ind.items()}
    keys = "acgtn-"

    print("\t"+"\t".join([char for char in keys]))
    for ii in range(6):
        print(keys[ii]+"\t"+"\t".join(["%d" % (int(xx)) for xx in mm[ii][:6]]))
            
        
################################
def MI( col1, col2 ):

    numkey=14
    key2ind = {"a":0, "c":1, "g":2, "t":3, "n":4, "-":5, "b":6, "d":7, "s":8, "w":9, "y":10, "h":11, "r":12, "v":13}
    counts = [ [0.01]*numkey for xx in range(numkey)]

    for ii in range(len(col1)):
        ind1 = key2ind[col1[ii]]
        ind2 = key2ind[col2[ii]]
        counts[ind1][ind2] +=1.0

    margrow = []
    for ii in range(numkey):
        margrow.append( sum( [ counts[ii][xx] for xx in range(numkey) ] ))
    Z = sum(margrow)
    margrow = [xx/Z for xx in margrow]

    margcol = []
    for ii in range(numkey):
        margcol.append( sum( [ counts[xx][ii] for xx in range(numkey) ] ))
    margcol = [xx/Z for xx in margcol]

    # mi = -p log(p/q)
    mi = []
    for ii in range(numkey):
        for jj in range(numkey):
            # indep
            q = margrow[ii]*margcol[jj]
            # joint
            p = counts[ii][jj]/Z
            mi.append(p*math.log2(p/q))
    return((sum(mi),counts))

    
################################
def entropy( mylist ):
    eps = 0.01
    cc = Counter( mylist ) # list to dict of counts

    # if ("b" in cc) or ("d" in cc) or ("s" in cc) or ("w" in cc) or ("y" in cc) or ("h" in cc) or ("r" in cc) or ("r" in cc):
    #     # is ambiguous
    #     # find read
    #     ambigidx = -1
    #     for ii in range(len(mylist)):
    #         if mylist[ii] not in "acgt-n": ambigidx=ii
    #     print("ambig", cc, ambigidx)

    Z = float(sum([xx+eps for xx in cc.values()]))
    pp = [(xx+eps)/Z for xx in cc.values()]
    entropy = sum([ -xx*math.log2(xx) for xx in pp])
    return(entropy)
    
if __name__ == "__main__":
          
    (ids,seqs) = readmsa(sys.argv[1])

    #### entropy
    myentropy = []
    for col in range(len(seqs[0])):
        myentropy.append(entropy(pullcol(seqs, col)))

    if True:
        ofp = open("covid.entropy.out","w")
        for col in range(len(seqs[0])):
            print(col, myentropy[col], file=ofp)
        ofp.close()

    #### high entropy
    if True:
        ethresh = 0.5
        highent = []
        for ii in range(len(myentropy)):
            if myentropy[ii]>ethresh:
                highent.append( (ii,myentropy[ii]))
        print("len(highent)", len(highent))

        #### pairwise MI
        ofp = open("covid.mi.out","w")

        for ii in range(len(highent)):
            for jj in range(ii+1,len(highent)):
                col1= highent[ii][0]
                col2= highent[jj][0]

                col1dat = pullcol(seqs, col1)
                col2dat = pullcol(seqs, col2)

                mymi = MI(col1dat, col2dat)

                print(col1,col2,mymi[0],file=ofp)
        ofp.close()

    #### top pairs
    if True:
        for xx in [ ( 397 , 3193 ) , ( 397 , 14583 ) , ( 397 , 23591 ) , ( 3193 , 14583 ) , ( 3193 , 23591 ) , ( 8953 , 28333 ) , ( 14583 , 23591 ) , ( 29070 , 29071 ) , ( 29070 , 29072 ) , ( 29071 , 29072 ) ]:
            col1=xx[0]
            col2=xx[1]
            res = MI(pullcol(seqs, col1),pullcol(seqs, col2))
            print("---- col1",col1,"col2",col2,"MI",res[0])
            printmatrix(res[1])

    # haplotypes at positions of top pairs
    allcols = []
    toppos=[397 , 3193 , 8953 , 14583 , 23591 , 28333 , 29070 , 29071 , 29072]
    for pos in toppos:
        allcols.append(pullcol(seqs,pos))

    hapcount = {}
    for ii in range(len(allcols[0])):
        hap = "-".join([ allcols[xx][ii] for xx in range(len(toppos))])
        if not hap in hapcount:
            hapcount[hap] = 1
        else:
            hapcount[hap] += 1
    ss = sorted(hapcount.items(), key=lambda xx: -xx[1])
    for (k,v) in ss:
        print(k,v)
