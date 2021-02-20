"""
|**********************************************************************;
* Project           : PacBio Covid Analysis
*
* Author            : mbrown
*
* Date created      : 2020-04-03
*
* Goal              : Take multiple sequence alignment of Covid genomes and compute entropy and mututal information to look for covarying positions
*
|**********************************************************************
"""

import sys
import collections
import math
import itertools

################################
class msaServer:

    ################################
    def __init__(self, msafilename):
        self.key2ind = {"a":0, "c":1, "g":2, "t":3, "n":4, "-":5, "b":6, "d":7, "s":8, "w":9, "y":10, "h":11, "r":12, "v":13, "m":14, "k":15}
        # read in MSA
        (self.msaids,self.msa) = self.readmsa(msafilename)
        # compute transpose lazily
        self.column = {}

        # compute reference to msa index mappings
        self.msaIndexToRefindex()

    ################################
    def readmsa( self, filename ):

        msaid = []
        msaseq = []

        thisid = None
        thisseq = []

        cc = 0

        for line in open(filename):
            if line[0] == ">":
                if thisid != None:
                    msaid.append(thisid.replace(">",""))
                    toadd = "".join(thisseq)
                    msaseq.append(toadd)
                    print("msa",cc, thisid,len(toadd),file=sys.stderr)
                thisid = line.strip()
                thisseq = []
                cc += 1
            else:
                thisseq.append(line.strip()) # .toupper()) mafft is all lowercase

        msaid.append(thisid.replace(">",""))
        toadd = "".join(thisseq)
        msaseq.append(toadd)
        print("msa",cc, thisid,len(toadd),file=sys.stderr)
        return( (msaid, msaseq) )

    ################################
    def msaIndexToRefindex( self ):
        # get the canonical reference to map msaIndexToRefindex

        refid = "NC_045512v2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome"

        if not refid in self.msaids:
            print("msaIndexToRefindex reference not in msa",refid,file=sys.stderr)
            return()

        refidx = self.msaids.index(refid)
        refseq = self.msa[refidx]

        cs = list(itertools.accumulate( [ int(xx!="-") for xx in refseq ] ))

        ofp = open("covid.msa2ref.out","w")
        self.ref2msa = {}
        self.msa2ref = {}
        for msaii in range(len(cs)):
            # msaii is 0-based, refii is 1-based per NCBI
            # annotations. 0based: refii = cs[msaii]-1 # -1 so 1based
            # back to 0based. msa before 0th ref has -1 index
            refii = cs[msaii]
            self.msa2ref[msaii] = refii
            self.ref2msa[refii] = msaii
            print(msaii,refii,file=ofp) 
        ofp.close()

    ################################
    def pullcol( self, col ):
        # if computed return it
        if col in self.column: return(self.column[col])

        coldat = []
        for row in range(len(self.msa)):
            coldat.append(self.msa[row][col])
        self.column[col] = coldat

        return(coldat)

    ################################
    def printmatrix( self, mm ):
        ind2key = {v: k for k,v in self.key2ind.items()}
        keys = "acgtn-"
        print("\t"+"\t".join([char for char in keys]))
        for ii in range(6):
            print(keys[ii]+"\t"+"\t".join(["%d" % (int(xx)) for xx in mm[ii][:6]]))
        
    ################################
    def MI( self, col1, col2 ):

        numkey=len(self.key2ind)
        counts = [ [0.01]*numkey for xx in range(numkey)]

        for ii in range(len(col1)):
            ind1 = self.key2ind[col1[ii]]
            ind2 = self.key2ind[col2[ii]]
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
    def computeMITopEntropy( self, topN ):

        self.entropysort = sorted( zip(range(len(server.myentropy)), server.myentropy),
                              key=lambda xx: xx[1],
                              reverse=True) # sorted returns new, leaving input unchanged
        highent = self.entropysort[:topN]

        print("computeMITopEntropy topN=",topN,"LowestEntropy(pos,entropy)=",highent[-1])

        #### pairwise MI
        self.pairwiseMI = []
        ofp = open("covid.mi.out","w")

        for ii in range(len(highent)):
            for jj in range(ii+1,len(highent)):
                col1idx= highent[ii][0]
                col2idx= highent[jj][0]

                col1dat = server.pullcol(col1idx)
                col2dat = server.pullcol(col2idx)

                mymi = server.MI(col1dat, col2dat)
                
                self.pairwiseMI.append( (col1idx,col2idx,mymi[0] ) )
                print(col1idx,col2idx,mymi[0],file=ofp)
        ofp.close()

    ################################
    def entropy( self, mylist ):
        eps = 0.01
        cc = collections.Counter( mylist ) # list to dict of counts

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

    ################################
    def computeAllEntropy( self ):
        self.myentropy = []
        for col in range(len(self.msa[0])):
            self.myentropy.append(self.entropy(self.pullcol( col)))

        if True:
            ofp = open("covid.entropy.out","w")
            for col in range(len(self.msa[0])):
                print(col, self.myentropy[col], file=ofp)
            ofp.close()

    ################################
    def pair( self, infoline ):
        #### top pairs ( 397 , 3193 ) , ( 397 , 14583 )
        print("================================", infoline)
        (pairstr, col1idx, col2idx) = infoline.split()
        col1idx = int(col1idx)
        col2idx = int(col2idx)
        res = self.MI(self.pullcol(col1idx),self.pullcol(col2idx))
        print("---- col1",col1idx,"col2",col2idx,"MI",res[0])
        self.printmatrix(res[1])
         
    ################################
    def haplotype( self, infoline ):
        #### haplotypes at positions [397 , 3193 , 8953 , 14583 , 23591 , 28333 , 29070 , 29071 , 29072]
        print("================================", infoline)
        (haplotypestr, posdat) = infoline.split(" ",1)
        toppos = [int(xx) for xx in posdat.split()]
        allcols = []
        for pos in toppos:
            allcols.append(self.pullcol(pos))

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

    ################################
    def getmsa2ref( self, infoline ):
        print("================================", infoline)
        (cmdstr, posdat) = infoline.split(" ",1)
        toppos = [int(xx) for xx in posdat.split()]

        print("msa\tref")
        for pos in toppos:
            print(pos, self.msa2ref[pos])

    def getref2msa( self, infoline ):
        print("================================", infoline)
        (cmdstr, posdat) = infoline.split(" ",1)
        toppos = [int(xx) for xx in posdat.split()]

        print("ref\tmsa")
        for pos in toppos:
            print(pos, self.ref2msa[pos])

################################    
if __name__ == "__main__":
          
    server = msaServer( sys.argv[1] )
    topN = int(sys.argv[2])

    #### entropy
    print("# computeAllEntropy")
    server.computeAllEntropy()
    myentropy = server.myentropy

    #### high entropy
    print("# computeMITopntropy")
    server.computeMITopEntropy(  topN ) # look at the top 200 entropy positions for pairwise mutual information

    print("# standard covid.*.out computation DONE! Listening for queries on stdin / 'nc' internet port")

    #### now handle requests for looking at top pairs and haplotypes based on outside analysis (R)
    for line in sys.stdin:

        #### top pairs "pair 397 3193"
        if "pair" in line:
            server.pair( line )

        #### haplotypes at positions "haplotype 397 3193 8953 14583 23591 28333 29070 29071 29072"
        if "haplotype" in line:
            server.haplotype( line )

        ### map ref<>msa "getmsa2ref 419   3215   8978   14608   17947   18058   18260   23616   28358   29095   29096   29097"
        if "getmsa2ref" in line:
            server.getmsa2ref( line )

        if "getref2msa" in line:
            server.getref2msa( line )

        ### quit
        if "quit" in line:
            print("==== Bye!")
            sys.exit(0)
            
