# ----------------------------------------
# Copyright 2009  Chip Stewart, Boston College
# All rights reserved.
# ----------------------------------------

__author__="chip stewart"
__date__ ="$Oct 29, 2009 1:26:36 PM$"
__version__="207"

import struct

class LibrarySpan:

    def __init__(self,f):
        self.V=0     # 207
        self.contigName=""
        self.setName=""
        self.typeName=""
        self.reclen=0
        self.light=0
        self.ReadGroupCode=[]
        self.MedianFragmentLength=[]
        self.SequencingTechnology=[]
        self.LM=[]
        self.LMlow=[]
        self.LMhigh=[]
        self.tailcut=[]
        self.LR=[]
        self.LRmin=[]
        self.LRmax=[]
        self.NPair=[]
        self.NSingle=[]
        self.NPairRedundant=[]
        self.NSingleRedundant=[]
        self.RGC=[]
        self.NLM=[]
        self.ReadGroupCodeX=[]
        self.ReadGroupID=[]
        self.CenterName=[]
        self.Description=[]
        self.LibraryName=[]
        self.PlatformUnit=[]
        self.SampleName=[]
        self.LFhist=[]
        self.LRhist=[]
        self.a=[];
        self.Na=0;

        # packed byte sizes
        Iz = struct.calcsize("i")
        Hz = struct.calcsize("H")
        Qz = struct.calcsize("Q")
        Dz = struct.calcsize("d")

        fid= open(f,"rb")
        data = fid.read()
        fid.close()

        # file pointer
        p=0
        # version (207)
        self.V = (struct.unpack("i", data[p:(p+Iz)]))[0]
        p+=Iz
        ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
        p+=Iz
        self.contigName   = data[p:(p+ns)].strip('\x00')
        p+=ns
        ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
        p+=Iz
        self.setName   = data[p:(p+ns)].strip('\x00')
        p+=ns
        ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
        p+=Iz
        self.typeName   = data[p:(p+ns)].strip('\x00')
        p+=ns
        self.reclen = (struct.unpack("i", data[p:(p+Iz)]))[0]
        p+=Iz
        self.light = (struct.unpack("Q", data[p:(p+Qz)]))[0]
        p+=Qz
        N = (struct.unpack("i", data[p:(p+Iz)]))[0]
        p+=Iz
        self.N=N

        # fixed record size section
        for n in range(N):
            self.ReadGroupCode.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            self.MedianFragmentLength.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            self.SequencingTechnology.append((struct.unpack("H", data[p:(p+Hz)]))[0])
            p+=Hz
            self.LM.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            self.LMlow.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            self.LMhigh.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            self.tailcut.append((struct.unpack("d", data[p:(p+Dz)]))[0])
            p+=Dz
            self.LR.append((struct.unpack("d", data[p:(p+Dz)]))[0])
            p+=Dz
            self.LRmin.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            self.LRmax.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            self.NPair.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            self.NSingle.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            self.NPairRedundant.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            self.NSingleRedundant.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz

            if(p>len(data)):
               print n

            # variable record size section

        for n in range(N):
            self.ReadGroupCodeX.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            self.NLM.append((struct.unpack("I", data[p:(p+Iz)]))[0])
            p+=Iz
            ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            self.ReadGroupID.append(data[p:(p+ns)].strip('\x00'))
            p+=ns
            ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            self.CenterName.append(data[p:(p+ns)].strip('\x00'))
            p+=ns
            ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            self.Description.append(data[p:(p+ns)].strip('\x00'))
            p+=ns
            ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            self.LibraryName.append(data[p:(p+ns)].strip('\x00'))
            p+=ns
            ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            self.PlatformUnit.append(data[p:(p+ns)].strip('\x00'))
            p+=ns
            ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            self.SampleName.append(data[p:(p+ns)].strip('\x00'))
            p+=ns

            # LFhist fragment length histogram
            hist={}
            ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            hist["title"]= data[p:(p+ns)].strip('\x00')
            p+=ns
            hist["Ntot"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["mean"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["std"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["median"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["Nin"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["Nunder"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["Nover"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            Nbin = (struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            hist["Nbin"]= Nbin

            # unpack histogram contents into x
            fmt = "%dd" % (Nbin*3)
            x=struct.unpack(fmt, data[p:(p+Nbin*3*Dz)])
            p+=Dz*3*Nbin
            hist["xc"]=x[0:(3*Nbin):3]
            hist["n"]=x[1:(3*Nbin):3]
            hist["c"]=x[2:(3*Nbin):3]

            self.LFhist.append(hist)

            # LR histogram
            hist={}
            ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            hist["title"]= data[p:(p+ns)].strip('\x00')
            p+=ns
            hist["Ntot"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["mean"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["std"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["median"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["Nin"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["Nunder"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            hist["Nover"]=(struct.unpack("d", data[p:(p+Dz)]))[0]
            p+=Dz
            Nbin = (struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            hist["Nbin"]= Nbin

            # unpack histogram contents into x
            fmt = "%dd" % (Nbin*3)
            x=struct.unpack(fmt, data[p:(p+Nbin*3*Dz)])
            p+=Dz*3*Nbin
            hist["xc"]=x[0:(3*Nbin):3]
            hist["n"]=x[1:(3*Nbin):3]
            hist["c"]=x[2:(3*Nbin):3]

            self.LRhist.append(hist)

            if(p>len(data)):
               print n

        Na = (struct.unpack("i", data[p:(p+Iz)]))[0]
        p+=Iz
        self.Na=Na

        for n in range(Na):

            if (p>(1+len(data))) :
               print n

#            if (n>34):
#               print n
            a1={}
            ns = (struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            a1["name"]= data[p:(p+ns)].strip('\x00')
            p+=ns
            a1["L"]=(struct.unpack("i", data[p:(p+Iz)]))[0]
            p+=Iz
            a1["use"]=((struct.unpack("c", data[p:(p+1)]))[0])=='\x01'
            p=p+1
#            print n,a1
            self.a.append(a1)


#        print self.ReadGroupID, self.N

def main():

    fname='/Users/stewardg/Projects/1000G/Pilot1/maq/build/NA07000.SLX.maq.SRP000031.2009_08.ERR000836.library.span'
    #fname='/Users/stewardg/tmp/NA19131.SLX.maq.SRP000031.2009_08.SRR005554.library.span'
    #fname='/share4/home/stewardg/Projects/1000G/Pilot1/data/build/NA07000/NA07000.SLX.maq.SRP000031.2009_08.ERR000836.library.span'

    print fname

    L=LibrarySpan(fname)

    print L.setName, L.SampleName[0], L.NPair[0], L.LM[0]

if __name__ == "__main__":
    print "test loadLibrarySpan"
    main()
