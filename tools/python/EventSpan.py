# ----------------------------------------
# Copyright 2009  Chip Stewart, Boston College
# All rights reserved.
# ----------------------------------------

__author__ = "chip stewart"
__date__ = "$Nov 4, 2009 1:26:36 PM$"
__version__ = "207"

import struct

class Event1Span():

    def __init__(self, f):
        self.V = 0     # 207
        self.contigName = ""
        self.setName = ""
        self.typeName = ""
        self.ReadGroupCode = []
        self.var = []
        self.p = []
        self.l = []
        self.a = []
        self.q = []
        self.Nr = []
        self.p0 = []
        self.p1 = []
        self.Ns = []
        self.Nm = []
        self.eNr = []
        self.pN = []
        self.Nf = []
        self.pc = []
        self.pcs = []
        self.pcL = []
        self.pcH = []
        self.lc = []
        self.lcs = []
        self.lcL = []
        self.lcH = []
        self.p5a = []
        self.p5b = []
        self.p3a = []
        self.p3b = []
        self.cn = []
        self.posU = []
        self.lenU = []
        self.outlier = []
        self.aberr = []
        self.source = []
        self.merge = []
        self.a5 = []
        self.a3 = []
        self.nr5 = []
        self.nr3 = []
        self.N5 = []
        self.N3 = []
        self.e = {}

        data = open(f).readlines() # open the file

        n = 0;
        for line in data:
            n += 1
            x = line.split()
            if (n == 1):
                self.typeName = x[0]
                self.setName = x[1]
                self.contigName = x[2]
                if (len(x) > 3):
                    self.setName2 = x[4]

            elif (n == 2):
                self.N = int(x[0])
                
            elif (n == 3):
                self.var = x
                for i in range(len(self.var)):
                    self.e[self.var[i]] = []

            else:
                for i in range(len(self.var)):
                    x1 = x[i]
                    self.e[self.var[i]].append(x1)

                self.p.append(int(x[0]))
                self.l.append(int(x[1]))
                self.a.append(int(x[2]))
                self.q.append(int(x[3]))
                self.Nr.append(int(x[4]))
                self.p0.append(int(x[5]))
                self.p1.append(int(x[6]))
                self.Ns.append(int(x[7]))
                self.Nm.append(int(x[8]))
                self.eNr.append(float(x[9]))
                self.pN.append(float(x[10]))
                self.Nf.append(int(float(x[11])))
                self.pc.append(float(x[12]))
                self.pcs.append(float(x[13]))
                self.pcL.append(int(float(x[14])))
                self.pcH.append(int(float(x[15])))
                self.lc.append(float(x[16]))
                self.lcs.append(float(x[17]))
                self.lcL.append(int(float(x[18])))
                self.lcH.append(int(float(x[19])))
                self.p5a.append(int(float(x[20])))
                self.p5b.append(int(float(x[21])))
                self.p3a.append(int(float(x[22])))
                self.p3b.append(int(float(x[23])))
                self.cn.append(int(float(x[24])))
                self.posU.append(float(x[25]))
                self.lenU.append(float(x[26]))
                self.outlier.append(float(x[27]))
                self.aberr.append(float(x[28]))
                self.source.append(x[29])
                self.merge.append(int(float(x[30])))
                self.a5.append(float(x[31]))
                self.a3.append(float(x[32]))
                self.nr5.append(float(x[33]))
                self.nr3.append(float(x[34]))
                self.N5.append(int(float(x[35])))
                self.N3.append(int(float(x[36])))

        print self.setName, self.typeName, self.contigName, self.N

class EventInv1Span():

    def __init__(self, f):
        self.V = 0     # 207
        self.contigName = ""
        self.setName = ""
        self.typeName = ""
        self.ReadGroupCode = []
        self.var = []
        self.p = []
        self.l = []
        self.a = []
        self.q = []

        self.Nf5 = []
        self.p5 = []
        self.sp5 = []
        self.p5L = []
        self.p5H = []
        self.l5 = []
        self.sl5 = []
        self.l5L = []
        self.l5H = []

        self.Nf3 = []
        self.p3 = []
        self.sp3 = []
        self.p3L = []
        self.p3H = []
        self.l3 = []
        self.sl3 = []
        self.l3L = []
        self.l3H = []

        self.Nf = []
        self.Nfpre = []
        self.Nfpost = []
        self.eNf = []

        self.p5a = []
        self.p5b = []
        self.p3a = []
        self.p3b = []

        self.Nr = []
        self.p0 = []
        self.p1 = []
        self.Ns = []
        self.Nm = []
        self.eNr = []
        self.cn = []

        self.posU = []
        self.lenU = []
        self.outlier = []
        self.aberr = []
        self.src5 = []
        self.src3 = []
        self.merge5 = []
        self.merge3 = []
        self.a5 = []
        self.a3 = []
        self.nr5 = []
        self.nr3 = []
        self.N5 = []
        self.N3 = []
        self.e = {}

        data = open(f).readlines() # open the file

        n = 0;
        for line in data:
            n += 1
            x = line.split()
            if (n == 1):
                self.typeName = x[0]
                self.setName = x[1]
                self.contigName = x[2]
                if (len(x) > 3):
                    self.setName2 = x[4]

            elif (n == 2):
                self.N = int(x[0])

            elif (n == 3):
                self.var = x
                for i in range(len(self.var)):
                    self.e[self.var[i]] = []

            else:
                for i in range(len(self.var)):
                    x1 = x[i]
                    self.e[self.var[i]].append(x1)

                self.p.append(int(x[0]))
                self.l.append(int(x[1]))
                self.a.append(int(x[2]))
                self.q.append(int(x[3]))

                self.Nf5.append(int(x[4]))
                self.p5.append(int(x[5]))
                self.sp5.append(float(x[6]))
                self.p5L.append(int(x[7]))
                self.p5H.append(int(x[8]))
                self.l5.append(int(x[9]))
                self.sl5.append(float(x[10]))
                self.l5L.append(int(x[11]))
                self.l5H.append(int(x[12]))

                self.Nf3.append(int(x[13]))
                self.p3.append(int(x[14]))
                self.sp3.append(float(x[15]))
                self.p3L.append(int(x[16]))
                self.p3H.append(int(x[17]))
                self.l3.append(int(x[18]))
                self.sl3.append(float(x[19]))
                self.l3L.append(int(x[20]))
                self.l3H.append(int(x[21]))

                self.Nf.append(int(x[22]))
                self.Nfpre.append(int(x[23]))
                self.Nfpost.append(int(x[24]))
                self.eNf.append(float(x[25]))
                self.p5a.append(int(float(x[26])))
                self.p5b.append(int(float(x[27])))
                self.p3a.append(int(float(x[28])))
                self.p3b.append(int(float(x[29])))

                self.Nr.append(int(x[30]))
                self.p0.append(int(x[31]))
                self.p1.append(int(x[32]))
                self.Ns.append(int(x[33]))
                self.Nm.append(int(x[34]))
                self.eNr.append(float(x[35]))
                self.cn.append(int(float(x[36])))

                self.posU.append(float(x[37]))
                self.lenU.append(float(x[38]))
                self.outlier.append(float(x[39]))
                self.aberr.append(float(x[40]))

                self.src5.append(x[41])
                self.src3.append(x[42])
                self.merge5.append(int(float(x[43])))
                self.merge3.append(int(float(x[44])))
                self.a5.append(float(x[45]))
                self.a3.append(float(x[46]))
                self.nr5.append(float(x[47]))
                self.nr3.append(float(x[48]))
                self.N5.append(int(float(x[49])))
                self.N3.append(int(float(x[50])))

        print self.setName, self.typeName, self.contigName, self.N

class EventMob1Span():

    def __init__(self, f):
        self.V = 0     # 207
        self.contigName = ""
        self.setName = ""
        self.typeName = ""
        self.ReadGroupCode = []
        self.var = []
        self.p = []
        self.l = []
        self.a = []
        self.q = []

        self.Nf5 = []
        self.p5 = []
        self.sp5 = []
        self.p5L = []
        self.p5H = []

        self.Nf3 = []
        self.p3 = []
        self.sp3 = []
        self.p3L = []
        self.p3H = []

        self.Nfuu = []
        self.Nf0 = []
        self.Nf1 = []
        self.NfE = []

        self.p5a = []
        self.p5b = []
        self.p3a = []
        self.p3b = []

        self.Nr = []
        self.p0 = []
        self.p1 = []
        self.Ns = []
        self.Nm = []
        self.eNr = []
        self.cn = []

        self.pmed = []
        self.posU = []
        self.lenU = []
        self.src5 = []
        self.src3 = []
        self.merge5 = []
        self.merge3 = []
        self.a5 = []
        self.a3 = []
        self.nr5 = []
        self.nr3 = []
        self.N5 = []
        self.N3 = []
        self.e = {}

        data = open(f).readlines() # open the file

        n = 0;
        for line in data:
            n += 1
            x = line.split()
            if (n == 1):
                self.typeName = x[0]
                self.setName = x[1]
                self.contigName = x[2]
                if (len(x) > 3):
                    self.setName2 = x[4]

            elif (n == 2):
                self.N = int(x[0])

            elif (n == 3):
                self.var = x
                for i in range(len(self.var)):
                    self.e[self.var[i]] = []

            else:
                for i in range(len(self.var)):
                    x1 = x[i]
                    self.e[self.var[i]].append(x1)

                self.p.append(int(x[0]))
                self.l.append(int(x[1]))
                self.a.append(int(x[2]))
                self.q.append(int(x[3]))

                self.Nf5.append(int(x[4]))
                self.p5.append(int(x[5]))
                self.sp5.append(float(x[6]))
                self.p5L.append(int(x[7]))
                self.p5H.append(int(x[8]))


                self.Nf3.append(int(x[9]))
                self.p3.append(int(x[10]))
                self.sp3.append(float(x[11]))
                self.p3L.append(int(x[12]))
                self.p3H.append(int(x[13]))

                self.Nfuu.append(int(x[14]))
                self.Nf0.append(int(x[15]))
                self.Nf1.append(int(x[16]))
                self.NfE.append(int(x[17]))
                self.p5a.append(int(float(x[18])))
                self.p5b.append(int(float(x[19])))
                self.p3a.append(int(float(x[20])))
                self.p3b.append(int(float(x[21])))

                self.Nr.append(int(x[22]))
                self.p0.append(int(x[23]))
                self.p1.append(int(x[24]))
                self.Ns.append(int(x[25]))
                self.Nm.append(int(x[26]))
                self.eNr.append(float(x[27]))
                self.cn.append(int(float(x[28])))

                #S.Ncon5=x{31};
                #S.Ncon3=x{32};

                self.posU.append(float(x[32]))
                self.lenU.append(float(x[33]))

                self.src5.append(x[34])
                self.src3.append(x[35])
                self.merge5.append(int(float(x[36])))
                self.merge3.append(int(float(x[37])))
                self.a5.append(float(x[38]))
                self.a3.append(float(x[39]))

        print self.setName, self.typeName, self.contigName, self.N

class EventCrx1Span():

    def __init__(self, f):
        self.V = 0     # 207
        self.contigName = ""
        self.setName = ""
        self.typeName = ""
        self.ReadGroupCode = []
        self.var = []
        self.p = []
        self.l = []
        self.a = []
        self.q = []

        self.p2 = []
        self.l2 = []
        self.a2 = []

        self.Nf5 = []
        self.p5 = []
        self.sp5 = []
        self.p5L = []
        self.p5H = []
        self.p5a2 = []
        self.sp5a2 = []
        self.p5La2 = []
        self.p5Ha2 = []

        self.Nf3 = []
        self.p3 = []
        self.sp3 = []
        self.p3L = []
        self.p3H = []
        self.p3a2 = []
        self.sp3a2 = []
        self.p3La2 = []
        self.p3Ha2 = []

        self.Nf = []
        self.Nfpre = []
        self.Nfpost = []
        self.eNf = []

        self.p5a = []
        self.p5b = []
        self.p3a = []
        self.p3b = []

        self.p5aT = []
        self.p5bT = []
        self.p3aT = []
        self.p3bT = []

        self.Nr = []
        self.p0 = []
        self.p1 = []
        self.Ns = []
        self.Nm = []
        self.eNr = []
        self.cn = []

        self.posU = []
        self.lenU = []
        self.outlier = []
        self.src5 = []
        self.src3 = []
        self.merge5 = []
        self.merge3 = []
        self.a5 = []
        self.a3 = []
        self.nr5 = []
        self.nr3 = []
        self.N5 = []
        self.N3 = []
        self.e = {}

        data = open(f).readlines() # open the file

        n = 0;
        for line in data:
            n += 1
            x = line.split()
            if (n == 1):
                self.typeName = x[0]
                self.setName = x[1]
                self.contigName = x[2]
                if (len(x) > 3):
                    self.setName2 = x[4]

            elif (n == 2):
                self.N = int(x[0])

            elif (n == 3):
                self.var = x
                for i in range(len(self.var)):
                    self.e[self.var[i]] = []

            else:
                for i in range(len(self.var)):
                    x1 = x[i]
                    self.e[self.var[i]].append(x1)

                self.p.append(int(x[0]))
                self.l.append(int(x[1]))
                self.a.append(int(x[2]))
                self.q.append(int(x[3]))

                self.p2.append(int(x[4]))
                self.l2.append(int(x[5]))
                self.a2.append(int(x[6]))


                self.Nf5.append(int(x[7]))
                self.p5.append(int(x[8]))
                self.sp5.append(float(x[9]))
                self.p5L.append(int(x[10]))
                self.p5H.append(int(x[11]))
                self.p5a2.append(int(x[12]))
                self.sp5a2.append(float(x[13]))
                self.p5La2.append(int(x[14]))
                self.p5Ha2.append(int(x[15]))

                self.Nf3.append(int(x[16]))
                self.p3.append(int(x[17]))
                self.sp3.append(float(x[18]))
                self.p3L.append(int(x[19]))
                self.p3H.append(int(x[20]))
                self.p3a2.append(int(x[21]))
                self.sp3a2.append(float(x[22]))
                self.p3La2.append(int(x[23]))
                self.p3Ha2.append(int(x[24]))

                self.Nf.append(int(x[25]))
                self.Nfpre.append(int(x[26]))
                self.Nfpost.append(int(x[27]))
                self.eNf.append(float(x[28]))
                self.p5a.append(int(float(x[29])))
                self.p5b.append(int(float(x[30])))
                self.p3a.append(int(float(x[31])))
                self.p3b.append(int(float(x[32])))
                self.p5aT.append(int(float(x[33])))
                self.p5bT.append(int(float(x[34])))
                self.p3aT.append(int(float(x[35])))
                self.p3bT.append(int(float(x[36])))

                self.Nr.append(int(x[37]))
                self.p0.append(int(x[38]))
                self.p1.append(int(x[39]))
                self.Ns.append(int(x[40]))
                self.Nm.append(int(x[41]))
                self.eNr.append(float(x[42]))
                self.cn.append(int(float(x[43])))

                self.posU.append(float(x[44]))
                self.lenU.append(float(x[45]))
                self.outlier.append(float(x[46]))

                self.src5.append(x[47])
                self.src3.append(x[48])
                self.merge5.append(int(float(x[49])))
                self.merge3.append(int(float(x[50])))
                self.a5.append(float(x[51]))
                self.a3.append(float(x[52]))
                self.nr5.append(float(x[53]))
                self.nr3.append(float(x[54]))
                self.N5.append(int(float(x[55])))
                self.N3.append(int(float(x[56])))

        print self.setName, self.typeName, self.contigName, self.N


def main():

    #fname = '/Users/stewardg/Projects/1000G/Pilot1/Spanner/detect/test/20.del.span.txt'
    fname= '/Users/stewardg/Projects/CSHL09/Spanner/detectB/LNCAP/1.cross.span.txt'
    print fname

    #E = Event1Span(fname)
    E = EventCrx1Span(fname)

    print E.typeName, E.setName, E.contigName, E.N
    for i in range(E.N):
       print "chr%d\t%d\t%d\t%s" % ( E.a[i],E.p0[i],E.p1[i],E.typeName)
       #print c, E.p[i], E.p[i]+E.l[i],E

if __name__ == "__main__":
    print "test EventSpan"
    main()
