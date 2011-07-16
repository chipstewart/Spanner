#! /usr/bin/python

# ----------------------------------------
# Copyright 2009  Chip Stewart, Boston College
# All rights reserved.
# ----------------------------------------

#import sys
import os
import os.path
import fnmatch
import optparse
import EventSpan

first=1

# inititalize event list
E=[]

__author__="chip stewart"
__date__ ="$Nov 5, 2009 5:02:13 PM$"

def mergeEvt(E1,E2):
        E3=E1;
        E3.p=E1.p+E2.p
        E3.l=E1.l+E2.l
        E3.a=E1.a+E2.a
        E3.q=E1.q+E2.q
        E3.Nr=E1.Nr+E2.Nr
        E3.p0=E1.p0+E2.p0
        E3.p1=E1.p1+E2.p1
        E3.Ns=E1.Ns+E2.Ns
        E3.Nm=E1.Nm+E2.Nm
        E3.eNr=E1.Nr+E2.eNr
        E3.pN=E1.pN+E2.pN
        E3.Nf=E1.Nf+E2.Nf
        E3.pc=E1.pc+E2.pc
        E3.pcs=E1.pcs+E2.pcs
        E3.pcL=E1.pcL+E2.pcL
        E3.pcH=E1.pcH+E2.pcH
        E3.lc=E1.lc+E2.lc
        E3.lcs=E1.lcs+E2.lcs
        E3.lcL=E1.lcL+E2.lcL
        E3.lcH=E1.lcH+E2.lcH
        E3.p5a=E1.p5a+E2.p5a
        E3.p5b=E1.p5b+E2.p5b
        E3.p3a=E1.p3a+E2.p3a
        E3.p3b=E1.p3b+E2.p3b
        E3.cn=E1.cn+E2.cn
        E3.posU=E1.posU+E2.posU
        E3.lenU=E1.lenU+E2.lenU
        E3.outlier=E1.outlier+E2.outlier
        E3.aberr=E1.aberr+E2.aberr
        E3.source=E1.source+E2.source
        E3.merge=E1.merge+E2.merge
        E3.a5=E1.a5+E2.a5
        E3.a3=E1.a3+E2.a3
        E3.nr5=E1.nr5+E2.nr5
        E3.nr3=E1.nr3+E2.nr3
        E3.N5=E1.N5+E2.N5
        E3.N3=E1.N3+E2.N3
        E3.contigName=""
        E3.N+=E2.N
        return E3

def mergeInvEvt(E1,E2):
        E3=E1;
        E3.p=E1.p+E2.p
        E3.l=E1.l+E2.l
        E3.a=E1.a+E2.a
        E3.q=E1.q+E2.q
        E3.Nr=E1.Nr+E2.Nr
        E3.p0=E1.p0+E2.p0
        E3.p1=E1.p1+E2.p1
        E3.Ns=E1.Ns+E2.Ns
        E3.Nm=E1.Nm+E2.Nm
        E3.eNr=E1.Nr+E2.eNr
        E3.Nf=E1.Nf+E2.Nf
        E3.Nf5=E1.Nf5+E2.Nf5
        E3.Nf3=E1.Nf3+E2.Nf3

        # missing some stuff..

        E3.p5a=E1.p5a+E2.p5a
        E3.p5b=E1.p5b+E2.p5b
        E3.p3a=E1.p3a+E2.p3a
        E3.p3b=E1.p3b+E2.p3b
        E3.cn=E1.cn+E2.cn
        E3.posU=E1.posU+E2.posU
        E3.lenU=E1.lenU+E2.lenU
        E3.outlier=E1.outlier+E2.outlier
        E3.aberr=E1.aberr+E2.aberr
        E3.src5=E1.src5+E2.src5
        E3.src3=E1.src3+E2.src3
        E3.merge5=E1.merge5+E2.merge5
        E3.merge3=E1.merge3+E2.merge3
        E3.a5=E1.a5+E2.a5
        E3.a3=E1.a3+E2.a3
        E3.nr5=E1.nr5+E2.nr5
        E3.nr3=E1.nr3+E2.nr3
        E3.N5=E1.N5+E2.N5
        E3.N3=E1.N3+E2.N3
        E3.contigName=""
        E3.N+=E2.N
        return E3

def mergeMobEvt(E1,E2):
        E3=E1;
        E3.p=E1.p+E2.p
        E3.l=E1.l+E2.l
        E3.a=E1.a+E2.a
        E3.q=E1.q+E2.q

        E3.Nf5=E1.Nf5+E2.Nf5
        E3.p5=E1.p5+E2.p5
        E3.sp5=E1.sp5+E2.sp5
        E3.p5L=E1.p5L+E2.p5L
        E3.p5H=E1.p5H+E2.p5H
        E3.Nf3=E1.Nf3+E2.Nf3
        E3.p3=E1.p3+E2.p3
        E3.sp3=E1.sp3+E2.sp3
        E3.p3L=E1.p3L+E2.p3L
        E3.p3H=E1.p3H+E2.p3H

        E3.Nfuu=E1.Nfuu+E2.Nfuu
        E3.Nf0=E1.Nf0+E2.Nf0
        E3.Nf1=E1.Nf1+E2.Nf1
        E3.NfE=E1.NfE+E2.NfE

        E3.p5a=E1.p5a+E2.p5a
        E3.p5b=E1.p5b+E2.p5b
        E3.p3a=E1.p3a+E2.p3a
        E3.p3b=E1.p3b+E2.p3b

        E3.Nr=E1.Nr+E2.Nr
        E3.p0=E1.p0+E2.p0
        E3.p1=E1.p1+E2.p1
        E3.Ns=E1.Ns+E2.Ns
        E3.Nm=E1.Nm+E2.Nm
        E3.eNr=E1.Nr+E2.eNr
        E3.cn=E1.cn+E2.cn

        # missing some stuff..

        E3.posU=E1.posU+E2.posU
        E3.lenU=E1.lenU+E2.lenU
        E3.src5=E1.src5+E2.src5
        E3.src3=E1.src3+E2.src3
        E3.merge5=E1.merge5+E2.merge5
        E3.merge3=E1.merge3+E2.merge3
        E3.a5=E1.a5+E2.a5
        E3.a3=E1.a3+E2.a3
        E3.contigName=""
        E3.N+=E2.N

        return E3

def mergeCrxEvt(E1,E2):
        E3=E1;
        E3.p=E1.p+E2.p
        E3.l=E1.l+E2.l
        E3.a=E1.a+E2.a
        E3.q=E1.q+E2.q
        E3.p2=E1.p2+E2.p2
        E3.l2=E1.l2+E2.l2
        E3.a2=E1.a2+E2.a2
        E3.Nr=E1.Nr+E2.Nr
        E3.p0=E1.p0+E2.p0
        E3.p1=E1.p1+E2.p1
        E3.Ns=E1.Ns+E2.Ns
        E3.Nm=E1.Nm+E2.Nm
        E3.eNr=E1.Nr+E2.eNr
        E3.Nf=E1.Nf+E2.Nf

        E3.Nf5=E1.Nf5+E2.Nf5
        E3.Nf3=E1.Nf3+E2.Nf3

        # missing some stuff..

        E3.p5a=E1.p5a+E2.p5a
        E3.p5b=E1.p5b+E2.p5b
        E3.p3a=E1.p3a+E2.p3a
        E3.p3b=E1.p3b+E2.p3b
        E3.p5aT=E1.p5aT+E2.p5aT
        E3.p5bT=E1.p5bT+E2.p5bT
        E3.p3aT=E1.p3aT+E2.p3aT
        E3.p3bT=E1.p3bT+E2.p3bT
        E3.cn=E1.cn+E2.cn
        E3.posU=E1.posU+E2.posU
        E3.lenU=E1.lenU+E2.lenU
        E3.outlier=E1.outlier+E2.outlier
        E3.src5=E1.src5+E2.src5
        E3.src3=E1.src3+E2.src3
        E3.merge5=E1.merge5+E2.merge5
        E3.merge3=E1.merge3+E2.merge3
        E3.a5=E1.a5+E2.a5
        E3.a3=E1.a3+E2.a3
        E3.nr5=E1.nr5+E2.nr5
        E3.nr3=E1.nr3+E2.nr3
        E3.N5=E1.N5+E2.N5
        E3.N3=E1.N3+E2.N3
        E3.contigName=""

        E3.N+=E2.N


        return E3


def printBed(params,E1):

    global E

    EventType=params["EventType"]

    outbed=params["outbed"].rstrip()

    f= open(outbed, 'w')

    outpath, outname=os.path.split(outbed)
    setName, extension = os.path.splitext(outname)


    print >> f, "track name=%s description=%s" % (setName, E.typeName)
 
    nbrk=2
    if ((EventType=='del')|(EventType=='dup')):
      nbrk=1

    for n in range(E.N):
        a='%d' % E.a[n]
        if (a=='23'):
            a='X'
        if (a=='24'):
            a='Y'

        NF=0
        if (nbrk==1):
            NF=E.Nf[n]
        else:
            NF=E.Nf3[n]+E.Nf5[n]


        print >> f,  "chr%s\t %d\t %d\t %s.L%d.cn%d %d" % (a,E.p0[n],E.p1[n],E.typeName[0:3],E.l[n],E.cn[n],NF)

    f.close()


def printCrossBed(params,E1):

    global E

    EventType=params["EventType"]

    outbed=params["outbed"].rstrip()

    f= open(outbed, 'w')

    outpath, outname=os.path.split(outbed)
    setName, extension = os.path.splitext(outname)


    print >> f, "track name=%s description=%s" % (setName, E.typeName)


    for n in range(E.N):
        a='%d' % E.a[n]
        if (a=='23'):
            a='X'
        elif (a=='24'):
            a='Y'

        a2='%d' % E.a2[n]
        if (a2=='23'):
            a2='X'
        elif (a2=='24'):
            a2='Y'
        elif (a2=='25'):
            a2='MT'



        NF=E.Nf3[n]+E.Nf5[n]

        s='b'
        if (E.Nf3[n]==0):
            s="+"
        elif (E.Nf5[n]==0):
            s="-"

        p1=E.p[n]+100+abs(E.l[n])

        print >> f,  "chr%s\t %d\t %d\t cross:chr%s:%d:%s:%d:%d %d" % (a,E.p[n],p1,a2,E.p2[n],s,E.l[n],E.cn[n],NF)

    f.close()

def filt(xx,ok):
    xok=[]
    for i in ok:
       xok.append(xx[i])
    return xok

def goodEvent1(params,E1):

    Nfragmin=params["Nfragmin"]
    cnslosh=params["cnslosh"]
    EventType=params["EventType"]

    delv= (EventType=='del')
    dupv=(EventType=='dup')

    if ((delv+dupv)<1):
         return E1

    CN=2
    if (delv):
       CN=1+cnslosh
    if (dupv):
       CN=3-cnslosh

    E2=E1

    ok=[]
    for e in range(E1.N):
        if ((delv)&(E1.cn[e]<=CN)):
                ok.append(int(e))
        if ((dupv)&(E1.cn[e]>=CN)):
                ok.append(int(e))

#    xxx=[0,1,2,3,4,5]
#    ok=[0,2,5]
#    xx=filt(xxx,ok)

    E2.p=filt(E1.p,ok)
    E2.l=filt(E1.l,ok)
    E2.a=filt(E1.a,ok)
    E2.q=filt(E1.q,ok)
    E2.Nr=filt(E1.Nr,ok)
    E2.p0=filt(E1.p0,ok)
    E2.p1=filt(E1.p1,ok)
    E2.Ns=filt(E1.Ns,ok)
    E2.Nm=filt(E1.Nm,ok)
    E2.eNr=filt(E1.eNr,ok)
    E2.pN=filt(E1.pN,ok)
    E2.Nf=filt(E1.Nf,ok)
    E2.pc=filt(E1.pc,ok)
    E2.pcs=filt(E1.pcs,ok)
    E2.pcL=filt(E1.pcL,ok)
    E2.pcH=filt(E1.pcH,ok)
    E2.lc=filt(E1.lc,ok)
    E2.lcs=filt(E1.lcs,ok)
    E2.lcL=filt(E1.lcL,ok)
    E2.lcH=filt(E1.lcH,ok)
    E2.p5a=filt(E1.p5a,ok)
    E2.p5b=filt(E1.p5b,ok)
    E2.p3a=filt(E1.p3a,ok)
    E2.p3b=filt(E1.p3b,ok)
    E2.cn=filt(E1.cn,ok)
    E2.posU=filt(E1.posU,ok)
    E2.lenU=filt(E1.lenU,ok)
    E2.outlier=filt(E1.outlier,ok)
    E2.aberr=filt(E1.aberr,ok)
    E2.source=filt(E1.source,ok)
    E2.merge=filt(E1.merge,ok)
    E2.a5=filt(E1.a5,ok)
    E2.a3=filt(E1.a3,ok)
    E2.nr5=filt(E1.nr5,ok)
    E2.nr3=filt(E1.nr3,ok)
    E2.N5=filt(E1.N5,ok)
    E2.N3=filt(E1.N3,ok)
    E2.N=len(E2.a)
    return E2

def goodEventInv(params,E1):

    bothsides=params["bothsides"]
    EventType=params["EventType"]

    inv= (EventType=='inv')

    if (inv<1):
         return E1

    if (bothsides<1):
         return E1


    E2=E1

    ok=[]

    for e in range(E1.N):

        if ((E1.Nf5[e]>0)&(E1.Nf3[e]>0)):
                ok.append(int(e))

#    xxx=[0,1,2,3,4,5]
#    ok=[0,2,5]
#    xx=filt(xxx,ok)

    E2.p=filt(E1.p,ok)
    E2.l=filt(E1.l,ok)
    E2.a=filt(E1.a,ok)
    E2.q=filt(E1.q,ok)
    E2.Nr=filt(E1.Nr,ok)
    E2.p0=filt(E1.p0,ok)
    E2.p1=filt(E1.p1,ok)
    E2.Ns=filt(E1.Ns,ok)
    E2.Nm=filt(E1.Nm,ok)
    E2.eNr=filt(E1.eNr,ok)
    E2.Nf5=filt(E1.Nf5,ok)
    E2.Nf3=filt(E1.Nf3,ok)
    E2.p5a=filt(E1.p5a,ok)
    E2.p5b=filt(E1.p5b,ok)
    E2.p3a=filt(E1.p3a,ok)
    E2.p3b=filt(E1.p3b,ok)
    E2.cn=filt(E1.cn,ok)
    E2.posU=filt(E1.posU,ok)
    E2.lenU=filt(E1.lenU,ok)
    E2.src5=filt(E1.src5,ok)
    E2.src3=filt(E1.src3,ok)
    E2.merge5=filt(E1.merge5,ok)
    E2.merge3=filt(E1.merge3,ok)
    E2.N=len(E2.a)
    return E2

def goodEventMob(params,E1):

    bothsides=params["bothsides"]
    EventType=params["EventType"]

    mob=((EventType=='ALU')|(EventType=='L1')|(EventType=='SVA'))

    if (mob<1):
         return E1

    if (bothsides<1):
         return E1


    E2=E1

    ok=[]
    for e in range(E1.N):
        if ((E1.Nf5[e]>0)&(E1.Nf3[e]>0)):
                ok.append(int(e))

#    xxx=[0,1,2,3,4,5]
#    ok=[0,2,5]
#    xx=filt(xxx,ok)

    E2.p=filt(E1.p,ok)
    E2.l=filt(E1.l,ok)
    E2.a=filt(E1.a,ok)
    E2.q=filt(E1.q,ok)
    E2.Nr=filt(E1.Nr,ok)
    E2.p0=filt(E1.p0,ok)
    E2.p1=filt(E1.p1,ok)
    E2.Ns=filt(E1.Ns,ok)
    E2.Nm=filt(E1.Nm,ok)
    E2.eNr=filt(E1.eNr,ok)
    E2.Nf5=filt(E1.Nf5,ok)
    E2.Nf3=filt(E1.Nf3,ok)
    E2.p5a=filt(E1.p5a,ok)
    E2.p5b=filt(E1.p5b,ok)
    E2.p3a=filt(E1.p3a,ok)
    E2.p3b=filt(E1.p3b,ok)
    E2.cn=filt(E1.cn,ok)
    E2.posU=filt(E1.posU,ok)
    E2.lenU=filt(E1.lenU,ok)
    E2.src5=filt(E1.src5,ok)
    E2.src3=filt(E1.src3,ok)
    E2.merge5=filt(E1.merge5,ok)
    E2.merge3=filt(E1.merge3,ok)
    E2.N=len(E2.a)
    return E2

def goodEventCrx(params,E1):

    bothsides=params["bothsides"]
    cnslosh=params["cnslosh"]
    EventType=params["EventType"]

    crx= (EventType=='cross')

    if (crx<1):
         return E1

    E2=E1

    ok=[]
    for e in range(E1.N):
        badsides= (bothsides & ((E1.Nf5[e]<1)|(E1.Nf3[e]<1)))
        badcn= (E1.cn[e]>(2+cnslosh))
        #if ((E1.Nf5[e]>0)&(E1.Nf3[e]>0)):
        if not (badsides|badcn):
            ok.append(int(e))

#    xxx=[0,1,2,3,4,5]
#    ok=[0,2,5]
#    xx=filt(xxx,ok)

    E2.p=filt(E1.p,ok)
    E2.l=filt(E1.l,ok)
    E2.a=filt(E1.a,ok)
    E2.q=filt(E1.q,ok)
    E2.p2=filt(E1.p2,ok)
    E2.l2=filt(E1.l2,ok)
    E2.a2=filt(E1.a2,ok)
    E2.Nr=filt(E1.Nr,ok)
    E2.p0=filt(E1.p0,ok)
    E2.p1=filt(E1.p1,ok)
    E2.Ns=filt(E1.Ns,ok)
    E2.Nm=filt(E1.Nm,ok)
    E2.eNr=filt(E1.eNr,ok)
    E2.Nf5=filt(E1.Nf5,ok)
    E2.Nf3=filt(E1.Nf3,ok)
    E2.p5a=filt(E1.p5a,ok)
    E2.p5b=filt(E1.p5b,ok)
    E2.p3a=filt(E1.p3a,ok)
    E2.p3b=filt(E1.p3b,ok)
    E2.p5aT=filt(E1.p5aT,ok)
    E2.p5bT=filt(E1.p5bT,ok)
    E2.p3aT=filt(E1.p3aT,ok)
    E2.p3bT=filt(E1.p3bT,ok)
    E2.cn=filt(E1.cn,ok)
    E2.posU=filt(E1.posU,ok)
    E2.lenU=filt(E1.lenU,ok)
    E2.src5=filt(E1.src5,ok)
    E2.src3=filt(E1.src3,ok)
    E2.merge5=filt(E1.merge5,ok)
    E2.merge3=filt(E1.merge3,ok)
    E2.N=len(E2.a)
    return E2


def filehook(params, directory, files):
    """os.walk hook for *.library.span file for checking to make Spanner detect links"""

    global E
    global first

    EventType=params["EventType"]

    pattern = '*.'+EventType+'.span.txt'


    for file in files:
        pathname = os.path.join(directory, file)
        #print file

        if fnmatch.fnmatch(file, pattern):
            if os.path.isfile(pathname):

                if ((EventType=='del')|(EventType=='dup')):
                    E1=EventSpan.Event1Span(pathname)
                    N1=E1.N
                    E2=goodEvent1(params,E1)
                    print '# keep %d of %d %s ' % (E2.N, N1, E1.typeName)
                    if (first>0):
                        E=E2;
                    else:
                        E1=mergeEvt(E,E2)
                        E=E1

                elif (EventType=='inv'):
                    E1=EventSpan.EventInv1Span(pathname)
                    N1=E1.N
                    E2=goodEventInv(params,E1)
                    print '# keep %d of %d %s ' % (E2.N, N1, E1.typeName)
                    if (first>0):
                        E=E2;
                    else:
                        E1=mergeInvEvt(E,E2)
                        E=E1

                elif ((EventType=='ALU')|(EventType=='L1')|(EventType=='SVA')):
                    E1=EventSpan.EventMob1Span(pathname)
                    N1=E1.N
                    E2=goodEventMob(params,E1)
                    print '# keep %d of %d %s ' % (E2.N, N1, E1.typeName)
                    if (first>0):
                        E=E2;
                    else:
                        E1=mergeMobEvt(E,E2)
                        E=E1

                elif (EventType=='cross'):
                    E1=EventSpan.EventCrx1Span(pathname)
                    N1=E1.N
                    E2=goodEventCrx(params,E1)
                    print '# keep %d of %d %s ' % (E2.N, N1, E1.typeName)
                    if (first>0):
                        E=E2;
                    else:
                        E1=mergeCrxEvt(E,E2)
                        E=E1


                if E:
                    first=0


def main():

    parser = optparse.OptionParser()
    parser.add_option("-i", "--inarea",
                  help="input event area")
    parser.add_option("-o", "--outbed",
                  help="output bed file")
    parser.add_option("-n", "--Nfragmin",default="2",
                  help="minimum number of supporting fragments")
    parser.add_option("-e", "--EventType",default="del",
                  help="event type: del, dup, ins, inv, alu, l1, sva")
    parser.add_option("-c", "--cnslosh",default="0",
                  help="cnslosh cnv coverage window: 0:none, 1:CN=2 ok, 3:wrong poloary ok ...")
    parser.add_option("-b", "--bothsides",default="1",
                  help=">0 forces each event to have support from bothsides")


#   parser.add_option("-m", "--mask",default="",
#                  help="masking bed file")


    (options, args) = parser.parse_args()

    if (options.inarea==None)|(options.outbed==None):
        parser.error("need at least -i and -o arguments")

    inarea = options.inarea # sys.argv[1]

    print "# Spanner event input area=",inarea

    outbed = options.outbed   #sys.argv[2]

    print "# output bed file=",outbed

    Nfragmin = int(float(options.Nfragmin))    #int(float(sys.argv[3]))
    bothsides= int(float(options.bothsides))
    cnslosh= int(float(options.cnslosh))

    EType = options.EventType

    if (EType=='alu'):
        EType='ALU'
    elif (EType=='l1'):
        EType='L1'
    elif (EType=='sva'):
        EType='SVA'

    print "EType %s\t Nfragmin>=%d\t CNslosh>=%d\t Bothsides>=%d" % ( EType,Nfragmin,cnslosh,bothsides)

#    MaskFile = options.mask
#
#    if (len(MaskFile)>0):
#        print "mask file =%s" % MaskFile


    origdir = os.getcwd()

    params={}
    params["Nfragmin"]=Nfragmin
    params["bothsides"]=bothsides
    params["cnslosh"]=cnslosh
    params["EventType"]=EType
    params["outbed"]=outbed
#    params["MaskFile"]=MaskFile
    


    print "# pwd=",origdir

    os.chdir(inarea)

    os.path.walk(os.curdir, filehook, params)

    os.chdir(origdir)

    if (EType=='cross'):
        printCrossBed(params,E)
    else:
        printBed(params,E)


if __name__ == "__main__":
    main()




