#! /usr/bin/python

# ----------------------------------------
# Copyright 2009  Chip Stewart, Boston College
# All rights reserved.
# ----------------------------------------

#import sys
import os
import os.path
import fnmatch
import re
import optparse
import LibrarySpan

first=1

__author__="chip stewart"
__date__ ="$Oct 29, 2009 5:02:13 PM$"

def makeLinks(params,L,subdirectory,pathname):
    """make soft links for each reference anchor"""

    global first
    linkarea=params["linkarea"].rstrip("/")
    makescript=params["makescript"].rstrip()
    verbose=params["verbose"]

    directory = os.getcwd()

    if makescript:
        if (first):
            f= open(makescript, 'w')
            print >>f, "mkdir -p ", linkarea

        else:
            f= open(makescript, 'a')

    elif not os.path.exists (linkarea):
        os.makedirs (linkarea)

    for ia in range(L.Na):
        if (L.a[ia]["use"]):
            A=L.a[ia]["name"].rstrip('\x00')
            linkareaA=linkarea+"/"+A

            filename=directory+pathname.strip('.')

            if os.path.exists (filename):
                if makescript:
                    if (first):
                        print >>f, "mkdir -p ", linkareaA
                
                    print >>f, "cd ", linkareaA

                    if ('special' in filename):                      
                        nospecial = filename.replace('special.', '')
                        (dirName, nospecialfilename) = os.path.split(nospecial)
                        print >>f, "ln -sf", filename, " ",nospecialfilename
                    else:
                        print >>f, "ln -s", filename

                else:
                    if not os.path.exists (linkareaA):
                        os.makedirs (linkareaA)

                    os.chdir(linkArea)

                    if ('special' in filename):
                        if (os.path.islink(filename)):
                            os.remove(filename)

                        nospecial = filename.replace('special.', '')
                        os.symlink(filename,nospecial)
                    else:
                        if not (os.path.islink(filename)):
                          os.symlink(filename)

            else:
                if(verbose):
                    print "# missing ",filename

            basename=re.sub("\.library\.span","" , filename)

            filename= basename+".anchors.txt"

            if os.path.exists (filename):
                if makescript:
                    if ('special' in filename):                     
                        nospecial = filename.replace('special.', '')
                        (dirName, nospecialfilename) = os.path.split(nospecial)
                        print >>f, "ln -sf", filename, " ",nospecialfilename
                    else:
                        print >>f, "ln -s", filename

                else:
                    if ('special' in filename):
                        if (os.path.islink(filename)):
                            os.remove(filename)
                        nospecial = filename.replace('special.', '')
                        os.symlink(filename,nospecial)
                    else:
                        if not (os.path.islink(filename)):
                            os.symlink(filename)
            else:
                if (verbose):
                    print "# missing ",filename

            #skip special pair,cross, multi span links
            if ('special' in filename):
                continue

            filename= basename+"."+A+".pair.span"

            if os.path.exists (filename):
                if makescript:
                    print >>f, "ln -s", filename
                else:
                    os.symlink(filename)
            else:
                if (verbose):
                    print "# missing ",filename


            filename= basename+"."+A+".cross.span"

            if os.path.exists (filename):
                if makescript:
                    print >>f, "ln -s", filename
                else:
                    os.symlink(filename)
            else:
                if (verbose):
                    print "# missing ",filename


            filename= basename+"."+A+".multi.span"

            if os.path.exists (filename):
                if makescript:
                    print >>f, "ln -s", filename
                else:
                    os.symlink(filename)
            else:
                if (verbose):
                    print "# missing ",filename



            filename= basename+"."+A+".dangle.span"

            if os.path.exists (filename):
                if makescript:
                    print >>f, "ln -s", filename
                else:
                    os.symlink(filename)
            else:
                if (verbose):
                    print "# missing ",filename



    if makescript:
        f.close()

    os.chdir(directory)

    first = 0

def goodLib(params,L):

   PairMin=params["PairMin"]
   ok =  (L.NPair[0]>PairMin)

   LFpeak=params["LFpeak"]
   peak=L.LFhist[0]["std"]/L.LFhist[0]["median"]
   ok =  ok&(peak<LFpeak)&(L.LFhist[0]["median"]>0)

   return ok



def filehook(params, directory, files):
    """os.walk hook for *.library.span file for checking to make Spanner detect links"""


    pattern = '*.library.span'

    for file in files:
        pathname = os.path.join(directory, file)
        #print file
        # skip special libraries (should be matched to regular)
        #if ('special.' in file):
        #    continue

        if fnmatch.fnmatch(file, pattern):
            if os.path.isfile(pathname):
                L=LibrarySpan.LibrarySpan(pathname)
                peak=L.LFhist[0]["std"]/L.LFhist[0]["median"]

                if goodLib(params,L):
                  print '# keep ', pathname, L.setName, L.SampleName[0], L.NPair[0], L.LM[0],peak
                  makeLinks(params,L,directory,pathname)

                else:
                  print '# skip ', pathname, L.setName, L.SampleName[0], L.NPair[0], L.LM[0],peak


def main():

    parser = optparse.OptionParser()
    parser.add_option("-i", "--inarea",
                  help="input build area")
    parser.add_option("-o", "--outarea",
                  help="output link area")
    parser.add_option("-n", "--Npairmin",default="-1",
                  help="minimum number of read pairs in file")
    parser.add_option("-v", "--verbose",dest='verbose', default=False,
                  help="detailed standard output ",action='store_true')
    parser.add_option("-p", "--peak",default="0.99",
                  help="max allow std/median fragment distribution")
    parser.add_option("-m", "--makescript",default="",
                  help="makelink script")



    (options, args) = parser.parse_args()

    if (options.inarea==None)|(options.outarea==None):
        parser.error("need at least -i and -o arguments")


    buildarea = options.inarea # sys.argv[1]

    print "# buildarea=",buildarea

    linkarea = options.outarea   #sys.argv[2]

    print "# linkarea=",linkarea

    PairMin = int(float(options.Npairmin))    #int(float(sys.argv[3]))

    LFpeak = float(options.peak)    #ifloat(sys.argv[4]))

    print "# NPairMin=",PairMin,"LFpeak=",LFpeak

    verbose = options.verbose    #ifloat(sys.argv[4]))

    if (verbose):
        print "# Verbose stdout"


    origdir = os.getcwd()

    params={}
    params["PairMin"]=PairMin
    params["LFpeak"]=LFpeak
    params["linkarea"]=linkarea
    params["makescript"]=origdir+"/"+options.makescript
    params["verbose"]=verbose



    print "# pwd=",origdir

    os.chdir(buildarea)

    os.path.walk(os.curdir, filehook, params)

    os.chdir(origdir)


if __name__ == "__main__":
    main()




