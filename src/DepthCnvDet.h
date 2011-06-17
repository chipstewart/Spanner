/*
 *  DepthCnvDet.h
 *  SpanDet
 *
 *  Created by Chip Stewart on 10/17/08.
 *  Copyright 2008 Boston College. All rights reserved.
 *
 */
/*
 *  SpanDet.h
 *  SpanDet
 *
 *  Created by Chip Stewart on 8/17/08.
 *  Copyright 2008 Boston College. All rights reserved.
 *
 */
#ifndef DEPTHCNVDET_H
#define DEPTHCNVDET_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <iterator>
#include <math.h>
//#include <boost/regex.hpp>
#include "Type-Hash.h"
#include "RunControlParameterFile.h"
#include "Histo.h"
#include "Function-Generic.h"
#include "PairedData.h"
#include "cluster.h"
#include "SpanDet.h"
#include "steps.h"


using namespace std;
//-----------------------------------------------------------------------------
// Unique data read coverage class
//-----------------------------------------------------------------------------
class C_UniqueCoverage {
  friend ostream &operator<<(ostream &, const C_UniqueCoverage &);
  public:
    C_UniqueCoverage();
    C_UniqueCoverage(C_contig &,  RunControlParameters &);
    // reads per binsize
    C_depth r;
    // RDreference model
    C_depth a;
    // expected reads per bin (null CNV hypothesis)   
    C_depth e;
    // GC content averaged per bin 
    C_depth gc;
    // GC content correction factor
    C_depth gcf;
    // read count bin size (bp)
    int binsize;
    //int abintotal;
    int slope;
    void write(C_contig &, string &);
    int rebin(RunControlParameters &);
    int gc_correctionFactor();
};

//-----------------------------------------------------------------------------
// single CNV event class
//-----------------------------------------------------------------------------
class C_CNV1 {
  friend ostream &operator<<(ostream &, const C_CNV1 &);
  public:
    C_CNV1() {};
    C_CNV1&operator=(const C_CNV1 &rhs);
    //int operator==(const C_CNV1 &rhs) const;
    //int operator<(const C_CNV1 &rhs) const;
    unsigned int pos;
    unsigned short anchor;
    unsigned int length;
    unsigned short copynumber;  
    unsigned char q;  
    double nr;                // total number of reads here
    double enr;               // expected number of reads here
    int nbi;                  // number of informative bins in this evt
    int nbin;                 // number of bins in this evt
    double anr;               // mean number of reads in bins
    double snr;               // stdev number of reads in bins
    double aenr;              // mean expected number of reads in bins
    double senr;              // stdev expected number of reads in bins
};    

//-----------------------------------------------------------------------------
// base class for CNV event lists
//-----------------------------------------------------------------------------
class C_CNV {
  friend ostream &operator<<(ostream &, C_CNV &);
  public:
    C_CNV() {};
    C_CNV(C_contig  &,  RunControlParameters &);
    list<C_CNV1>  evt;
    string typeName;
    string contigName;
    string setName;    
    void print(string &);
}; 

//-----------------------------------------------------------------------------
// Depth CNV detector main class 
//-----------------------------------------------------------------------------
class C_DepthCNV {
  friend ostream &operator<<(ostream &, const C_DepthCNV &);
  public:
    C_DepthCNV() {};  
    C_DepthCNV(C_pairedfiles &,  RunControlParameters &);
    C_anchorinfo anchor;
    C_UniqueCoverage u;
    C_steps steps;
    vector<int> bb;
    int Ndel;
    C_CNV del;
    int Ndup;    
    C_CNV dup;
    vector<int> mx;
  //private:
    int findDel(C_contig  &, RunControlParameters &);
    int findDup(C_contig  &, RunControlParameters &);
    vector<C_CNV1> getinfo(C_contig  &, RunControlParameters &, int);
    int MINBINS;
    int BMAX;
    int KBMAX;
    int ALLOW;
    int GAP;
    //char p2q(double);
    void printSummary( string & ,C_contig &);
}; 

#endif


