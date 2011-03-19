/*
 *  SpanDet.h
 *  SpanDet
 *
 *  Created by Chip Stewart on 8/17/08.
 *  Copyright 2008 Boston College. All rights reserved.
 *
 */
#ifndef SPANDET_H
#define SPANDET_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <time.h>
//#include <algorithm>
#include <map>
#include <list>
#include <iterator>
#include <math.h>
// #include <boost/regex.hpp>
#include "Type-Hash.h"
#include "RunControlParameterFile.h"
#include "Histo.h"
#include "Function-Generic.h"
#include "PairedData.h"
#include "cluster.h"
#include "SpannerVersion.h"
#include "BedFile.h"

using namespace std;


//-----------------------------------------------------------------------------
// Nominal coverage class
//-----------------------------------------------------------------------------
typedef std::map<int, HistObj, less<int> >  C_HistObjs;
typedef std::vector<double>  C_vectorDouble;

//------------------------------------------------------------------------------
// sorts keys of a hash in order of associated value
// Why the hell doesn't this work from Function-Generic.cpp ???
//------------------------------------------------------------------------------
template< typename keyType, typename valueType >
vector<keyType> sortKeysByValue(map<keyType, valueType, std::less<keyType> > hash, bool descend) {

  // instantiate inverse hash as a multimap
  multimap<valueType, keyType, std::less<valueType> > inverseHash;

  // load elements of hash into inverseHash
  for (typename map<keyType, valueType, std::less<keyType> >::const_iterator iter = hash.begin(); iter != hash.end(); iter++) {
    keyType key = iter->first;
    valueType value = iter->second;
    inverseHash.insert(typename multimap<valueType, keyType, std::less<valueType> >::value_type(value, key));
  }

  // compose vector of original keys sorted by original values
  vector<keyType> sortedKeys;
  for(typename multimap<valueType, keyType, std::less<valueType> >::const_iterator iter = inverseHash.begin(); 
      iter != inverseHash.end(); iter++) {
    keyType key = iter->second;
    sortedKeys.push_back(key);
  }

  // reverse if descending order was required
  if (descend) {
    reverse(sortedKeys.begin(), sortedKeys.end());
  }

  // return
  return sortedKeys;
}


class C_NominalCov {
  friend ostream &operator<<(ostream &, C_NominalCov &);
  public:
    C_NominalCov();    
    C_NominalCov(C_set &, vector<int> &);
    C_HistObjs hist;
    double getplow(int,double);
    int N;
    vector<int> nu;
  private:
    double interp(int,double,int,double,int);
};
   
//-----------------------------------------------------------------------------
// single event read coverage class
//-----------------------------------------------------------------------------
class C_SVcoverage1 {
  friend ostream &operator<<(ostream &, const C_SVcoverage1 &);
  //friend class C_SV1;
  public:
    C_SVcoverage1();
    C_SVcoverage1(C_contig &, int, int, C_NominalCov &);
    C_SVcoverage1&operator=(const C_SVcoverage1 &rhs);
    int operator==(const C_SVcoverage1 &rhs) const;
    int operator<(const C_SVcoverage1 &rhs) const;
    int p0;               // start base position of region
    int p1;               // end base position
    int Nsite;            // number of available start sites in region
    int Nrepeat;          // number of repeat start sites in region
    int N;                // number of reads 
    float eN;             // expected number of reads
    double p;             // p value for null CNV (poisson)
};

typedef std::map<unsigned int,unsigned int, less<unsigned int> >  C_RGmap;

class C_SVspanfrags1 {
  friend ostream &operator<<(ostream &, const C_SVspanfrags1 &);
  public:
    C_SVspanfrags1();
    C_SVspanfrags1(unsigned int, unsigned int, unsigned int, unsigned int, float);
    C_SVspanfrags1&operator=(const C_SVspanfrags1 &rhs);
    int operator==(const C_SVspanfrags1 &rhs) const;
    unsigned int N;       // count of supporting fragments 
    unsigned int NN;      // count of non-supporting spanning fragments 
    unsigned int N5;      // count of  5' supporting fragments 
    unsigned int N3;      // count of  4' supporting fragments 
    unsigned int NR;      // count of  reads mapped to SV region 
    float        ER;      // estimated count of  read mapped to SV region     
};


class C_SVCF {
  friend ostream &operator<<(ostream &, C_SVCF &);
  public:
    C_SVCF();    
    C_SVCF(RunControlParameters &, vector <string> & Info1,vector <string> & Format1, vector <string> & Samples1);
    string Version;
    string Date;
    string Source;
    string Reference;
    string EventType;
    int NEvent;
    vector<string> Info;
    vector<string> Format;
    vector<string> Samples;
    
};

typedef std::map<string, C_SVspanfrags1, std::less<string> >   C_SAMmap;

//-----------------------------------------------------------------------------
// single event class
//-----------------------------------------------------------------------------
class C_SV1 {
  friend ostream &operator<<(ostream &, const C_SV1 &);
  //friend class C_SVcoverage1;
  //friend class C_cluster2d_element1;
  public:
    C_SV1();
    C_SV1&operator=(const C_SV1 &rhs);
    int operator==(const C_SV1 &rhs) const;
    int operator<(const C_SV1 &rhs) const;
    string type;
    unsigned int pos;
    unsigned short anchor;
    unsigned int length;
    unsigned char q;            // -10log (p nominal coverage)
    unsigned int p5[2];
    unsigned int p3[2];    
    unsigned short copynumber;  
    C_cluster2d_element1 cls;
    vector<C_localpair> pair;
    C_SVcoverage1 cov;          // coverage inside event
    C_SVcoverage1 cov5;         // coverage at 5' end of event
    C_SVcoverage1 cov3;         // coverage at 5' end of event
    float a5;                   // alignability estimate at 5' end of event
    float a3;                   // alignability estimate at 3' end of event
    unsigned int posU;          // position uncertainty
    unsigned int lenU;          // length uncertainty
    unsigned int qOutlier;      // probability of outlier in cluster
    unsigned int qAberrantLM;   // probability of aberrant frag length
    C_RGmap ReadGroupMap;       // count of fragments for each read in event
    C_SAMmap SampleMap;         // count of fragments for each sample in library
    bool merge;                 // merged clusters flag
    int id;

};    


//-----------------------------------------------------------------------------
// base class for SV event lists
//-----------------------------------------------------------------------------
class C_SV {
  friend ostream &operator<<(ostream &, const C_SV &);
  public:
    C_SV();
    C_SV(C_contig  &,  RunControlParameters &);
    void print(string &);
    list<C_SV1>  evt;
    string typeName;
    string contigName;
    string setName; 
    void finalize(C_contig  &, C_libraries &, RunControlParameters &);
    void genotype(C_contig  &, C_libraries &, RunControlParameters &);
    C_SVCF SVCF; 
    vector<string> samples;
    C_BedChr Mask;      
        
}; 

//-----------------------------------------------------------------------------
// single retro event class
//-----------------------------------------------------------------------------
class C_SVR1 {
  friend ostream &operator<<(ostream &, const C_SVR1 &);
  public:
    C_SVR1();
    C_SVR1&operator=(const C_SVR1 &rhs);
    int operator==(const C_SVR1 &rhs) const;
    int operator<(const C_SVR1 &rhs) const;
    string type;
    unsigned int pos;
    unsigned short anchor;
    unsigned int length;
    unsigned char q;  
    unsigned int p5[2];
    unsigned int p3[2];    
    unsigned short copynumber;  
    C_SVcoverage1 cov;  
    // retro event info
    C_cluster2d_element1 cls5;    
    C_cluster2d_element1 cls3;
    vector<C_umpair> retro5;
    vector<C_umpair> retro3;
    unsigned int NfragCluster[2];
    unsigned int NfragCov;
    unsigned int NfragCovOut[2];
    unsigned int NfragCovExp;
    unsigned int pmedian;
    unsigned int Nconstrain[2];  
    // 8/2009
    float a5;                   // alignability estimate at 5' end of event
    float a3;                   // alignability estimate at 3' end of event
    unsigned int posU;          // position uncertainty
    unsigned int lenU;          // length uncertainty
    C_RGmap ReadGroupMap5;      // count of fragments for each read in 5' frags
    C_RGmap ReadGroupMap3;      // count of fragments for each read in 3' frags
    C_SAMmap SampleMap;         // count of fragments for each sample in library
    bool merge5;                // merged clusters flag
    bool merge3;                // merged clusters flag
    string subtype;             // subtype most likely
    float Qsubtype;             // phred Q subtype
    float QsubtypeI;            // phred Q subtype level I (first character after type)
    string subtype5;            // 5' end subtype most likely
    float Qsubtype5;            // 5' end phred Q subtype
    string subtype3;            // 3' end subtype most likely
    float Qsubtype3;            // 3' end phred Q subtype    
    int gap;                    // length of element inserted
    float Esense;               // sense (F/R) of element
    int id;                     // id for this event
};    

//-----------------------------------------------------------------------------
// Retro class for SV event lists
//-----------------------------------------------------------------------------
class C_SVR {
  friend ostream &operator<<(ostream &, const C_SVR &);
  public:
    C_SVR();
    C_SVR(C_contig  &,  RunControlParameters &);
    void print(string &);
    void finalize(C_contig  &, C_libraries &, RunControlParameters &);
    void genotype(C_contig  &, C_libraries &, RunControlParameters &);
    C_SVCF SVCF; 
    list<C_SVR1>  evt;
    string typeName;
    string contigName;
    string setName;    
    vector<string> samples;
    C_BedChr Mask;       
}; 

//-----------------------------------------------------------------------------
// single inversion event class
//-----------------------------------------------------------------------------
class C_SVV1 {
  friend ostream &operator<<(ostream &, const C_SVV1 &);
  public:
    C_SVV1();
    C_SVV1&operator=(const C_SVV1 &rhs);
    int operator==(const C_SVV1 &rhs) const;
    int operator<(const C_SVV1 &rhs) const;
    string type;
    unsigned int pos;
    unsigned short anchor;
    unsigned int length;
    unsigned char q;  
    unsigned int p5[2];
    unsigned int p3[2];    
    unsigned short copynumber;  
    C_SVcoverage1 cov;  
    // retro event info
    C_cluster2d_element1 cls5;    
    C_cluster2d_element1 cls3;
    vector<C_localpair> pair5;
    vector<C_localpair> pair3;
    unsigned int NfragCluster[2];
    unsigned int NfragCov;
    unsigned int NfragCovOut[2];
    unsigned int NfragCovExp;
    // August 2009 vars ...
    C_SVcoverage1 cov5;         // coverage at 5' end of event
    C_SVcoverage1 cov3;         // coverage at 5' end of event
    float a5[2];                // alignability estimate at 5' end of event
    float a3[2];                // alignability estimate at 3' end of event
    unsigned int posU;          // position uncertainty
    unsigned int lenU;          // length uncertainty
    unsigned int qOutlier;      // probability of outlier in cluster
    unsigned int qAberrantLM;   // probability of aberrant frag length
    C_RGmap ReadGroupMap5;      // count of fragments for each read in 5' frags
    C_RGmap ReadGroupMap3;      // count of fragments for each read in 3' frags
    C_SAMmap SampleMap;         // count of fragments for each sample in library
    bool merge5;                 // merged clusters flag
    bool merge3;                 // merged clusters flag
    int id;
};    

//-----------------------------------------------------------------------------
//  class for SV inversion event lists
//-----------------------------------------------------------------------------
class C_SVV {
  friend ostream &operator<<(ostream &, const C_SVV &);
  public:
    C_SVV();
    C_SVV(C_contig  &,  RunControlParameters &);
    void print(string &);
    list<C_SVV1>  evt;
    string typeName;
    string contigName;
    string setName;    
    void finalize(C_contig  &, C_libraries &, RunControlParameters &);
    void genotype(C_contig  &, C_libraries &, RunControlParameters &);
    C_SVCF SVCF; 
    vector<string> samples;
    C_BedChr Mask;      
}; 

//-----------------------------------------------------------------------------
// inter-chromosomal link event class
//-----------------------------------------------------------------------------
class C_SVX1 {
  friend ostream &operator<<(ostream &, const C_SVX1 &);
  public:
    C_SVX1();
    C_SVX1&operator=(const C_SVX1 &rhs);
    int operator==(const C_SVX1 &rhs) const;
    int operator<(const C_SVX1 &rhs) const;
    string type;
    unsigned int pos;
    unsigned short anchor;
    unsigned int length;
    unsigned char q;  
    unsigned int p5[2];
    unsigned int p3[2];    
    unsigned short anchor2;    
    unsigned int pos2;
    unsigned int length2;
    unsigned int p5a2[2];
    unsigned int p3a2[2];    
    unsigned short copynumber;  
    C_SVcoverage1 cov;  
    // retro event info
    C_cluster2d_element1 cls5;    
    C_cluster2d_element1 cls3;
    vector<C_crosspair> cross5;
    vector<C_crosspair> cross3;
    unsigned int NfragCluster[2];
    unsigned int NfragCov;
    unsigned int NfragCovOut[2];
    unsigned int NfragCovExp;
    // August 2009 vars ...
    C_SVcoverage1 cov5;         // coverage at 5' end of event
    C_SVcoverage1 cov3;         // coverage at 5' end of event
    float a5[2];                // alignability estimate at 5' end of event
    float a3[2];                // alignability estimate at 3' end of event
    unsigned int posU;          // position uncertainty
    unsigned int lenU;          // length uncertainty
    unsigned int qOutlier;      // probability of outlier in cluster
    C_RGmap ReadGroupMap5;      // count of fragments for each read in 5' frags
    C_RGmap ReadGroupMap3;      // count of fragments for each read in 3' frags
    C_SAMmap SampleMap;         // count of fragments for each sample in library
    bool merge5;                // merged clusters flag
    bool merge3;                // merged clusters flag
    int id;                     // id for this event

};    

//-----------------------------------------------------------------------------
//  class for SV cross linked event lists
//-----------------------------------------------------------------------------
class C_SVX {
  friend ostream &operator<<(ostream &, const C_SVX &);
  public:
    C_SVX();
    C_SVX(C_contig  &,  RunControlParameters &);
    void print(string &);
    C_SVCF SVCF; 
    void finalize(C_contig  &, C_libraries &, RunControlParameters &);
    void genotype(C_contig  &, C_libraries &, RunControlParameters &);
    list<C_SVX1>  evt;
    string typeName;
    string contigName;
    string setName; 
    vector<string> samples;
    C_BedChr Mask;   
}; 



//-----------------------------------------------------------------------------
// clustering
//-----------------------------------------------------------------------------
class C_SpannerCluster {
  friend ostream &operator<<(ostream &, const C_SpannerCluster &);
  public:
    C_SpannerCluster();
    C_SpannerCluster(C_contig &, C_libraries &, RunControlParameters &);    
    
    void writeall();
    
    // library info
    C_libraries libraries;
    
    string typeName;
    string contigName;
    string setName;    
    C_NNcluster2d longpairc;
    C_NNcluster2d shortpairc;
    C_NNcluster2d invert5c;
    C_NNcluster2d invert3c;

    // novel insertion clusters
    C_NNcluster1d dangle5c;
    C_NNcluster1d dangle3c;
    
    // select abberant localpair read vectors
    vector<C_localpair>  invert5;
    vector<C_localpair>  invert3;
    vector<C_localpair>  longpair;
    vector<C_localpair>  shortpair;

    // cross linked fragment clusters
    C_NNcluster2d cross5c;
    C_NNcluster2d cross3c;

    // selected cross pair read vectors
    vector<C_crosspair>  cross5;
    vector<C_crosspair>  cross3;


 private:
    int setClusterWindow(C_contig &,  RunControlParameters &);
    int makeDangleX(C_contig &, char, vector<double> & );
    void selectPairs(C_contig &, int );    
    void selectCross(C_contig &, int );    
    int makepairP(vector<C_localpair> &, char , vector<vector<double> > &, 
        vector<vector<double> > &, vector<int> &);
     int makepairX(vector<C_crosspair> &, char , vector<vector<double> > &, 
        vector<vector<double> > &, vector<int> &);
    void write(string &, C_NNcluster2d &, vector<C_localpair> &);
    void write(string &, C_NNcluster2d &, vector<C_crosspair> &);
    int L;
    double PairDensity;
    double EndDensity;
    RunControlParameters pars;    
}; 

//-----------------------------------------------------------------------------
// clustering
//-----------------------------------------------------------------------------
class C_SpannerRetroCluster {
  friend ostream &operator<<(ostream &, const C_SpannerRetroCluster &);
  public:
    C_SpannerRetroCluster();
    C_SpannerRetroCluster(C_contig &,  C_libraries & libs1, RunControlParameters &, int, string &);    
    void writeall();
    int Mask(C_BedChr &, int);
    string typeName;
    string contigName;
    string setName;       
    // retro element type bit number 
    int etype;
    // mobile element insertions clusters
    C_NNcluster2d e5c;
    C_NNcluster2d e3c;

    // select mobile element umpairs
    vector<C_umpair>  e5;
    vector<C_umpair>  e3;
    
    // library info
    C_libraries libraries;
    
 private:
    int setClusterWindow(C_contig &,  RunControlParameters &);
    void selectRetro(C_contig &, int);  // mobile elements
    void write(string &, C_NNcluster2d &, vector<C_umpair> &);
    int L;
    double PairDensity;
    double RetroDensity;
    // obsolete Aug 2009
    int makeRetroX(vector<C_umpair> &, char, vector<double> & );
    // with lib info aug 2009
    int makepairP(vector<C_umpair> &, char , vector<vector<double> > &, 
        vector<vector<double> > &, vector<int> &);

    RunControlParameters pars;    
}; 

//-----------------------------------------------------------------------------
// get Retro elements from anchor info
//-----------------------------------------------------------------------------
class C_retroElements {
  public:
  C_retroElements(C_anchorinfo &);
  vector <int> e;
  vector <string> name;
  int N;
};

//-----------------------------------------------------------------------------
// Spanner SV detection main class
//-----------------------------------------------------------------------------
class C_SpannerSV {
  friend ostream &operator<<(ostream &, const C_SpannerSV &);
  public:
    C_SpannerSV();  
    C_SpannerSV(C_pairedfiles &,  RunControlParameters &);
    C_anchorinfo anchor;
    C_libraries libraries;
    C_NominalCov nomcov;
    C_SV del;    
    C_SV dup;
    C_SVV inv;
    C_SVR ret;
    C_SVX crx;
  private:
    int findDel(C_contig  &, C_SpannerCluster &, RunControlParameters &);
    C_SV1 merge(C_SV1 &, C_SV1 &, C_contig &, RunControlParameters &, int);
    int findDup(C_contig  &, C_SpannerCluster &, RunControlParameters &);
    int findInv(C_contig  &, C_SpannerCluster &, RunControlParameters &);
    C_SVV findInvDir(C_contig  &, C_SpannerCluster &, RunControlParameters &, int);
    C_SVV1 merge(C_SVV1 &, C_SVV1 &, C_contig & );
    int findCross(C_contig  &, C_SpannerCluster &, RunControlParameters &);
    C_SVX findCrossDir(C_contig  &, C_SpannerCluster &, RunControlParameters &, int);
    C_SVX1 merge(C_SVX1 &, C_SVX1 &, C_contig & );
    int findRet(C_contig  &, C_SpannerRetroCluster &, RunControlParameters & );
    int findRet0(C_contig  &, C_SpannerRetroCluster &,  RunControlParameters & ); // old way
    C_SVR1  makeRetEvent(C_contig  &,  C_cluster2d_element1 &,   vector<C_umpair> &,double);
    C_SVR1 merge(C_SVR1 &, C_SVR1 &, C_contig &);
    C_NNcluster2d cullRet(C_contig  &, C_NNcluster2d & ,  vector<C_umpair> & , RunControlParameters & );
    
    char p2q(double);
    void printSummary( string & ,C_contig &);
    
}; 

#endif

