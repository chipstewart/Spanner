/*
 *  steps.h
 *  SpanDet
 *
 *  Created by Chip Stewart on 10/24/08.
 *  Copyright 2008 Boston College. All rights reserved.
 *
 */
#ifndef STEPS_H
#define STEPS_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <list>
#include <algorithm>
#include <math.h>
#include <dirent.h> 
#include <stdio.h> 
// local includes
#include "Histo.h"
#include "Function-Generic.h"

// float HUGE = 1e30;


//------------------------------------------------------------------------------
// class for vector x outlier removal
//------------------------------------------------------------------------------
class C_trim {
  public:
    C_trim(vector<float> & x0,float STD=0.4f,int W=3, int CUT=2);
    float std;
    float w;
    float cut;  
    vector<float> x;  //		trimmed list  (size equal or smaller than original x list)
    vector<int>   m;  //    indices of x  (1:N with holes for cut x)
    vector<float> xx;	//  	list of removed original x
    vector<int>   mx; //   	indices of removed original x
};

class C_trimstat {
  public:
    C_trimstat(vector<float> & x0,float trim=0.05f);
    int N0;
    int N;
    float std;
    float mean;  
    float median;  
};

//------------------------------------------------------------------------------
// function to gaussian prob of sampling value > x given mu and sigma
//------------------------------------------------------------------------------
double pnorm(const double x0, const double mu=0, const double sigma=1.0f);

//------------------------------------------------------------------------------
// function to calcuate trimmed mean 
//------------------------------------------------------------------------------
float meantrimmed(vector<float> & X,float trim=0.05f);

//------------------------------------------------------------------------------
// function to convolute a with tap - smoothing filter
//------------------------------------------------------------------------------
vector <float> convolve(vector <float> & a, vector <float> & tap); 

//------------------------------------------------------------------------------
// function to find a set of break candidates around peaks 
//------------------------------------------------------------------------------
list <int> peakBounds(vector <float> & y, float t, int flip);

//------------------------------------------------------------------------------
// template functions to do quick & dirty printout of vectors or lists
//------------------------------------------------------------------------------
template< typename vType>
void vdump(vector <vType> & , const string & );
template< typename wType>
void ldump(list <wType> & , const string & );

//------------------------------------------------------------------------------
// indexlist for extracting order of sort
//------------------------------------------------------------------------------
class C_indexlist {
  public:
    float x;
    int   i;
    int operator==(const C_indexlist &rhs) const;
    int operator<(const C_indexlist &rhs) const;  
};

//------------------------------------------------------------------------------
// matrix class
//------------------------------------------------------------------------------
class C_NxM {
  public:
    C_NxM(int,float);
    C_NxM(int,int,float);
    C_NxM(C_NxM &);  // 
    int    N; // first index
    int    M; // second index
    vector<vector<float> > x;
    void write(const string &);
};

//------------------------------------------------------------------------------
// int matrix class
//------------------------------------------------------------------------------
class C_NxMi {
  public:
    C_NxMi(int,int);
    C_NxMi(int,int,int);
    C_NxMi(C_NxMi &);  // 
    int    N; // first index
    int    M; // second index
    vector<vector<int> > x;
    void write(const string &);
};

//------------------------------------------------------------------------------
// stepset for a set number of breaks
//------------------------------------------------------------------------------
class C_stepset1 {
  public:
    int    nb; // number of breaks
    float  J; // chisquare estimate
    vector<int> b;
    C_stepset1 & operator=(const C_stepset1 &rhs);    
    int operator==(const C_stepset1 &rhs) const;
    int operator<(const C_stepset1 &rhs) const;  
    void print(const string &);
};



//------------------------------------------------------------------------------
// main steps class 
//------------------------------------------------------------------------------
class C_steps {
  friend ostream &operator<<(ostream &, const C_steps &);
  public:
    C_steps() {};
    C_steps(vector<float> & X,float STD,int MINBINS,int BMAX,int KBMAX, int DBG, int ALLOW);
    C_steps(vector<float> & X,float STD,int MINBINS,int BMAX,float FP,float SLOSH,int KBMAX, int DBG);
    //C_steps(vector<float> & X,float STD=0.25,int MINBINS=5,int BMAX=50,float FP=0.1,float SLOSH=0.05,int KBMAX=500) {
    // output values
    vector<int> b;      // best estimate of break points defined as starting point of new region
    vector<int> n;      // number of readcount bins for each region
    vector<float> cn;   // stepped CN for each breakpoint region
    vector<float> ax;   // average of x for each region
    vector<float> sx;   // stdev of x for each region
    vector<float> qx;   // quantized x for each region
    float chisq;        // fit chisq
    float Jmin;         // best fit J 
    int Nb;             // number of breaks for best fit
    // parameters
    float stdev;        // characteritic stdev of raw signal w/o CNV
    float mean;         // raw signal mean
    int minbins;        // no detected CNV can have fewer readcount bins
    int bmax;           // max number of breakpoints allowed
    //float fp;           // max fraction of bins allowed for candidate breakpoints
    //float chisqslosh;   // "close-enough" slosh for minimum chisquared
    int allow;
    int kbmax;          // maximum possible number of candidate breakpoints 
    int kbmin;          // minimum number of candidate breakpoints from prestep 
    float tdup;           // min dup x level for prestep kb selection 
    float tdel;           // max del x level for prestep kb selection 
    // internal values - kept public because I'm that sort of programmer...
    vector<float> x;    // raw signal: log10( read counts / expected counts) 
    vector<float> xb;   // signal space edges between quantized CN bins .. CN = 0,1,2,3...
    vector<float> xm;   // mean signal given fixed CN=0,1,2,3 ...
    vector<int> kb;     // filtered candidate breakpoints from prestep
    int L;              // number of x bins 
    int dbg;   
     // result of steplist()
    vector<C_stepset1> set; // stepsets for up to bmax breaks
    // inloglike vector
    vector<double> nloglike;
  //private:
    void defaultXB();
    void prestep();  // simpler way to get kb   
    void steplist();
    void stepjiggle();
    void stepset(int); 
    // function to return chisq for alternate b breakpoint vectors
    double stepset2chisq(vector <int> &);
    void calc_qx(vector <float> &);     
    float calc_qx1(float); 
    void trim(); 
};
#endif
