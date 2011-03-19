/*
 *  steps.cpp
 *  SpanDet
 *
 *  Created by Chip Stewart on 10/24/08.
 *  Copyright 2008 Boston College. All rights reserved.
 *
 */

#include "steps.h"

//------------------------------------------------------------------------------
// steps segmentation algorithm - translated from steps.R
//------------------------------------------------------------------------------
C_steps::C_steps(vector <float> & X, float STD,int MINBINS,int BMAX,int KBMAX, int DBG, int ALLOW) {
//C_steps::C_steps(vector <float> & X, float STD,int MINBINS,int BMAX,float FP,float SLOSH,int KBMAX, int DBG) {
/*
#
steps <- function (x,STD=0.4,SHORT=4,BMAX=20,XB=NULL,FP=0.1,LESS=0.05,KBMAX=50) {
#
# Find list of breaks in list of x's based on chisquare to stepped local average
# for all optimal intervals. Hypothesis that x points are normally distributed
# about local mean with a standard deviation of std (default 0.4). There can be up to 
# BMAX regions with BMAX-1  breaks of differents means. 
#
#  Input:  x list of values - log2 (count of reads/expected count of reads)
#			  STD:    assumed standard deviation (0.4)
#			SHORT:   	cost parameter for minimum break region size  (4)
#               no extra cost term if SHORT=0
#               also used to smooth x for calculating candidate break list
#			 BMAX:    maximum number of breaks per x (20)
#        XB:	 	bin edges in x space (low edges of bins starting with 1...)
#			   FP:    Fraction of points to make break candidate positions
#			 LESS:    Fraction less than minimum of chisquare for acceptable steps
#     KBMAX:	 	Max number of filtered breaks from prestep
#
#  Output: b:  	best estimate of break points defined as starting point of new region 
#         n:    number of SNP's for each region
#         x:    average of x for each region
#         s:    standard deviation for each region
#     chisq:    fit chisq
#		 chisqb:    fit chisq's for each level of breakpoints from 0 to BMAX
#
#  test: x=randn(1000); STD=0.4; CHISQL=0.5; SHORT=4; BMAX=20; XB=NULL; FP=0.1
#
*/
//------------------------------------------------------------------------------
	// total number of bins
  //cout << HUGE << endl;
  x=X;
  L=x.size();
	// fraction of L to allow candidate breaks  
  //fp=float(KBMAX)/L;
  //if (FP<fp) fp=FP;
  allow=ALLOW;
  // detection X signal stdev  - fixed 
  C_trimstat trimstats(x,0.05f);
  mean=trimstats.mean;
  if (STD>0.0) {
    stdev = STD;
  } else { 
    stdev = trimstats.std;
  }
  // Debug level 
  dbg =DBG;
  //
  // minium number of bins in a CNV region
  minbins = MINBINS;
  // maximum number of breaks in a chromosome 
  bmax = BMAX;
  // chisquare slosh
  //chisqslosh = SLOSH;
  // maximum number of candidate breakpoints
  kbmax = KBMAX;
  if (kbmax>L) kbmax=L;
  // minimum number of candidate breakpoints
  kbmin = KBMAX/10;
  if (kbmin>bmax) kbmax=bmax;
  // set xb quantized  x calibration scale 
  defaultXB();
  // set default prestep spike scale 
  tdel =log2(1.f/2.f);
  tdup =log2(3.f/2.f); 

  //------------------------------------------------------------------------------
  // list of candidate breaks	 
  //------------------------------------------------------------------------------
	prestep();
  //
  //vdump(kb,"prestep.kb.txt"); 
  
  //------------------------------------------------------------------------------
  // allow up to kb.size()-1 breaks
  //------------------------------------------------------------------------------
  bmax=kb.size()-1;

  //------------------------------------------------------------------------------
  // determine breaks for 1 to BMAX regions
  //------------------------------------------------------------------------------
	steplist();
  
  //------------------------------------------------------------------------------
  // jiggle breakpoint edges for best fit
  //------------------------------------------------------------------------------
  stepjiggle();
  
  //------------------------------------------------------------------------------
  // find the best break set
  //------------------------------------------------------------------------------
	stepset(-1);
}

//------------------------------------------------------------------------------
// function to estimate CN break values of log2 (normalized readcounts)
//------------------------------------------------------------------------------
void C_steps::defaultXB() {
  xb.clear();  // bin edges
  xm.clear();  // bin middles
  for (int i=0; i<101; i++) {
    float cnb = float(i)+0.5f;
    xb.push_back( float(log2(cnb/2)) );
    if (i==0) cnb=0.51f;
    xm.push_back( float(log2((cnb-0.5)/2)) );
  }
}

//------------------------------------------------------------------------------
// prestep2  filter low prob breakpoints for steps segmentation algorithm 
//------------------------------------------------------------------------------
void C_steps::prestep() {
// matlab function [kc,kb]=prebreaks(x1);
// peaks must be consistent with CNV levels
// 2/15/2009 CS
//----------------------------------------------------------------------------
	// count data
  // ? vector<int> n(L,0);
  // ? for (int i=0; i<L; i++)  	n[i]=i;
  // Smoothing coefficient
  int W = minbins;
	if ( (W%2)==0) {
    cerr << "need odd number of smoothing bins in steps " << minbins << endl;
    exit(-1);
  }
  // smoothing window w
  vector<float> w(W,float(1/float(W)));  
 	// Smooth data -> y
  vector<float> y=convolve(x, w);
      
  // peaks above CN~3 level
  float t3=tdup;
  list<int> pb3=peakBounds(y,t3,1);
  // peaks below CN~1 level
  float t1=tdel;
  list<int> pb1=peakBounds(y,t1,-1);
  list<int> pb;
  //
  if (dbg>0) {
    if ((dbg&1)>0) vdump(y,"prestep.y.txt");
    if ((dbg&2)>0) ldump(pb1,"prestep.pb1.txt"); 
    if ((dbg&4)>0) ldump(pb3,"prestep.pb3.txt"); 
  }
  //
  pb.merge(pb1);
  //
  if ((dbg>0)&&((dbg&8)>0)) ldump(pb,"prestep.pb1a.txt"); 
  //
  pb.merge(pb3);
  //
  if ((dbg>0)&&((dbg&16)>0))  ldump(pb,"prestep.pb3a.txt"); 
  //
  pb.sort();
  pb.unique();
  //
  if ((dbg>0)&&((dbg&32)>0))  ldump(pb,"prestep.pbsorted.txt"); 
  //
  for (list<int>::iterator ipb=pb.begin(); ipb!=pb.end(); ++ipb) {
    kb.push_back(*ipb);
  }
  if ((dbg>0)&&((dbg&64)>0)) {
    vdump(x,"prestep.x.txt"); 
    vdump(y,"prestep.y.txt"); 
    vdump(kb,"prestep.kb.txt"); 
  }
}  

//-------------------------------------------------------------------------
// calculates candidate break points around peaks
//-------------------------------------------------------------------------
list <int> peakBounds(vector <float> & y, float t, int flip) {
  // flip*y peaks above t level 
  // ... "flip" allows peak and trough bounding with the same code....
  vector<int> pp;  // peak positions
  vector<float> pv;  // peak values
  size_t L = y.size();
  float y1=0;
  int i1;
  for (int i=0; i<int(L); i++) { 	
     if ((flip*y[i])>(flip*t)) {
        pp.push_back(i);
        pv.push_back(flip*y[i]);
     }
     if ((flip*y[i])>y1) {
        y1=flip*y[i];
        i1=i;
     }
  }  
  //
  // use biggest peak if none are over threshold
  if (pp.size()<1)  {      
        pp.push_back(i1);
        pv.push_back(flip*y[i1]);
  }
  // cluster by contiguous k's
  size_t Np=pp.size();
  list<int> group;  // group index for peaks 
  if (Np==0) {      // if there are no peaks, then return empty list
     return group;
  }
  group.push_back(1);
  int dp=0, gap=0, pg=1;
  for (int i=0; i<int(Np-1); i++)  {
      dp=pp[i+1]-pp[i];          // distance to next peak value
      gap=(dp>2? 1: 0);         // mark gaps dx>2
      pg+=gap;
      group.push_back(pg);  // bump group index if gap
  }
  list<int> gu=group;
  gu.unique();              // list of unique groups
  //
  /*
  vdump(y,"prestep.y.txt"); 
  vdump(pp,"prestep.pp.txt"); 
  vdump(pv,"prestep.pv.txt"); 
  ldump(group,"prestep.group.txt");
  ldump(gu,"prestep.gu.txt"); 
  */
  //
  list<int> ppc;            // "central" index for each group 
  list<int> ppb;            // several "breaks" indices for each group 
  // loop over groups
  for (list<int>::iterator ig=gu.begin(); ig!=gu.end(); ++ig) {
    float ymx = 0;
    int   pmx = 0;
    int gu1=*ig;
    int i=0; // index into pp pv and group
    // loop over pp & pv within group to find "central" pp
    for (list<int>::iterator im=group.begin(); im!=group.end(); ++im) {
      i++;
      int g1=*im;
      if (g1<gu1) continue;
      if (g1>gu1) break;
      if (pv[i-1]>ymx) {
         ymx=pv[i-1]; // start index at 0
         pmx=pp[i-1];
      }
    }
    // center of peak
    ppc.push_back(pmx);
    // find peak low boundary
    int p1=pmx, p2=0;
    while (p1>1) {
        p1--;
        if ((flip*y[p1])<(0.5*(flip*y[pmx]))) p2=(p2==0? p1 : p2);
        if ((flip*y[p1])<0) {
          ppb.push_back(p1);
          ppb.push_back(p2);
          break;
        }
    }
    // find peak upper boundary
    p1=pmx, p2=0;
    while (p1<int(L)) {
        if ((flip*y[p1])<(0.5*(flip*y[pmx]))) p2=(p2==0? p1 : p2);
        if ((flip*y[p1])<0) {
          ppb.push_back(p1);
          ppb.push_back(p2);
          break;
        }
        p1++;
    }
  }
  ppb.sort();
  ppb.unique();
  return ppb;
}


//------------------------------------------------------------------------------
// steplist main ugly crunching code for steps
//------------------------------------------------------------------------------
void C_steps::steplist() {
/*
# Find list of breaks in list of x's based on chisquare to quantized local average
# for every possible interval. Hypothesis that x points are normally distributed
# about quantized local mean with a standard deviation of std (default 0.4). 
# There can be up to BMAX regions with BMAX-1 breaks of different quantized means. 
# The quantized mean is quantized to midpoints of binned average response. 
# 
#  Input:  x: list of values
#         KB: candidate break points
#        STD: assumed standard deviation (0.4)
#		 MINBINS: cost parameter for minimum break region size  (4)
#			  BMAX: maximum number of breaks   (20)
#         XB:	bin edges in x space (low edges of bins starting with 1...)
#
#  Output: J: list of chisq values for each hypothetical number of breaks (1-20)
#          b: matrix of break points for each hypothetical number of breaks (1-20)
#
#  test:  STD=0.4; MINBINS=5; BMAX=20; XB=NULL;
#   KB=NULL
*/
//------------------------------------------------------------------------------
	// candidate break points ... use them all if no kb specified
	if  (kb.size()<1) {
    cerr << " prestep failure - use spaced kb " << endl;
    for (int i=2; i<(L-1); i+=minbins) {
      kb.push_back(i);
    }
  }
  
  // remove ends of kb if they include the start or end of the chromosome
  // or even close to ends (within minbins)
  int ikb=0;
  while (kb[0]<minbins) {
     kb.erase(kb.begin());
     ikb++;
     if (ikb>=minbins) {
       cerr << "runaway bin trimming kb begin in steplist " << endl;
       break;
     }
  }
  ikb=0;
  while (kb[kb.size()-1]>(L-minbins-1)) {
     kb.pop_back();
     ikb++;
     if (ikb>=minbins) {
       cerr << "runaway bin trimming kb end in steplist " << endl;
       break;
     }
  }
  //if (kb[0]==0) kb.erase(kb.begin());
  //if (kb[kb.size()-1]==L-1) kb.pop_back();
  
  // check x scale bin edges & mid-points
  if (xb.size()==0)  {
    	defaultXB();
  }
  if ((dbg>0)&((dbg&128)>0)) {
    vdump(x,"steplist.x.txt"); 
    vdump(kb,"steplist.kb.txt"); 
    vdump(xb,"steplist.xb.txt"); 
    vdump(xm,"steplist.xm.txt");
  }
  
  // number of candidate breakpoints - including 0
  int Nk=kb.size()+1;
  // k1 region starts
  vector <int> k1(1,0);
  k1[0]=0;
  for (int i=0; i<(Nk-1); i++) {
      k1.push_back(kb[i]);
  }
  // k2 region ends
  vector <int> k2(kb);
  for (int i=0; i<(Nk-1); i++) k2[i]--;
 
  k2.push_back(L-1);
	// cumulative sum x and x**2 
  vector <float> x1(L,0);
  vector <float> x2(L,0);
  for (int i=0; i<L; i++) {
    x1[i]=x[i];
    x2[i]=x[i]*x[i];
    if (i>0) {
      x1[i]+=x1[i-1];
      x2[i]+=x2[i-1];
    }
  }  	
	// running bin count
  //vector<float> n(L,0);
  // for (int i=0; i<L; i++)  	this->n[i]=i+1;
  // cumulative running average
  vector<float> xn(L,0);  
  for (int i=0; i<L; i++)  	xn[i]=x1[i]/(i+1);
  // quantized cumulative running average
	calc_qx(xn);  	
  /*
  vdump(x1,"steplist.x1.txt"); 
  vdump(x2,"steplist.x2.txt"); 
  vdump(n,"steplist.n.txt"); 
  vdump(xn,"steplist.xn.txt"); 
  vdump(qx,"steplist.qx.txt"); 
  */
	// declare X and reserve array space for all possible intervals
	C_NxM X(Nk,HUGE);	
	// X starting at pos = 1, ending at k2
  for (int i=0; i<Nk; i++) {
    X.x[0][i] = x2[k2[i]] -2*x1[k2[i]]*qx[k2[i]] + k2[i]*(qx[k2[i]]*qx[k2[i]]);
    if (minbins>0) X.x[0][i]+=float(2*log(1+exp(minbins-k2[i])));
  }

  //X.write("steplist.X.x.txt"); 
	
	// start at k1[k]+1, end at kb[k]
	for (int k=1; k<Nk; k++) { 
    	int i = k1[k]-1;       // last point in preceeding region 
      vector<float> X2(Nk,0);
      vector<float> X1(Nk,0);
      vector<int> np(Nk,0);
      vector<float> xna(Nk,0);
      for (int j=0; j<Nk; j++) {     	
        X2[j]=x2[k2[j]]-x2[i];	// each segment ending at k2 - everything up to k1[k]
     		X1[j]=x1[k2[j]]-x1[i];  
        np[j] = k2[j]-i;
        if (np[j]==0) np[j]=-1;   // avoid /0
        xna[j]=X1[j]/np[j];
      }
      calc_qx(xna);
      for (int j=0; j<Nk; j++) {     	
        X.x[k][j] = X2[j] -2*X1[j]*qx[j] + np[j]*(qx[j]*qx[j]);
        if (minbins>0) X.x[k][j]+=float(2*log(1+exp(minbins-np[j])));
   	  	// set nonsensical intervals to infinite variance
	  	  if (np[j]<1) X.x[k][j] = HUGE;
      }
	}
	// normalize for chisquare[i,j] of (x-qx)^2/sigma from kb[i] to kb[j] 
  C_NxM chisq(X);  
  for (int i=0; i<Nk; i++) {     	
    for (int j=0; j<Nk; j++) {     	
      chisq.x[i][j] = X.x[i][j]/(stdev*stdev);
    }
  }	
  if ((dbg>0)&((dbg&256)>0)) {
    X.write("steplist.X.x.txt"); 
    chisq.write("steplist.chisq.txt"); 
	}
  
  // max allowed breakpoints
	int Km = bmax;
	// declare Jn (candidate break set chisq) and p (positions)
	C_NxM Jn(Km,Nk,HUGE);
	C_NxMi p(Km-1,Nk,0);
	// Jn with no breaks...
  for (int i=0; i<Nk; i++) {     	
   	Jn.x[0][i] = chisq.x[0][i];
  }
  // loop from 1 to Km-1 candidate break regions
	for (int k=1; k<(Km-1); k++) {
	   // loop to regions up to N
	   for (int j=k; j<Nk; j++) {
        int p1 = 0;
        for (int i=0; i<=j-1; i++) {  // i<j-1 ??
          float q = Jn.x[k-1][i]+chisq.x[i+1][j];
          if (q<=Jn.x[k][j]) {
            Jn.x[k][j] = q*(1.0f+1e-10f);
            p1=i;
          }
        }
        if (Jn.x[k][j]<HUGE) {
  	  	 	// last break point at minimum chisquare
         	p.x[k-1][j] = p1;
			  }
	   }
	}
	
  if ((dbg>0)&((dbg&512)>0)) {
    Jn.write("steplist.Jn.txt"); 
    p.write("steplist.p.txt"); 
  }
  
  // calculate Km max
  int p1 = 0;
  for (int i=0; i<Nk-1; i++) {
      if ((Jn.x[Km-2][i]+chisq.x[i+1][Nk-1])<=Jn.x[Km-1][Nk-1]) {
          Jn.x[Km-1][Nk-1] = Jn.x[Km-2][i]+chisq.x[i+1][Nk-1];
          p1=i;
      }
  }
  if (Jn.x[Km-1][Nk-1]<HUGE) {
      // last break point at minimum chisquare
      p.x[Km-2][Nk-1] = p1;
  }
  
  /*
  Jn.write("steplist.Jn.2.txt"); 
	p.write("steplist.p.2.txt"); 
	*/
  
  // best chisquare for each k
  C_NxMi bp(Km,-1);
  for (int i=0; i<Km;i++) {
    bp.x[i][i]=Nk-1;
  }
  for (int i=0; i<Km;i++) {
    C_stepset1 set1;
    set1.J=Jn.x[i][Nk-1];
    set.push_back(set1);
  }
  //set boundary positions
	for (int K=1; K<Km; K++) {
    for (int k1=K-1; k1>=0; k1--) {
    		bp.x[K][k1] = p.x[k1][bp.x[K][k1+1]];
    }
    for (int k1=0; k1<Km; k1++) {
      if (bp.x[K][k1]>=0) {
        set[K].b.push_back(k2[bp.x[K][k1]]);
      }
    }
    // trim off last element of b if b[b.size()-1] == x.size()=L-1
    if (set[K].b.size()>0) {
      if ( set[K].b[set[K].b.size()-1]==(L-1) ) {
        set[K].b.pop_back();
      }
    }
    /*
    bp.write("steplist.bp.txt"); 
    stringstream strm;
    strm << K;
    string sk = strm.str();
    set[K].print("steplist.set."+sk+".txt");
    */
  }
}

//------------------------------------------------------------------------------ 
// find best set[k].J and calculate fitted parameters & levels for breakset k
//------------------------------------------------------------------------------
void C_steps::stepjiggle() {
  // range of jiggling
  int BJ=int(round(float(minbins)*2.0f));
  // small event limit 
  int DB=int(round(float(minbins)*0.5f));
  //---------------------------------------------------------------------------- 
  // loop over all step sets - increasing by number of breakpoints
  //---------------------------------------------------------------------------- 
  for (int k=0; k<int(set.size()); k++) {
    // set member data for this set of breaks b, xa,qx,
    stepset(k);
    // original chisq 
    double chisq0=set[k].J;
    // best chisq so far
    double chisq1=HUGE;
    //---------------------------------------------------------------------------- 
    // loop over each break
    //----------------------------------------------------------------------------     
    for (int q1=0; q1<k; q1++) {
      // original break
      int b0=b[q1];
      // preceding breakpoint
      int bbefore=(q1>0? b[q1-1]: 0);
      // following breakpoint
      int bafter=(q1<(k-1)? b[q1+1]: L);
      // save best
      int best=b0;
      // loop over jiggle vector of breakpoints
      for (int db=-BJ; db<=BJ; db++) {
        // test break 
        int bt=b0+db;
        // not too close to previous break
        if (bt<(bbefore+DB)) { continue; }
        // not too close to next break
        if (bt>(bafter-DB)) { continue; }
        // try this break
        b[q1]=bt;
        double chisqt=stepset2chisq(b);
        // was this an improvement? 
        if (chisqt<chisq1) {
           chisq1 = chisqt;
           best = bt;
        }
      }
      // was this an overall improvement? 
      if (chisq1<chisq0) {
        b[q1]=best;
        set[k].b=b;
        stepset(k);
        chisq0=chisq1;
        set[k].J=chisq1;
      } else {
        b = set[k].b;
      }
    }
  }
}

double C_steps::stepset2chisq(vector<int> & b1) {
    // calculate chisq for this particular set of breaks b1 
    int k=b1.size();
    // local vectors
    vector<int> n1(k+1,0);
    vector <double> ax1(k+1,0);
    vector <double> sx1(k+1,0);
    vector <double> qx1(k+1,0);
    vector <double> cs1(k+1,0);
    int pp=0;
    int p1,p2; 
    // main chisqr
    double chisq=0;       
    // loop over all breaks to calculate local chisq
	  for (size_t i =0; i<=b1.size(); i++)  {
      p1 = pp;
      if (i<b1.size()) {
        p2 = b1[i];
      } else  { 
        p2 = L;
      }
      pp = p2;
      int na = p2-p1;
      float x1=0;
      float x2=0;
      for (int p=p1; p<p2; p++) {
        x1+=x[p];
        x2+=x[p]*x[p];
      }
      //local ave
      x1=x1/na;
      // local stdev
      float xs1=float(sqrt(fabs(x2/na-x1*x1)));
      n1[i]=na;
      ax1[i]=x1;
      sx1[i]=xs1;
      // original quantized levels
      qx1[i]=calc_qx1(x1);
      cs1[i]=0;
      for (int p=p1; p<p2; p++) {
        double c1=(x[p]-qx[i])/stdev;
        c1*=c1;
        cs1[i]+=c1;
      }
      chisq+=cs1[i];
    }
  return chisq;
}


//------------------------------------------------------------------------------ 
// find best set[k].J and calculate fitted parameters & levels for breakset k
//------------------------------------------------------------------------------
void C_steps::stepset(int k1) {
  //---------------------------------------------------------------------------- 
  // scan for best set[k].J
  //---------------------------------------------------------------------------- 
  int k;
  if (k1<0) {
    Jmin = HUGE;
    Nb=0;
    // scan for minimum J
    for (k=0; k<bmax; k++) {
      if (set[k].J<Jmin) {
         Jmin = set[k].J;
         Nb=k;
      }
    }
    
    //---------------------------------------------------------------------------- 
    // renormalize stdev such that the min chisq is equal to dof
    // then find minimum Nb with chisq negative log likelihood > 1 (about 1/e probable)  
    //----------------------------------------------------------------------------     
    // dof for Jmin
    int dofL=L-Nb;
    // effective stdev normalizing chi2 for  overfitted CNV function: 1=chi2./dof
    double stdev0=stdev;
    stdev=sqrt(Jmin/dofL)*stdev;
    // init nloglike vector
    nloglike.clear();
    double nloglike1=HUGE;
    // renormalize all chisq
    for (k=0; k<bmax; k++) {
      set[k].J=set[k].J*(stdev0/stdev)*(stdev0/stdev);
      // dof for this Nb
      double dof=L-k;
      // running chisq prob that data fits model 
      double pr = 1-pnorm( set[k].J,dof,sqrt(2*(dof)));
      // negative log like 
      nloglike.push_back(-log(pr));
      // find the lowest nloglike - should be ~0
      if (nloglike[k]<nloglike1) {
         nloglike1=nloglike[k];
         Nb=k;
      }
    }
    
    //-------------------------------------------------------------------------- 
    // allow more slosh to get minimum set of breakpoints 
    //-------------------------------------------------------------------------- 
    double Allow=(0.4f)+nloglike1;
    if (Nb<allow) {
      Allow=(1e-10)+nloglike1;
    }        

    //-------------------------------------------------------------------------- 
    // scan nloglike for lowest nloglike within 1 unit of global min
    //-------------------------------------------------------------------------- 
    for (k=0; k<bmax; k++) {
       //   pick the lowest within 1 loglike of min
       if (nloglike[k]<=Allow) {
         Jmin = set[k].J;
         Nb=k;
         break;
       }
    }
  } else {
    k=k1;
    Nb=k;
  }
  //---------------------------------------------------------------------------- 
  // calculate stuff for the best break set k
  //---------------------------------------------------------------------------- 
  b=set[k].b;
  n.clear();
  ax.clear();
  sx.clear();
  int pp=0;
  int p1,p2; 
	for (size_t i =0; i<=b.size(); i++)  {
    p1 = pp;
    if (i<b.size()) {
      p2 = b[i];
      pp=b[i];
    } else  { 
      p2 = L;
    }
    pp=p2;
    int na = p2-p1;
    float x1=0;
    float x2=0;
    for (int p=p1; p<p2; p++) {
      x1+=x[p];
      x2+=x[p]*x[p];
    }
    //local ave
    x1=x1/na;
    // local stdev
    float xs=float(sqrt(fabs(x2/na-x1*x1)));
    n.push_back(na);
    ax.push_back(x1);
    sx.push_back(xs);
	}
  calc_qx(ax);
  
  if ((dbg>0)&((dbg&1024)>0)) {
    vdump(b,"stepset.b.txt");
    vdump(x,"stepset.x.txt");
    vdump(n,"stepset.n.txt");
    vdump(ax,"stepset.ax.txt");
    vdump(sx,"stepset.sx.txt");
    vdump(sx,"stepset.qx.txt");
    vdump(sx,"stepset.cn.txt");
    set[k].print("steplist.finalset.txt"); 
  }
}

float C_steps::calc_qx1(float x1) {
  //----------------------------------------------------------------------------
  // quantized copy number function given XB break positions in X
  //----------------------------------------------------------------------------
  int b1=0;
  while (x1>=xb[b1]) {
      b1++;
      if (b1==int(xb.size())) break;
  }
  // bloated duplication trim to 100x
  if ((b1+1)>int(xm.size())) {
    b1=xm.size()-1;
  }
  return xm[b1];
}  

void C_steps::calc_qx(vector <float> & xn) {
  //----------------------------------------------------------------------------
  // quantized copy number function given XB break positions in X
  //----------------------------------------------------------------------------
	if (xb.size()==0) qx=xn;
  if (xb.size()==1) {
    if (xb[0]==0)  qx=xn;
  }
  qx.clear();
  cn.clear();
  // counter N  - is there a better way to do this? 
  vector <int> N;
  for (size_t i=1; i<xb.size(); i++)  N.push_back(i);
  
  size_t Nb=xb.size();
  if (xm.size()!=Nb) {
     cerr << " problem in qx with xm " << endl;
     exit(-1);
  }
  //----------------------------------------------------------------------------
  // loop over xn 
  //----------------------------------------------------------------------------
  int n = xn.size();
  for (int i=0; i<n; i++) {
      int b1=0;
      while (xn[i]>=xb[b1]) {
        b1++;
        if (b1==int(Nb)) break;
      }
      // bloated duplication trim to 100x
      if ((b1+1)>int(xm.size())) {
         b1=xm.size()-1;
      }
      //------------------------------------------------------------------------
      // quantized mean response for this level
      //------------------------------------------------------------------------
      qx.push_back(xm[b1]);
      //------------------------------------------------------------------------
      // corresponding copy number
      //------------------------------------------------------------------------
      cn.push_back(b1);
  } 
  // set qx==0 to a small number
  /*
  for (int i=0; i<n; i++) {
      if (qx[i]<1.0) qx[i]=float(0.1);
  }
  */
}  

C_stepset1 & C_stepset1::operator=(const C_stepset1 &rhs) { 
   this->nb=rhs.nb; // number of breaks
   this->J=rhs.nb; // chisquare estimate
   this->b=rhs.b;
   return *this;
}
  
void C_stepset1::print(const string & f) {
 // open output file. bomb if unable to open
  fstream output(f.c_str(), ios::out );
  if (!output) {
      cerr << "Unable to open file: " << f << endl;
      return;
  }
  output << f << " " << J << endl;
  for (size_t i=0; i<b.size(); i++)  {
    output << b[i] << endl;
  }
  output.close(); 
}

//------------------------------------------------------------------------------
// function to calcuate trimmed mean 
// trimmed mean with trim/2 tails removed on both sides 
// 0 < trim < 1
//------------------------------------------------------------------------------
float meantrimmed(vector<float> & X,float trim) { 
   int N=X.size();
   list<float> x;
   for (size_t i=0; i<X.size(); i++) {
      x.push_back(X[i]);
   }
   int Ntail = int(N*trim/2.0);
   list<float>::iterator it,it1,it2;
   x.sort();
   /*
   for (it=x.begin(); it!=x.end(); ++it) { 
     cout << (*it) << endl;
   }
   */
   int i=0;
   it1=x.begin();
   it2=x.end();
   while (i<Ntail) {
      it1++;
      it2--;
      i++;
   }
   //x.erase(it2,x.end());
   //x.erase(x.begin(),it1);
   float a=0;
   float n=0;
   for (it=it1; it!=it2; ++it) { 
    a+=(*it);
    n++;
   }
   a = a/n;
   return a;
}

//------------------------------------------------------------------------------
// class to provide trimmed mean std for vector X
// 0 < trim < 1
//------------------------------------------------------------------------------
C_trimstat::C_trimstat(vector<float> & X,float trim) {
   N0=X.size();
   list<float> x;
   for (size_t i=0; i<X.size(); i++) {
      x.push_back(X[i]);
   }
   int Ntail = int(N0*trim/2.0);
   list<float>::iterator it,it1,it2;
   x.sort();
   int i=0;
   it1=x.begin();
   it2=x.end();
   while (i<Ntail) {
      it1++;
      it2--;
      i++;
   }
   N=0;
   double x1=0;
   double x2=0;
   for (it=it1; it!=it2; ++it) { 
      x1+=(*it);
      x2+=(*it)*(*it);
      N++;
   }
   mean = x1/N;
   std = sqrt(fabs(x2/N-(mean*mean)));
   x2=0;
   for (it=it1; it!=it2; ++it) { 
      x2+=(*it);
      if (x2>(x1/2.0f)) {
        median=(*it);
        break;
      }
   } 
}


//------------------------------------------------------------------------------
// simple convolution function for smoothing ... no phase shift
//------------------------------------------------------------------------------
vector <float> convolve(vector <float> & a, vector <float> & tap) {
// tap count
  int ntap=tap.size();
  if (ntap%2==0) {
   cerr << " use odd tap length in convolve" << endl;
   exit(-1);
  }
  //reserve output vector of same size as a (my convention)
  vector <float> x(a.size(),0);
  // side taps
  int W=(ntap-1)/2;
  // zero pad 
  vector <float>  z(W,0);
  // padded a
  vector <float> x0(z);
  x0.insert (x0.end(),a.begin(),a.end());
  x0.insert (x0.end(),z.begin(),z.end());  
  for (size_t i=0; i< a.size(); i++) {
    for (size_t j=0; j< tap.size(); j++) {
       x[i]+=x0[i+j]*tap[ntap-1-j];
    }
  }
  return x;
}

//------------------------------------------------------------------------------
// indexed list for sorting things with order recorded
//------------------------------------------------------------------------------
int C_indexlist::operator==(const C_indexlist &rhs) const
{
   if( this->x != rhs.x) return 0;
   if( this->i != rhs.i) return 0;
   return 1;
}

// This function is required for built-in STL list functions like sort
int C_indexlist::operator<(const C_indexlist &rhs) const
{
   if( this->x < rhs.x ) return 1;
   if( this->x == rhs.x && this->i < rhs.i ) return 1;
   return 0;
}



//------------------------------------------------------------------------------
// inititalize NxM matrix
//------------------------------------------------------------------------------
C_NxM::C_NxM(int N1, float x0) {
    N=N1; // number of breaks
    M=N1; // number of breaks
    x.clear();
    vector<float> x1(N,x0);
    for (int i=0; i<N; i++) {
      x.push_back(x1);
    };
}
//------------------------------------------------------------------------------
// inititalize NxM matrix
//------------------------------------------------------------------------------
C_NxM::C_NxM(int N1, int M1, float x0) {
    N=N1; // number of breaks
    M=M1; // number of breaks
    x.clear();
    vector<float> x1(M,x0);
    for (int i=0; i<N; i++) {
      x.push_back(x1);
    };
}

//------------------------------------------------------------------------------
// copy constructor NxM matrix
//------------------------------------------------------------------------------
C_NxM::C_NxM(C_NxM & x0) {
    N=x0.N; // number of breaks
    M=x0.M; // number of breaks
    x.clear();
    for (int i=0; i<N; i++) {
      x.push_back(x0.x[i]);
    };
}

void C_NxM::write(const string & f) {
 // open output file. bomb if unable to open
  fstream output(f.c_str(), ios::out );
  if (!output) {
      cerr << "Unable to open file: " << f << endl;
      return;
  }
  output << f << endl;
  for (int i=0; i<N; i++)  {
    for (int j=0; j<M; j++)  {     
      output << x[i][j] << " ";
    }
    output << endl;
  }
  output.close(); 
}


//------------------------------------------------------------------------------
// inititalize NxMi matrix
//------------------------------------------------------------------------------
C_NxMi::C_NxMi(int N1, int x0) {
    N=N1; // number of breaks
    M=N1; // number of breaks
    x.clear();
    vector<int> x1(N,x0);
    for (int i=0; i<N; i++) {
      x.push_back(x1);
    };
}
//------------------------------------------------------------------------------
// inititalize NxM matrix
//------------------------------------------------------------------------------
C_NxMi::C_NxMi(int N1, int M1, int x0) {
    N=N1; // number of breaks
    M=M1; // number of breaks
    x.clear();
    vector<int> x1(M,x0);
    for (int i=0; i<N; i++) {
      x.push_back(x1);
    };
}

//------------------------------------------------------------------------------
// copy constructor NxM matrix
//------------------------------------------------------------------------------
C_NxMi::C_NxMi(C_NxMi & x0) {
    N=x0.N; // number of breaks
    M=x0.M; // number of breaks
    x.clear();
    for (int i=0; i<N; i++) {
      x.push_back(x0.x[i]);
    };
}
void C_NxMi::write(const string & f) {
 // open output file. bomb if unable to open
  fstream output(f.c_str(), ios::out );
  if (!output) {
      cerr << "Unable to open file: " << f << endl;
      return;
  }
  output << f << endl;
  for (int i=0; i<N; i++)  {
    for (int j=0; j<M; j++)  {     
      output << x[i][j] << " ";
    }
    output << endl;
  }
  output.close(); 
}

template< typename vType>
void vdump(vector <vType> & v1, const string & f) {
 // open output file. bomb if unable to open
  fstream output(f.c_str(), ios::out );
  if (!output) {
      cerr << "Unable to open file: " << f << endl;
      return;
  }
  output << f << endl;
  for (size_t i=0; i<v1.size(); i++)  output << v1[i] << endl;
  output.close(); 
}

template< typename wType>
void ldump(list<wType> & w1, const string & f) {
 // open output file. bomb if unable to open
  fstream output(f.c_str(), ios::out );
  if (!output) {
      cerr << "Unable to open file: " << f << endl;
      return;
  }
  output << f << endl;
  typename list<wType>::iterator i;
  for( i = w1.begin(); i != w1.end(); i++ ) output << *i << endl;
  output.close(); 
}

//------------------------------------------------------------------------------
// constructor trim to remove outliers from vector x
//------------------------------------------------------------------------------
/*  trim  <- function (x,STD=0.4,W=3,CUT=2) {
# trim single points from list of x with large residuals 
# from local mean (excluding this point)
# INPUT: 	x 	data list to be trimmed
#       STD 	expected standard deviation of x
#     		W 	length of averaging window
#       CUT		number of std deviations to trim from x
#
# OUTPUT:	structure
#         x		trimmed list  (size equal or smaller than original x list)
#         m   indices of x  (1:N with holes for cut x)
#        xx   list of removed original x
#        mx   indices of removed original x
#
#  debug x = x1; std=0.4; window=3; cut=2
*/
C_trim::C_trim(vector<float> & X,float STD,int W, int CUT) { 
  std=STD;
  w=W;
  cut=CUT;
  // first remove all -infs and +infs
  x.clear();  //		trimmed list  (size equal or smaller than original x list)
  m.clear();  //    indices of x  (1:N with holes for cut x)
  xx.clear();	//  	list of removed original x
  mx.clear(); //   	indices of removed original x
  vector <float> x0;
  vector <int> m0;
  for (int i=0; i< int(X.size()); i++) {
    bool pass= fabs(X[i])<1e30;
    if (pass) {
      m.push_back(i);
      x.push_back(X[i]);
    } else {
      mx.push_back(i);
      xx.push_back(X[i]);
    }
  }
  /*
  // tap smoothing vector in window of length 2W+1: k = 0*c(1:(2*W+1))+1/(2*W)
  vector <float>  tap(2*W+1,1/float(2*W));
  tap[W]=0;  // center element zeroed out: k[W+1]=0
  // filter:     x2=convolve(x0, k, type = "filter")
  vector <float> x2 = convolve(x0,tap); 
  // absolute deviation:   r2=abs(x-x2)
  vector <float> r2(x0.size(),0);
  for (size_t i=0; i< x2.size(); i++) {
    r2[i]=float(fabs(x0[i]-x2[i]));
  }
  vdump(x0,"x0.txt"); 
  vdump(x2,"x2.txt"); 
  vdump(r2,"r2.txt"); 
  for (int i=0; i< int(x0.size()); i++) {
    bool pass=r2[i]<float(CUT*STD);
    if (pass) {
      m.push_back(m0[i]);
      x.push_back(x0[i]);
    } else {
      mx.push_back(m0[i]);
      xx.push_back(x0[i]);
    }
  }
  vdump(x,"pass.x.txt"); 
  vdump(m,"pass.m.txt"); 
  vdump(xx,"trim.x.txt"); 
  vdump(mx,"trim.m.txt"); 
  */
  //sort(mx.begin(),mx.end());
  //sort(xx.begin(),xx.end());
}  	

double pnorm(const double x0, const double mu, const double sigma)
{
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;
  
  double x = (x0-mu)/sigma;

  if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  else {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}




