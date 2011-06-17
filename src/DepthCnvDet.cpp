/*
 *  DepthCnvDet.cpp
 *  SpanDet
 *
 *  Created by Chip Stewart on 10/17/08.
 *  Copyright 2008 Boston College. All rights reserved.
 *
 */

#include "DepthCnvDet.h"


C_DepthCNV::C_DepthCNV(C_pairedfiles & data,  RunControlParameters & pars) {
  C_contigs::const_iterator iterContig;   
  int Nset = data.set.size();
  // debug level  
  int dbg=-pars.getDbg();
  // max gap of measured bin in one CNV
  GAP=pars.getRDmaxGap(); // 100;    
  // loop over sets...
  for (int iset = 0; iset< Nset; iset++) {
    string setName = data.set[iset].getSetName();
    // loop over chromosomes/contigs
    for (iterContig = data.set[iset].contig.begin();	iterContig != data.set[iset].contig.end(); iterContig++) {
      // determine set/chromosome basename for this iteration
      string name = iterContig->first;
      string area = pars.getOutputDir();
      int k = 1+area.find_last_of("/");
      string subdir = area.substr(k);
      if (area.size()>0) { area = area+"/";}
      string prefix = pars.getPrefix();
      if (prefix.size()>0) { prefix = prefix+".";}
      string sn1 = setName;
      // don't need redundant setName if subdirectory is already the setName
      //if (sn1==subdir) sn1="";
      //if (sn1.size()>0) { sn1 = sn1+".";}
      string basename = area+prefix+sn1+"."+name;
      if (setName==name) basename = area+prefix+sn1;
    
      // replace evil character "|" with benign "_" 
      size_t found=basename.find("|");
      while (found!=string::npos) {
        basename.replace(found,1,"_");
        found=basename.find("|");
      }
      
      //-----------------------------------------------------------------------
      // Get unique coverage
      //-----------------------------------------------------------------------
      C_UniqueCoverage u1(data.set[iset].contig[name],  pars);
      u=u1;
      u.write(data.set[iset].contig[name],basename);
      //-----------------------------------------------------------------------
      // threshold  count  e.n at 0.20 
      //-----------------------------------------------------------------------
      // float low=0.2f;
      float etoolow = u.e.Stats.h.mode*pars.getRDminExp();   //u.e.Stats.h.p2xTrim(low);
      vector <float> x(u.r.n.size(),HUGE);
      float nx=0,sx=0,sxx=0;
      //-----------------------------------------------------------------------
      // remove poorly measured bins...
      //-----------------------------------------------------------------------
      for (size_t i=0; i<u.r.n.size(); i++) {
         if (u.e.n[i]>etoolow)  {
            if  (u.r.n[i]>0) {
              x[i]=float(log2(fabs(u.r.n[i]/u.e.n[i]))) - u.gcf.n[i];
              if (fabs(x[i])<2.0f) {
                nx++;
                sx+=x[i];
                sxx+=x[i]*x[i];
              }
            } else { 
              x[i]=-3.0;
            }
          }
      }
      float xa=sx/nx;
      float S1=sqrt( (sxx/nx)-(xa*xa)); //0.4f;
      // recenter on zero      
      for (size_t i=0; i<u.r.n.size(); i++) {
         if (x[i]>-2.99) {
            x[i]=x[i]-xa;
         }
      }
      int W1=4;
      int CUT1=4;
      
      //-----------------------------------------------------------------------
      // trimmed list based on x
      //-----------------------------------------------------------------------
      C_trim xt(x,S1,W1,CUT1);
      mx=xt.m;
      
      //-----------------------------------------------------------------------
      // these parameters should be in pars but for now just set them here
      //-----------------------------------------------------------------------
      float STD=0.0f;
      MINBINS=pars.getRDminbin();    // 5
      BMAX=150;
      //float FP=1000.1f;
      //float SLOSH=0.01f;
      KBMAX=500;
      ALLOW=pars.getRDallow();
      //------------------------------------------------------------------------ 
      // steps segmentation algorithm
      //------------------------------------------------------------------------
      C_steps steps1(xt.x,STD,MINBINS,BMAX,KBMAX, dbg,ALLOW);
      this->steps = steps1;
      //------------------------------------------------------------------------
      // back to untrimmed coordinates
      //------------------------------------------------------------------------
      int nb = steps.b.size();
      bb.clear();
      for (int ib = 0; ib<nb; ib++) {
         bb.push_back(xt.m[steps.b[ib]]);
      }

      //------------------------------------------------------------------------
      //  dump raw steps results
      //------------------------------------------------------------------------
      /*
      vdump(x,"steps.x.txt"); 
      vdump(bb,"steps.bb.txt"); 
      vdump(xt.x,"steps.xt.x.txt"); 
      vdump(xt.m,"steps.xt.m.txt"); 
      vdump(steps.b,"steps.b.txt"); 
      vdump(steps.ax,"steps.ax.txt"); 
      vdump(steps.sx,"steps.sx.txt"); 
      vdump(steps.qx,"steps.qx.txt"); 
      vdump(steps.cn,"steps.cn.txt"); 
      */
      //-----------------------------------------------------------------------
      // Deletions
      //-----------------------------------------------------------------------
      int Ndel = findDel(data.set[iset].contig[name], pars);
      del.typeName="deletions";
      del.contigName=name;
      del.setName=setName;
      cout << " CNV find " << Ndel << " deletions" << endl;
      string fname = basename+".del.u.span.txt";
      cout << "write deletions : " << fname << endl;
      del.print(fname);      
      //-----------------------------------------------------------------------
      // Duplications
      //-----------------------------------------------------------------------
      int Ndup = findDup(data.set[iset].contig[name], pars);
      dup.typeName="duplications";
      dup.contigName=name;
      dup.setName=setName;
      cout << " CNV find " << Ndup << " duplications" << endl;
      fname = basename+".dup.u.span.txt";
      cout << "write duplications : " << fname << endl;
      dup.print(fname); 
      //-----------------------------------------------------------------------
      // Done - print summary 
      //-----------------------------------------------------------------------
      fname = basename+".summary.u.span.txt";
      cout << " CNV write summary : " << fname << endl;                    
      printSummary(fname, data.set[iset].contig[name]);
    }
  }
}

int C_DepthCNV::findDel(C_contig  & contig, RunControlParameters & pars) {
  vector<C_CNV1> d1;
  del.evt.clear();
  Ndel=0;
  int na = steps.ax.size();
  for (int ia = 0; ia<na; ia++) {
     if (steps.cn[ia]<2) {  
        d1 = getinfo(contig,pars,ia);
        int nd=d1.size();
        for (int i=0; i<nd; i++) {
          del.evt.push_back(d1[i]);
          Ndel++;
        }
     }
  }
  return Ndel;
}
int C_DepthCNV::findDup(C_contig  & contig, RunControlParameters & pars) {
  Ndup =0;
  dup.evt.clear();
  vector<C_CNV1> d1;
  int na = steps.ax.size();
  for (int ia = 0; ia<na; ia++) {
     if (steps.cn[ia]>2) {  
        d1 = getinfo(contig,pars,ia);
        int nd=d1.size();
        for (int i=0; i<nd; i++) {
          dup.evt.push_back(d1[i]);
          Ndup++;
        }
        //cout << d1;
    }
  }
  return Ndup;
}

vector<C_CNV1> C_DepthCNV::getinfo(C_contig  & contig, RunControlParameters & pars, int ia) {
    vector<C_CNV1> cnv;
    C_CNV1 cnv1;
    //--------------------------------------------------------------------------
    // split CNV if trimmed x has gap > GAP
    //--------------------------------------------------------------------------
    
    // trim position units
    int P1=(ia>0? steps.b[ia-1]: 0);
    int P2=(ia<int(bb.size()-1)? steps.b[ia]-1: int(steps.x.size()-1));
    
    // trim position vector
    vector<int> p1(1,P1);
    vector<int> p2(1,P2);
    int dp=0;     
    for (int i = P1; i<P2; i++) {
       dp++;
       int gap=mx[i+1]-mx[i];
       if (gap>GAP) {
          if (dp>(MINBINS-1) ) {
             p1.push_back(i+1);
             p2.push_back(i);
          }
       }
    }
    //--------------------------------------------------------------------------
    // number of split CNVs
    //--------------------------------------------------------------------------
    int N1=p1.size();
    for (int i=0; i<N1; i++) {
        int q1=p1[i];
        int q2=p2[N1-1-i];
        int c1=mx[q1];
        int c2=mx[q2];
        cnv1.pos=c1*u.binsize;
        cnv1.anchor=contig.getAnchorIndex();
        cnv1.length=(c2-c1)*u.binsize;
        cnv1.copynumber=int(steps.cn[ia]);
        // informative bins in this region  
        cnv1.nbi=q2-q1;   
        cnv1.q = 0;  
        cnv1.nbin=c2-c1;
        cnv1.nr = 0;
        float r2=0;
        float e2=0;
        for (int ib=c1; ib<=c2; ib++) {
          cnv1.nr+=u.r.n[ib];
          cnv1.enr+=u.e.n[ib];
          r2+=u.r.n[ib]*u.r.n[ib];
          e2+=u.e.n[ib]*u.e.n[ib];
        }
        cnv1.anr=float(cnv1.nr)/cnv1.nbin;
        cnv1.snr=sqrt(fabs(r2/cnv1.nbin-cnv1.anr*cnv1.anr));
        cnv1.aenr=float(cnv1.enr)/cnv1.nbin;
        cnv1.senr=sqrt(fabs(e2/cnv1.nbin-cnv1.aenr*cnv1.aenr));
        cnv.push_back(cnv1);
    }
    return cnv;
}

void C_DepthCNV::printSummary(string & outfilename, C_contig & c1) {
  fstream output(outfilename.c_str(), ios::out);
  if (!output) {
      cerr << "Unable to open summary file: " << outfilename << endl;
      return;
  }
  output << del.typeName << "\t " <<del.setName << "\t "<<del.contigName << endl;
  output << del.evt.size() << "events " << endl;
  output << dup.typeName << "\t " <<dup.setName << "\t "<<dup.contigName << endl;
  output << dup.evt.size() << "events " << endl;
  output << "Set " << c1.setName << "\tContig " << c1.getContigName() << endl;
  output.close();
}

 
//------------------------------------------------------------------------------
// CNV event list container class      
//------------------------------------------------------------------------------
C_CNV::C_CNV(C_contig  & c1,  RunControlParameters & pars) {
    contigName=c1.getContigName();
    setName=c1.setName;    
}

void C_CNV::print(string & outfilename) {
  fstream output(outfilename.c_str(), ios::out);
  if (!output) {
      cerr << "Unable to open CNV event output file: " << outfilename << endl;
      return;
  }
  output << *this;
  output.close();
}

ostream &operator<<(ostream &output, C_CNV & cnv)
{
  output << cnv.typeName << "\t " << cnv.setName << "\t " << cnv.contigName << endl;
  int Nevt = cnv.evt.size();
  output << Nevt << endl;
  char b[100];
  sprintf(b,"%10s %10s %4s %4s", "pos","length","anc","q");
  string s = b;
  output << s;
  sprintf(b,"%4s %10s %10s %10s %10s ", "CN","nr","enr","nbi","nbin");
  s = b;
  output << s;
  sprintf(b,"%10s %10s %10s %10s", "anr","snr","aenr","senr");
  s = b;
  output << s << endl; 
  //
  list<C_CNV1>::iterator i;
  for(i=cnv.evt.begin(); i != cnv.evt.end(); ++i) {
      output << (*i) << endl;  
  }  
  return output;
}

ostream &operator<<(ostream &output, const C_CNV1 & e1)
{
  char b [100];
  sprintf(b,"%10d %10d %4d %4d ", e1.pos,e1.length,e1.anchor,e1.q);
  string s = b;
  output << s;
  sprintf(b,"%4d %10.0f %10.2f %10d %10d ", e1.copynumber,e1.nr,e1.enr,e1.nbi,e1.nbin);
  s = b;
  output << s;
  sprintf(b,"%10.3f %10.3f %10.3f %10.3f", e1.anr,e1.snr,e1.aenr,e1.senr);
  s = b;
  output << s;
  return output;
}


C_CNV1& C_CNV1::operator=(const C_CNV1 &rhs)
{
   this->pos = rhs.pos;
   this->anchor = rhs.anchor;
   this->length = rhs.length;
   this->q = rhs.q;
   this->copynumber = rhs.copynumber;
   this->nr = rhs.nr;
   this->enr = rhs.enr;
   this->nbi   = rhs.nbi;  
   this->nbin   = rhs.nbin;  
   this->anr   = rhs.anr;  
   this->snr   = rhs.snr;  
   this->aenr   = rhs.aenr;  
   this->senr   = rhs.senr;  
   return *this;
}

//-----------------------------------------------------------------------
// Unique read coverage in big bins
//-----------------------------------------------------------------------
C_UniqueCoverage::C_UniqueCoverage() {
  binsize=1000;
}

C_UniqueCoverage::C_UniqueCoverage(C_contig & c1,  RunControlParameters & pars) {
    // RDreference model
    string RDreferencefile=pars.getRDreferenceFile();
    a = c1.loadDepth(RDreferencefile) ;
    light_t info;
    bool depthRef=a.pos1>a.pos0;
    if (depthRef) {
      info.i64=a.light;
    } else {
      // if no depth ref file then use 1kb bins 
      info.i64=1000L;
    }
    binsize=info.i16[0];
    //abintotal = info.i16[1];  
    a.calcStats();
    // 
    //cout << a.Stats.h << endl;
    //
    // read count bin size (bp)
    r = c1.countReads(binsize);
    // expected reads based on RDreference
    int N = r.n.size();
    // if no depth ref then make expected reads flat.
    if (!depthRef) {
      a.n.resize(N);
      for (int b = 0; b<N; b++) {
        a.n[b]=1;
      }
    }    
    //-------------------------------------------------------------
    // if coverage is too low, rebin until r.nzMedian exceeds pars.RDminMedian
    //-------------------------------------------------------------
    rebin(pars);
    // fetch binzie from r readcount object
    info.i64=r.light;
    binsize=info.i16[0];
    //-------------------------------------------------------------
    // expected counts 
    // calculate simple line slope with cutoff for high r.n
    //-------------------------------------------------------------
    N = r.n.size();
    e.n.resize(N);
    float num=0, den=0;
    double high=0.95;
    double toomany = r.Stats.h.p2xTrim(high)+1; 
    for (int p=0; p<N; p++) {
      if (r.n[p]<toomany) {
        num = num+a.n[p]*r.n[p];
        den = den+a.n[p]*a.n[p];
      }
    }
    double slope=num/den;    
    for (int b = 0; b<N; b++) {
        e.n[b]=slope*a.n[b];
    }
    // stats
    e.calcStats();
    e.Stats.h.setTitle("ER count of starts/base ");

   //-----------------------------------------------------------------------
   // GC content
   //-----------------------------------------------------------------------
    gc.n.resize(N);
    gc.pos0=r.pos0;
    gc.pos1=r.pos1;
    gc.light=r.light;
    string fastafile = pars.getReferenceFastaFile();
    string contigName=c1.getContigName();
    gc.GCcontent(fastafile, contigName);
    gc.calcStats(101,-0.005,1.005);
    gc.Stats.h.setTitle("GC content average/bin  ");
    slope=gc_correctionFactor();
}
//-----------------------------------------------------------------------
// Collapse bins until the median non-zero 
// count of reads exceeds RDminMedian
//-----------------------------------------------------------------------
int C_UniqueCoverage::rebin(RunControlParameters & pars) {
  // minimum number of expected reads per bin
  int minExp=pars.getRDminMedian();
  // array of rebinning factors
  double W[10] = {2,5,10,20,50,100,200,500,1000,2000};
  // original binsize from a
  light_t info;
  info.i64=r.light;
  binsize=info.i16[0];
  double W1=1;
  if (r.nzMedian>minExp) return int(W1);
  for (int iw=0; iw<10; iw++) {
    W1=W[iw];
    if ((r.nzMedian*W1)>minExp)  break;
  }
  // get out if no rebinning needed
  if (W1<1.1) return int(W1);
  
  // rebin a - alignability input
  vector<float> n0 = a.n;
  info.i16[0]=binsize*int(W1);
  a.light=info.i64;
  int N0= int(n0.size());
  int N=  int(floor(N0/W1));
  a.n.clear();
  a.n.resize(N);
  for (int b0 = 0; b0<N0; b0++) {
    int b1=int(floor((0.01+(float(b0)/W1) )));
    a.n[b1]+=n0[b0];
  }     
  a.calcStats();

  // rebin r - read counts 
  n0 = r.n;
  r.light=info.i64;
  r.n.clear();
  r.n.resize(N);
  for (int b0 = 0; b0<N0; b0++) {
    int b1=int(floor((0.01+(float(b0)/W1) )));
    r.n[b1]+=n0[b0];
  }     
  r.calcStats();
  return int(W1);
} 

void C_UniqueCoverage::write(C_contig & c1, string & basename) {
   string filename = basename+".u.readcount.span";
   c1.writeDepth(filename, r, binsize, 0);
   filename = basename+".u.expectedcount.span";
   c1.writeDepth(filename, e, binsize, 0);    
   filename = basename+".u.gc.span";
   c1.writeDepth(filename, gc, binsize, 0);    
}

// calc gc content parameter from data
int C_UniqueCoverage::gc_correctionFactor() {

    int N = gc.n.size();
    // allocate gcf size
    gcf.n.resize(N);
    // init correlation coef
    int n1=0;
    double sxy=0,  sx1=0, sy1=0, sx2=0, sy2=0, cc=0;
    // init slope 
    double num=0, den=1e-5, slope0=0; 
    // select only 20% highest expected coverage bins
    double highexp=0.8;
    double highexpectations = e.Stats.h.p2xTrim(highexp)+1; 
    // reject 5% highest read coverage bins
    double high=0.95;
    double toomany = r.Stats.h.p2xTrim(high)+1; 
    // select bins with more than 5 counts 
    double toolow = 5.0;
    // loop over all bins
    for (int p=0; p<N; p++) {
      gcf.n[p]=(gc.n[p]-gc.Stats.h.median);
      gcf.n[p]*=gcf.n[p];
      // select best bins
      if ( (r.n[p]<toomany)&&(r.n[p]>toolow)&&(e.n[p]>highexpectations )) {
        // log2 ratio
        double x1=log2(r.n[p]/e.n[p]);
        // gc correction factor form  (gc-gc0)^2
        if (gcf.n[p]>0.02) {  // increase leverage of fit to weight gc tails
          n1++;
          sx1+=gcf.n[p];
          sx2+=gcf.n[p]*gcf.n[p];
          sy1+=x1;
          sy2+=x1*x1;
          sxy+=x1*gcf.n[p];
          num = num+x1*gcf.n[p];
          den = den+gcf.n[p]*gcf.n[p];
        }
      }
    }
    if (n1>10) { 
      double xa=sx1/n1;
      double ya=sy1/n1;
      double x2a=sx2/n1;
      double y2a=sy2/n1;
      double ys=sqrt(fabs(y2a-ya*ya));
      double xs=sqrt(fabs(x2a-xa*xa));
      cc=(sxy/n1-xa*ya)/(ys*xs);
      if (fabs(cc)>0.5) {
        // gc correction slope
        double slope1=num/den;    
        // use slope only within order of magnitude....
        slope0=10*round(slope1/10);
      }
    }
    for (int p = 0; p<N; p++) {
        gcf.n[p]=slope0*gcf.n[p];
    }
    return int(slope0);
}
