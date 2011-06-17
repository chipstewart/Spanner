/*
 *  cluster.cpp
 *  SpanDet
 *
 *  Created by Chip Stewart on 8/15/08.
 *  Copyright 2008 Boston College. All rights reserved.
 *
 */

#include "cluster.h"

// default constructor for single cluster element
C_cluster1d_element1::C_cluster1d_element1() {  
    N=0;  // number of hits in this element
    mean=0; 
    std=0;  
    median=0; 
    low=0; 
    high=0; 
}     

C_cluster1d_element1& C_cluster1d_element1::operator=(const C_cluster1d_element1 &rhs)
{
    this->N = rhs.N;
    this->mean = rhs.mean;
    this->std = rhs.std;
    this->median = rhs.median;
    this->low = rhs.low;
    this->high = rhs.high;
   this->inp.assign(rhs.inp.begin(), rhs.inp.end() );  
   return *this;
}

int C_cluster1d_element1::operator==(const C_cluster1d_element1 &rhs) const
{
    if( this->N != rhs.N) return 0;
    if( this->mean != rhs.mean) return 0;
    if( this->std != rhs.std) return 0;
    if( this->median != rhs.median) return 0;
    if( this->low != rhs.low) return 0;
    if( this->high != rhs.high) return 0;
    if( this->inp != rhs.inp) return 0;
    return 1;
}

int C_cluster1d_element1::operator<(const C_cluster1d_element1 &rhs) const
{
   // sort on low pos
   double x0 = this->low;
   double x1 = rhs.low;
   if( x0 < x1 ) return 1;
   if( x0 > x1 ) return 0;
   // then sort on low lm 
   // then sort on high pos 
   x0 = this->high;
   x1 = rhs.high;
   if( x0 < x1 ) return 1;
   if( x0 > x1 ) return 0;
   // then sort on mean pos
   x0 = this->mean;
   x1 = rhs.mean;
   if( x0 < x1 ) return 1;
   if( x0 > x1 ) return 0;
   return 0;
}

C_cluster1d_element1::~C_cluster1d_element1() {}
 

// default constructor
C_NNcluster1d::C_NNcluster1d() {
  Nmin=1;
  Smin=0;
}     

// default destructor
C_NNcluster1d::~C_NNcluster1d() {
}

// full constructor
C_NNcluster1d::C_NNcluster1d(const vector<double> & x1, const double DX 
        , const string & tn1, const string & sn1, const string & cn1) {
  typeName = tn1;
  setName = sn1;
  contigName = cn1;
  NC=0;
  cluster.clear();
  dx=DX; // neighborhood scale
  x=x1;
  N=x.size();
  init();
  Nmin=1;
  Smin=0;
  makeConnections();
  makeClusters();
}

// constructor w/o label
C_NNcluster1d::C_NNcluster1d(const vector<double> & x1, const double DX) {
  NC=0;
  cluster.clear();
  dx=DX; // neighborhood scale
  x=x1;
  N=x.size(); 
  init();
  Nmin=1;
  Smin=0;
  makeConnections();
  makeClusters();
}
    

void C_NNcluster1d::init() {
  cls.clear();
  for (int i=0; i<N; i++) {
    cls.push_back(i);
  }
  nxt=cls;
  // initialize neighbor counter y
  y.resize(N,0);
  int j0=0;
  for (int i=0; i<N; i++) {
    int j=j0;
    // find lowest x[j] within dx of x[i]
    while ((j<N)&((x[i]-x[j])>dx)) { j++; }
    j0 = j;
    // find lowest x[j] within dx of x[i]
    while (fabs(x[i]-x[j])<=dx) { 
      y[i]+=1;
      j++;
      if (j==N) break; 
    }    
  }

  
}


void C_NNcluster1d::makeConnections() {
  int j0=0; 
  for (int i=0; i<N; i++) {
    int j=j0;
    // find lowest x[j] within dx of x[i]
    while ((j<N)&((x[i]-x[j])>dx)) { j++; }
    j0 = j;
    // find lowest x[j] within dx of x[i]
    double ymax = 0;
    int ij = -1;
    while (fabs(x[i]-x[j])<=dx) { 
      if (y[j]>ymax) {
        ij = j;
        ymax = y[j];
      }
      j++; 
      if (j==N) break; 
    }  
    connect(i,ij);  
  }
}

void C_NNcluster1d::connect(int i1,int i2) {
  if (i1==i2) return; 
  int nxt1 = nxt[i1];
  //int nxt2 = nxt[i2];
  int cls1 = cls[i1];
  nxt[i1]=i2;
  int i = i2;
  while (nxt[i]!=i2) {
    cls[i]=cls1;
    i = nxt[i];
  }
  nxt[i]=nxt1;
  cls[i]=cls1;      
}

void C_NNcluster1d::makeClusters() {
  for (int i=0; i<N; i++) {
    if (cluster.count(cls[i])==0) {
       cluster[cls[i]].N=1;
       cluster[cls[i]].mean=x[i];
       cluster[cls[i]].std=x[i]*x[i];
       cluster[cls[i]].low=x[i];
       cluster[cls[i]].high=x[i];
       cluster[cls[i]].inp.clear();
       cluster[cls[i]].inp.push_back(i);
    } else {
       cluster[cls[i]].N+=1;
       cluster[cls[i]].mean+=x[i];
       cluster[cls[i]].std+=x[i]*x[i];
       cluster[cls[i]].low=(x[i]>cluster[cls[i]].low?cluster[cls[i]].low: x[i]);
       cluster[cls[i]].high=(x[i]<cluster[cls[i]].high?cluster[cls[i]].high: x[i]);
       cluster[cls[i]].inp.push_back(i);
    }
  }
  C_cluster1d_elements::iterator it;
  C_cluster1d_element1 c1;
  for ( it=cluster.begin() ; it != cluster.end(); it++ ) {
      int i = (*it).first;
      cluster[i].mean=cluster[i].mean/cluster[i].N;
      cluster[i].std=sqrt(cluster[i].std/cluster[i].N-(cluster[i].mean*cluster[i].mean));
  }
  NC = cluster.size();
}

void C_NNcluster1d::cleanClusters() {
  C_cluster1d_elements::iterator it;
  C_cluster1d_element1 c1;
  vector<int> itoss;
  //int N0 = cluster.size();
  for ( it=cluster.begin() ; it != cluster.end(); it++ ) {
      int i = (*it).first;
      bool toss = cluster[i].N<=Nmin;
      toss = toss | (cluster[i].std<=Smin);
      if (toss) itoss.push_back(i);
  }
  int NT = itoss.size();
  for (int i=0; i<NT; i++) {
     cluster.erase(itoss[i]);
  }
  NC = cluster.size();
  // printf(" cleanCluster removed %d of %d clusters leaving NC %d\n",NT,N0,NC);
}        

//------------------------------------------------------------------------------
// merge near clusters
//------------------------------------------------------------------------------
void C_NNcluster1d::mergeClusters() {
  C_cluster1d_elements::iterator it,jt;
  int N0 = cluster.size();
  int NM=0;
  for ( it=cluster.begin() ; it != cluster.end(); it++ ) {
      if (it != cluster.begin()) {
        jt = it;
        jt--;  
        int i = (*it).first;
        int j = (*jt).first;
        bool morethan1 = cluster[j].N>1;
        bool thin = cluster[j].std<dx;
        bool close = (cluster[i].low-cluster[j].high)<(dx);
        close = close & (fabs(cluster[i].mean-cluster[j].mean)<(2*dx));
        if (thin&close&morethan1)  {
           int ix=cluster[i].inp[0];
           int jx=cluster[j].inp[0];
           connect(ix,jx); 
           NM++;
        }
      }
  }
  cluster.clear();
  makeClusters();
  NC = cluster.size();
  printf(" mergeCluster merged %d of %d clusters leaving NC %d\n",NM,N0,NC);
}        

        
ostream &operator<<(ostream &output,  C_cluster1d_element1 & c1)
{
    output << c1.mean << "\t " << c1.N << "\t " << c1.std << "\t " << c1.low << "\t " << c1.high << endl ;      
    return output;
}

ostream &operator<<(ostream &output,  C_cluster1d_elements & c1)
{
  C_cluster1d_elements::iterator it;
  //C_cluster1d_element1 c1;
  output << "mean" << "\t " << "N" << "\t " << "std" << "\t " << "low" << "\t " << "high" << endl ;      
  for ( it=c1.begin() ; it != c1.end(); it++ ) {
      int i = (*it).first;
      output << c1[i];
  }
  return output;
}

ostream &operator<<(ostream &output,  C_NNcluster1d & nn1)
{
 output << "Number of clusters" << "\t " << nn1.NC << endl ;
 output << nn1.cluster;
 return output;
}


void C_NNcluster1d::write(string & outfilename)
{
  // format: optimized to loadFragments.m matlab script  
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out | ios::binary);
  if (!output) {
      cerr << "Unable to open file: " << outfilename << endl;
      return;
  }
  C_headerSpan h;
  h.V = 1203;
  h.setName    = setName;
  h.contigName = contigName;
  h.typeName   = typeName;
  h.reclen =   4*sizeof(double)+sizeof(int);
  h.N = this->cluster.size();
  h.write(output);
  C_cluster1d_elements::iterator it;
  C_cluster1d_element1 c1;
  for ( it=cluster.begin() ; it != cluster.end(); it++ ) {
    int i = (*it).first;
    c1 = cluster[i];
    double mean1 = c1.mean;  
    int N1 = c1.N;  
    double std1 = c1.std;  
    double low1 = c1.low;  
    double high1 = c1.high;  
    //position
    output.write(reinterpret_cast<const char *>(&mean1), sizeof(double));
    //length
    output.write(reinterpret_cast<const char *>(&N1), sizeof(int));
    output.write(reinterpret_cast<const char *>(&std1), sizeof(double));
    output.write(reinterpret_cast<const char *>(&low1), sizeof(double));
    output.write(reinterpret_cast<const char *>(&high1), sizeof(double));
  }
  output.close();
}


// default constructor for single cluster element
C_cluster2d_element1::C_cluster2d_element1() {  
    N=0;  // number of hits in this element
    for (int j=0; j<2; j++) {
      mean[j]=0;
      std[j]=0;
      median[j]=0; 
      low[j]=0; 
      high[j]=0;
    }
}     

C_cluster2d_element1::~C_cluster2d_element1() {}
 
C_cluster2d_element1& C_cluster2d_element1::operator=(const C_cluster2d_element1 &rhs)
{
   this->N = rhs.N;
   for (int i=0; i<2; i++) {
     this->mean[i] = rhs.mean[i];
     this->std[i] = rhs.std[i];
     this->median[i] = rhs.median[i];
     this->low[i] = rhs.low[i];
     this->high[i] = rhs.high[i];
   }
   this->inp.assign(rhs.inp.begin(), rhs.inp.end() );  
//   this->inp =  rhs.inp;
   return *this;
}

int C_cluster2d_element1::operator==(const C_cluster2d_element1 &rhs) const
{
   if( this->N != rhs.N) return 0;
   for (int i=0; i<2; i++) {
     if( this->mean[i] != rhs.mean[i]) return 0;
     if( this->std[i] != rhs.std[i]) return 0;
     if( this->median[i] != rhs.median[i]) return 0;
     if( this->low[i] != rhs.low[i]) return 0;
     if( this->high[i] != rhs.high[i]) return 0;
   }
   if( this->inp != rhs.inp) return 0;
   return 1;
}

int C_cluster2d_element1::operator<(const C_cluster2d_element1 &rhs) const
{
   // sort on low pos
   double x0 = this->low[0];
   double x1 = rhs.low[0];
   if( x0 < x1 ) return 1;
   if( x0 > x1 ) return 0;
   // then sort on low lm 
   x0 = this->low[1];
   x1 = rhs.low[1];
   if( x0 < x1 ) return 1;
   if( x0 > x1 ) return 0;
   // then sort on high pos 
   x0 = this->high[0];
   x1 = rhs.high[0];
   if( x0 < x1 ) return 1;
   if( x0 > x1 ) return 0;
   // then sort on high lm
   x0 = this->high[1];
   x1 = rhs.high[1];
   if( x0 < x1 ) return 1;
   if( x0 > x1 ) return 0;
   // then sort on mean pos
   x0 = this->mean[0];
   x1 = rhs.mean[0];
   if( x0 < x1 ) return 1;
   if( x0 > x1 ) return 0;
   // then sort on mean lm
   x0 = this->mean[1];
   x1 = rhs.mean[1];
   if( x0 < x1 ) return 1;
   if( x0 > x1 ) return 0;
   return 0;
}


// default constructor
C_NNcluster2d::C_NNcluster2d() {
  Nmin=1;
  Smin[0]=0;
  Smin[1]=0;
  NC=0;
  N=0;
}     

// default destructor
C_NNcluster2d::~C_NNcluster2d() {
}

// full constructor
C_NNcluster2d::C_NNcluster2d(const vector<vector<double> > & x1,const vector<vector<double> > & wx1
        , const vector<int> & ip1
        , const double fx1[2] , const string & tn1, const string & sn1, const string & cn1) {
  typeName = tn1;
  setName = sn1;
  contigName = cn1;
  Nmin=1;
  Smin[0]=-1;
  Smin[1]=-1;
  NC=0;
  cluster.clear();
  fx[0]=fx1[0]; // neighborhood scale
  fx[1]=fx1[1]; 
  x=x1;
  wx=wx1;
  ip=ip1;
  N=x.size();
  dx[0]=0;
  dx[1]=0;
  for (int i=0; i<N; i++) {
      wx[i][0]*=fx[0];
      wx[i][1]*=fx[1];
      if ( wx[i][0]>dx[0] ) { dx[0]=wx[i][0]; }  
      if ( wx[i][1]>dx[1] ) { dx[1]=wx[i][1]; }  
  }
  init();
  makeConnections();
  makeClusters();
}
    

void C_NNcluster2d::init() {
  cls.clear();
  for (int i=0; i<N; i++) {
    cls.push_back(i);
  }
  nxt=cls;
  // initialize neighbor counter y
  n.resize(N,0);
  int j0=0;
  for (int i=0; i<N; i++) {
    int j=j0;
    // find lowest x[j] within dx of x[i]
    while ((j<N)&((x[i][0]-x[j][0])>wx[i][0])) { j++; }
    j0 = j;
    // find highest x[j] within dx of x[i]
    while (fabs(x[i][0]-x[j][0])<=wx[i][0]) { 
      if (fabs(x[i][1]-x[j][1])<=wx[i][1]) n[i]+=1;
      j++; 
      if (j==N) break;
    }    
  }  
}


void C_NNcluster2d::makeConnections() {
  int j0=0; 
  for (int i=0; i<N; i++) {
    int j=j0;
    // find lowest x[j] within dx of x[i]
    //while ((x[i][0]-x[j][0])>dx[0]) { j++; }
    while ((j<N)&((x[i][0]-x[j][0])>wx[i][0])) { j++; }
    j0 = j;
    // find lowest x[j] within dx of x[i]
    int nmax = 0;
    int ij = -1;
    while (fabs(x[i][0]-x[j][0])<=wx[i][0]) { 
      if (fabs(x[i][1]-x[j][1])<=wx[i][1]) { 
        if (n[j]>nmax) {
          ij = j;
          nmax = n[j];
        }
      }
      j++; 
      if (j==N) break;
    }  
    connect(i,ij);  
  }
}
 
void C_NNcluster2d::connect(int i1,int i2) {
  if (i1==i2) return; 
  int nxt1 = nxt[i1];
  //int nxt2 = nxt[i2];
  int cls1 = cls[i1];
  nxt[i1]=i2;
  int i = i2;
  while (nxt[i]!=i2) {
    cls[i]=cls1;
    i = nxt[i];
  }
  nxt[i]=nxt1;
  cls[i]=cls1;      
}

void C_NNcluster2d::makeClusters() {
  for (int i=0; i<N; i++) {
    if (cluster.count(cls[i])==0) {
      cluster[cls[i]].N=1;
      for (int j=0; j<2; j++) {
        cluster[cls[i]].mean[j]=x[i][j];
        cluster[cls[i]].std[j]=x[i][j]*x[i][j];
        cluster[cls[i]].low[j]=x[i][j];
        cluster[cls[i]].high[j]=x[i][j];
      }
      cluster[cls[i]].inp.clear();
      cluster[cls[i]].inp.push_back(ip[i]);
    } else {
      cluster[cls[i]].N+=1;
      for (int j=0; j<2; j++) {
        cluster[cls[i]].mean[j]+=x[i][j];
        cluster[cls[i]].std[j]+=x[i][j]*x[i][j];
        cluster[cls[i]].low[j]=(x[i][j]>cluster[cls[i]].low[j]?cluster[cls[i]].low[j]: x[i][j]);
        cluster[cls[i]].high[j]=(x[i][j]<cluster[cls[i]].high[j]?cluster[cls[i]].high[j]: x[i][j]);
      }
      cluster[cls[i]].inp.push_back(ip[i]);
    }
  }
  C_cluster2d_elements::iterator it;
  C_cluster2d_element1 c1;
  for ( it=cluster.begin() ; it != cluster.end(); it++ ) {
      int i = (*it).first;
      for (int j=0; j<2; j++) {
        cluster[i].mean[j]=cluster[i].mean[j]/cluster[i].N;
        cluster[i].std[j]=sqrt(cluster[i].std[j]/cluster[i].N-(cluster[i].mean[j]*cluster[i].mean[j]));
      }
  }
  NC = cluster.size();
}

void C_NNcluster2d::cleanClusters() {
  C_cluster2d_elements::iterator it;
  C_cluster2d_element1 c1;
  vector<int> itoss;
  //int N0 = cluster.size();
  for ( it=cluster.begin() ; it != cluster.end(); it++ ) {
      int i = (*it).first;
      bool toss = cluster[i].N<=Nmin;
      for (int j=0; j<2; j++) {
         toss = toss | (cluster[i].std[j]<=Smin[j]);
      }
      if (toss) itoss.push_back(i);
  }
  int NT = itoss.size();
  for (int i=0; i<NT; i++) {
     cluster.erase(itoss[i]);
  }
  NC = cluster.size();
  //printf(" cleanCluster removed %d of %d clusters leaving NC %d\n",NT,N0,NC);
}        

void C_NNcluster2d::mergeClusters() {
  C_cluster2d_elements::iterator it,jt;
  int N0 = cluster.size();
  int NM=0;
  for ( it=cluster.begin() ; it != cluster.end(); it++ ) {
      if (it != cluster.begin()) {
        jt = it;
        jt--;  
        int i = (*it).first;
        int j = (*jt).first;
        bool morethan1 = cluster[j].N>1;
        bool thin = cluster[j].std[0]<dx[0];
        bool close = (cluster[i].low[0]-cluster[j].high[0])<2*dx[0];
        close = close & (fabs(cluster[i].mean[1]-cluster[j].mean[1])<dx[0]);
        if (thin&close&morethan1)  {
          //  int ix=cluster[i].inp[0];  why inp? 
          //  int jx=cluster[j].inp[0];
          //  connect(ix,jx); 
           connect(i,j); 
           NM++;

           cerr << "merge " << NM << endl;

        }
      }
  }
  cluster.clear();
  makeClusters();
  NC = cluster.size();
  printf(" mergeCluster merged %d of %d clusters leaving NC %d\n",NM,N0,NC);
}        


ostream &operator<<(ostream &output,  C_cluster2d_element1 & c1)
{
    output << c1.N << "\t ";
    for (int j=0; j<2; j++) {
      output << c1.mean[j] << "\t " << c1.std[j] << "\t " << c1.low[j] << "\t " << c1.high[j]  ;      
    }
    output << endl;
    return output;
}

ostream &operator<<(ostream &output,  C_cluster2d_elements & c1)
{
  C_cluster2d_elements::iterator it;
  //C_cluster1d_element1 c1;
  output << "N" << "\t " << "mean1" << "\t " << "std1" << "\t " << "low1" << "\t " << "high1";
  output << "\t " << "mean2" << "\t " << "std2" << "\t " << "low2" << "\t " << "high2" << endl ;      
  for ( it=c1.begin() ; it != c1.end(); it++ ) {
      int i = (*it).first;
      output << c1[i];
  }
  return output;
}

ostream &operator<<(ostream &output,  C_NNcluster2d & nn2)
{
 output << "Number of clusters" << "\t " << nn2.NC << endl ;
 output << nn2.cluster;
 return output;
}

void C_NNcluster2d::setSmin(double Smin1, double Smin2) {
  Smin[0]=Smin1;
  Smin[1]=Smin2;
}

void C_NNcluster2d::setNmin(int Nmin1) {
  Nmin=Nmin1;
}

void C_NNcluster2d::write(string & outfilename)
{
  // format: optimized to loadCluster2Span.m matlab script  
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out | ios::binary);
  if (!output) {
      cerr << "Unable to open file: " << outfilename << endl;
      return;
  }
  C_headerSpan h;
  h.V=1203;
  h.setName    = setName;
  h.contigName = contigName;
  h.typeName   = typeName;
  h.reclen =   8*sizeof(double)+sizeof(int);
  h.N = this->cluster.size();
  h.write(output);
  C_cluster2d_elements::iterator it;
  C_cluster2d_element1 c1;
  int Npair = 0;
  for ( it=cluster.begin() ; it != cluster.end(); it++ ) {
    int i = (*it).first;
    c1 = cluster[i];
    output.write(reinterpret_cast<const char *>(&c1.N), sizeof(int));
    Npair+=c1.N;
    for (int j=0; j<2; j++) {
      output.write(reinterpret_cast<const char *>(&c1.mean[j]), sizeof(double));
      output.write(reinterpret_cast<const char *>(&c1.std[j]), sizeof(double));
      output.write(reinterpret_cast<const char *>(&c1.low[j]), sizeof(double));
      output.write(reinterpret_cast<const char *>(&c1.high[j]), sizeof(double));
    }
  }
  // corresponding pairs data
  h.setName    = setName;
  h.contigName = contigName;
  h.typeName   = "localpairs";
  h.reclen =   sizeof(int)+2*sizeof(int)+2*sizeof(short)+2*sizeof(char);
  h.N = Npair;
  h.write(output);
  
  
  output.close();
}

/*
  list<C_localpair>::iterator i;
  for(i=localpairs.begin(); i != localpairs.end(); ++i) {
    unsigned int pos = (*i).pos;  
    int lm = (*i).lm;  
    char o = (*i).orient;  
    char q = (*i).q;  
    short len1 = (*i).len1;  
    short len2 = (*i).len2;  
    //position
    output.write(reinterpret_cast<const char *>(&pos), sizeof(int));
    //length
    output.write(reinterpret_cast<const char *>(&lm), sizeof(int));
    output.write(reinterpret_cast<const char *>(&o), sizeof(char));
    output.write(reinterpret_cast<const char *>(&len1), sizeof(short));
    output.write(reinterpret_cast<const char *>(&len2), sizeof(short));
    output.write(reinterpret_cast<const char *>(&q), sizeof(char));    
  }
  */

