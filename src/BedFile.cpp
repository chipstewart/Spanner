/*
 *  BedFile.cpp
 *  Spanner
 *
 *  Created by Chip Stewart on 1/28/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *
 */

#include "BedFile.h"

bool comparePair (pair<int,short> p1, pair<int, short> p2)
{
  if (p1.first<p2.first) {return true;}
  if (p1.second<p2.second) {return true;}
  return false;
}


//------------------------------------------------------------------------------
// Bed class constructor 
//------------------------------------------------------------------------------
C_BedChr::C_BedChr() {
    N = 0;
    plim[0] =INT_MAX;
    plim[1] =0;
    name= "";
    header= "";
}



C_BedChr::C_BedChr(string & fn, string &chr1)  { // load  bed record for chr    
	
  bool sort=false;
  
  this->filename=fn;
  if (fn.length()==0) {
     return;
  }

	this->chr=chr1;
  if (chr1.length()==0) {
     return;
  }
  
  // input line
  string line;

  // open
  ifstream bed(this->filename.c_str(), ios::in);
  
  if (!bed) {
    cerr << "Unable to open file: " << this->filename<< endl;
    exit(1);
  }
  
  this->header = "";
  
  int nl=0;
  while (getline(bed, line)) {

    nl++;
    
    if (line.compare(0,3,"chr") != 0) {
       this->header=this->header+"\n"+line;
    } else {
    
      // construct a stream from the string
      stringstream strstr(line);
      // use stream iterators to copy the stream to the vector as whitespace separated strings
      istream_iterator<string> it(strstr);
      istream_iterator<string> end;
      vector<string> f(it, end);
      
      if (f.size()<3) {
          cerr << fn << " too few fields:\n " << line << endl;
          continue;
      }
  
      string c=f[0];
      if (c.compare(3,c.size()-3,chr1)==0) {
        int p1=  atoi(f[1].c_str());
        int p2=  atoi(f[2].c_str());
        int l1 = p2-p1;
        if (p.size()>0) {
          if (p1<p[p.size()-1]) { sort=true;};
        }
        this->p.push_back(p1);
        this->l.push_back(l1);
        if (p1<plim[0]) { plim[0]=p1;};
        if (p2>plim[1]) { plim[1]=p2;};
      }
    }
  }
  
  bed.close();
 
  if (sort) {
     
     pair <int,short>  q1;
     list < pair <int,short> > q;
     list < pair <int,short> >::iterator iq;
     for (int i=0; i<int(p.size()); i++) {
       q1.first=p[i];
       q1.second=l[i];
       q.push_back(q1);
     }
  
     q.sort(comparePair);
     
     int i=0;
     for (iq=q.begin(); iq!=q.end() ; iq++) {
       p[i]=(*iq).first;
       l[i]=(*iq).second;
       i++;
     }
  }
  N=p.size();     
}


C_BedChr& C_BedChr::operator=(const C_BedChr &rhs) {
    N=rhs.N ;
    chr=rhs.chr;
    name=rhs.name;
    header=rhs.header;
    filename=rhs.header;
    plim[0]=rhs.plim[0];
    plim[1]=rhs.plim[1];
    p=rhs.p ;
    l=rhs.l ;
    return *this;
}


