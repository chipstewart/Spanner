/*
 *  Histo.h
 *  Spanner
 *
 *  Created by Chip Stewart on 11/3/07.
 *  Copyright 2007 Boston College. All rights reserved.
 *
 */
#ifndef HISTO_H
#define HISTO_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>

using namespace std;


// structure for histogram 
class HistObj {
  friend ostream &operator<<(ostream &, const HistObj &);
  public:
    HistObj(const vector<double> & , const int , const double , const double); // constructor (do it all)
    HistObj(); // empty constructor 
    HistObj(string &); // load histogram from file 
    ~HistObj();
    HistObj&operator=(const HistObj &rhs);
    void Initialize(const int , const double , const double);  // constructor with bins only 
    void Fill(const vector<double> & , const int , const double , const double);  
    void Fill(const vector<int> & , const int , const double , const double);  
    void Fill(const vector<short> & , const int , const double , const double);
    void Fill1(const double);  
    void Fill1(const int);  
    void Fill1(const short);  
    void setTitle(const string &);
    void setXlabel(const string &);
    void setBinLabels(vector<string> &);
    void  Finalize();  
    double x2p(double); 
    double p2x(double); 
    double p2xTrim(double); 
    double x2pTrim(double); 
  	HistObj collapse( int);
	  HistObj expand();
    int Nbin;
    double xlow;
    double xhigh;
    double dx;
    double Ntot;
    double Nin;
    double Nover;
    double Nunder;
    double mode;
    double median;
    vector<double> n;
    vector<double> xc;
    vector<double> c;
    double mean;
    double std;
    double sumx;
    double sumxx;
    string title;
    string xlabel;
    vector<string> binlabels;    
    bool normalize;
    double mode1;
		bool collapsed;
		bool expanded;
};

class C_Histos {
	friend ostream &operator<<(ostream &, const C_Histos &);
public:
    C_Histos(string &); // constructor 
    C_Histos(ifstream &); // constructor 
    map<string, HistObj, less<string> > h;
    vector<string> ReadGroupTag;	
};


class C_HistoGroups {
	friend ostream &operator<<(ostream &, const C_HistoGroups &);
public:
	C_HistoGroups(); // constructor 
	C_HistoGroups(string &);
  C_Histos getHistos(string &);
 	vector<C_Histos> Groups;	
  map<string,int,less<string> > ReadGroupIndex; 
};



class StatObj {
  friend ostream &operator<<(ostream &, const StatObj &);
  public:
    StatObj(const vector<double> & , const int , const double , const double ); // constructor 
    StatObj();      // constructor 
    StatObj(string &); // load stats from file 
    void Fill(const vector<double> & , const int , const double , const double );  
    void Fill(const vector<int> & , const int , const double , const double );  
    void Fill(const vector<short> & , const int , const double , const double );  
    void Initialize(const int , const double , const double );  
    void Fill1(const double);  
    void Fill1(const int);  
    void Fill1(const short);  
    void Finalize();  
    int N;
    double mean;
    double std;
    double sumx;
    double sumxx;
    double threshold(const double);
    HistObj h;
};

#endif


