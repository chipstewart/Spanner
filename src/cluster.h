/*
 *  cluster.h
 *  SpanDet
 *
 *  Created by Chip Stewart on 8/15/08.
 *  Copyright 2008 Boston College. All rights reserved.
 *
 */
#ifndef CLUSTER_H
#define CLUSTER_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <map>
#include <list>
#include <utility>
#include <math.h>
#include "headerSpan.h"

using namespace std;

// structure for 1d NN cluster object 
class C_cluster1d_element1 {
  friend ostream &operator<<(ostream &, const C_cluster1d_element1 &);
  public:
    C_cluster1d_element1(); // constructor 
    ~C_cluster1d_element1();
    C_cluster1d_element1&operator=(const C_cluster1d_element1 &rhs);
    int operator==(const C_cluster1d_element1 &rhs) const;
    int operator<(const C_cluster1d_element1 &rhs) const;
    // data
    int N;  // number of hits in this element
    double mean; 
    double std;  
    double median; 
    double low; 
    double high; 
    //list<int> ind;    
    vector<int> inp;    
};

class C_cluster2d_element1 {
  friend ostream &operator<<(ostream &, const C_cluster2d_element1 &);
  public:
    C_cluster2d_element1(); // constructor 
    ~C_cluster2d_element1();
    C_cluster2d_element1&operator=(const C_cluster2d_element1 &rhs);
    int operator==(const C_cluster2d_element1 &rhs) const;
    int operator<(const C_cluster2d_element1 &rhs) const;
    // data
    int N;  // number of hits in this element
    double mean[2]; 
    double std[2];  
    double median[2]; 
    double low[2]; 
    double high[2]; 
    vector<int> inp;    
};


typedef std::map<int, C_cluster1d_element1, std::less<int> >  C_cluster1d_elements;
typedef std::map<int, C_cluster2d_element1, std::less<int> >  C_cluster2d_elements;

// structure for 1d NN cluster object 
class C_NNcluster1d {
  friend ostream &operator<<(ostream &, const C_NNcluster1d &);
  public:
    C_NNcluster1d(const vector<double> & , const double ); // constructor (input list and NN scale)
    C_NNcluster1d(const vector<double> & , const double,const string &, const string &,const string &); // constructor (input list and NN scale)
    C_NNcluster1d(); // empty constructor 
    ~C_NNcluster1d();
    void write(string &);
    C_cluster1d_elements  cluster;
    int NC;
    string typeName;
    string contigName;
    string setName;    
    void cleanClusters(); 
    void mergeClusters();   
  private:
    // alg steps
    void init();
    void makeConnections();
    void makeClusters();
    void connect(int,int);
    // data
    int N;              // number of elements to NN cluster
    double dx;          // neighborhood scale
    vector<double> x;   // points in 1d
    vector<double> y;   // local density of points
    vector<int> cls;    // class label
    vector<int> nxt;    // next pointer
    int Nmin;
    double Smin;
};


// structure for 1d NN cluster object 
class C_NNcluster2d {
  friend ostream &operator<<(ostream &, const C_NNcluster2d &);
  public:
    // constructor (input list and NN scale)
    C_NNcluster2d(const vector<vector<double> > & , const vector<vector<double> > & 
        , const vector<int> & 
        , const double[2], const string &, const string &,const string &); 
    C_NNcluster2d(); // empty constructor 
    ~C_NNcluster2d();
    void write(string &);
    C_cluster2d_elements  cluster;
    int NC;
    string typeName;
    string contigName;
    string setName;    
    void setNmin(int);
    void setSmin(double,double);
    void cleanClusters();    
    void mergeClusters();
  private:
    // alg steps
    void init();
    void makeConnections();
    void makeClusters();
    void connect(int,int);
    // data
    int N;                        // number of elements to NN cluster
    double dx[2];                 // max absolute neighborhood scale 
    double fx[2];                 // relative neighborhood scale 
    vector<vector<double> > x;    // data to cluster
    vector<vector<double> > wx;   // local window (for each library)
    vector<int> ip;               // index to elements
    vector<int> n;                // local density
    vector<int> cls;
    vector<int> nxt;
    int Nmin;
    double Smin[2];

};

#endif

