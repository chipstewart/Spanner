/*
 *  headerInfo.h
 *  SpanDet
 *
 *  Created by Chip Stewart on 8/19/08.
 *  Copyright 2008 Boston College. All rights reserved.
 *
 */
#ifndef HEADERSPAN_H
#define HEADERSPAN_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <list>
//#include <boost/regex.hpp>
// private
#include "Type-Hash.h"
#include "RunControlParameterFile.h"
#include "Function-Generic.h"

using namespace std;

//------------------------------------------------------------------------------
// header  info container class
//------------------------------------------------------------------------------
class C_headerSpan {
  friend ostream &operator<<(ostream &, C_headerSpan &);
  public:
    C_headerSpan();                                 // constructor
    C_headerSpan(fstream &); 
    ~C_headerSpan(){};                              // destructor
    bool write(fstream &); 
    int V; 
    string contigName;                              // contig/anchor name
    string setName;                                 // set name
    string typeName;                                // info type name
    long long light;                                // 8 bytes of info (type specific)
    unsigned int reclen;
    unsigned int N;
    //--------------------------------------------------------------------------
    // Spanner file types: 
    // 0=local pair, 1=cross pair, 2=singleton, 3=dstarts, 4=dEnds, 5=uEnds, 6=rpt
    // clusters: 10=local pair, 11=cross pair, 12=singleton, 13=dstarts, 14=dEnds, 15=uEnds
    //--------------------------------------------------------------------------
    unsigned int type;
    // type file extensions for paired-end analysis
    vector<string> spanext; 
    // type file extensions for clusters
    vector<string> spanextc; 
    // type file extensions for single-end read analysis
    // vector<string> spanext1; 
};

#endif

