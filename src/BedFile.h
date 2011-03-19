/*
 *  BedFile.h
 *  Spanner
 *
 *  Created by Chip Stewart on 1/28/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *
 */

#ifndef BEDFILE_H
#define BEDFILE_H

#include <iostream>
#include <ostream>
#include <istream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <limits.h>


using namespace std;

// structure for Bed File
class C_BedChr {
  friend ostream &operator<<(ostream &, const C_BedChr &);
  public:
    C_BedChr(); // empty constructor
    C_BedChr(string &); // create bed record for chr
    C_BedChr(string &, string &); // load  bed record for chr
   // ~C_BedChr();
    C_BedChr&operator=(const C_BedChr &rhs);
    C_BedChr SelectRegion(int, int);  // return records within chr window p1 p2
    void setChr(const string &);
    void setHeader(const string &);
    int N;
    string chr;
    string name;
    string header;
    string filename;
    int plim[2];
    vector<int> p;
    vector<short> l;
};    

/*
// structure for Bed File
class C_BedFile {
  friend ostream &operator<<(ostream &, const C_BedFile &);
  public:
    C_BedFile(); // empty constructor
    C_BedFile(string &); // load bed info from file 
    C_BedFile(string &, string &); // load info from file within chr 
    ~C_BedFile();
    C_BedFile&operator=(const C_BedFile &rhs);
    C_BedChr SelectChr(string &);  // return records within chr
    C_BedChr SelectRegion(string &, int, int);  // return records within chr window p1 p2
    void setNames(const string &);
    void setHeader(const string &);
    int Nchr;
    map<string, C_BedChr, less<string> > c;
    string name;
    string header;
};
*/

#endif


 

