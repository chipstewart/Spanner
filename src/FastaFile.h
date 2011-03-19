/*
 *  FastaFile.h
 *  cnv
 *
 *  Created by Chip Stewart on 11/5/07.
 *  Copyright 2007 Boston College. All rights reserved.
 *
 */
#ifndef FASTAFILE_H
#define FASTAFILE_H

// standard includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
// "boost" regular expression library
// #include <boost/regex.hpp>

//"RE2" regular expression library  (replaces boost...)
#include <re2/re2.h>
using namespace re2;

// "hash_map" true hashes
#include <ext/hash_map>

using namespace std;

// type vector of for map of sequences for this ace file
// map of reads indexed by readName
typedef std::map<string, string, std::less<string> > FastaMap;

// Contig class 
class FastaObj
{
public:
   FastaObj( const string &, const string &); // constructor with specific header string
   string getSeq(const string &);
   string getSubSeq(const string &, size_t p, size_t n);
   size_t getSeqLength(const string &);
   int getNumberSeq();                                 //  number of sequences
   string getFastaFile();
   vector<string> seqNames;
   void addSeq(const string &,const string &);
   vector<size_t> getLocationsWithChar(char);
//private:
   FastaMap seq;                                  //  sequence maps
   int numberSeq;                                 //  number of sequences
   string FastaFile;
}; // end class 

#endif

