//
//  UtilityFunctions.h
//  SpannerScan
//
//  Created by Chip Stewart on 5/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef UTILTYFUNCTIONS_H
#define UTILTYFUNCTIONS_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <errno.h> 
#include <dirent.h> 
#include <sys/stat.h> 
// private
//#include "Type-Hash.h"

// list of stuff to trim at ends of  strings 
//string SPACES=" \t\r\n\"";
//#define SPACES ' \n\t\r\"'
#define SPACES " \n\t\r\""

using namespace std;

//------------------------------------------------------------------------------
// prints a vector
//------------------------------------------------------------------------------
template < class T >
void printVector(const vector< T > &, const bool, const bool, ostream &);

//------------------------------------------------------------------------------
// converts a string class string to a short
//------------------------------------------------------------------------------
short string2Short (string);

//------------------------------------------------------------------------------
// converts a string class string to an int
//------------------------------------------------------------------------------
int string2Int (string);

//------------------------------------------------------------------------------
// converts a string class string to a long
//------------------------------------------------------------------------------
long string2Long (string);

//------------------------------------------------------------------------------
// converts a string class string to a long long
//------------------------------------------------------------------------------
long long string2LongLong (string);

//------------------------------------------------------------------------------
// converts a string class string to a double
//------------------------------------------------------------------------------
double string2Double (string);

//------------------------------------------------------------------------------
// converts a string class string to a float
//------------------------------------------------------------------------------
float string2Float (string);

//------------------------------------------------------------------------------
// converts a string class string to a long double
//------------------------------------------------------------------------------
long double string2LongDouble (string);


//------------------------------------------------------------------------------
// returns keys of a hash in order of associated value
//------------------------------------------------------------------------------
template< typename keyType, typename valueType >
vector<keyType> sortKeysByValue(map<keyType, valueType, less<keyType> >, bool);

//=========================================================================
// splits a string class string to a vector of strings using char delimeter
//=========================================================================
int split(vector<string>& , const string& , const string& );

//=========================================================================
// looks at a path/file stub for a list of files with a given pattern
//=========================================================================
int selectfiles(vector<string>& , const string & );

//=========================================================================
// extract path from full file stub 
//=========================================================================
string extractpath(const string & );

//=========================================================================
// extract file name from file stub 
//=========================================================================
string extractfilename(const string & );

//=========================================================================
// extract extension from file stub 
//=========================================================================
string extractfileext(const string & );

//=========================================================================
// convert int to binary number string 
//=========================================================================
string int2binary(const int);

//------------------------------------------------------------------------------
// file check
//------------------------------------------------------------------------------
//bool FileExists(string strFilename);

//------------------------------------------------------------------------------
// string to lowercase 
//------------------------------------------------------------------------------
string lowerCase(const string s);
//------------------------------------------------------------------------------
// string to Titlecase 
//------------------------------------------------------------------------------
string upperCase(const string s);
//------------------------------------------------------------------------------
// string to Titlecase 
//------------------------------------------------------------------------------
string Titlecase(const string s);
//------------------------------------------------------------------------------
// string trim utility functions
//------------------------------------------------------------------------------
string trim_right(const string & , string t = SPACES );
string trim_left(const string & s, string t = SPACES );
string trim(const string & s, string t = SPACES );

//=========================================================================
// extract text following Bam tag 
//=========================================================================
string extractBamTag(const string & line, const string & tag);
//=========================================================================
// extract text following Bam subtag 
//=========================================================================
string extractBamSubTag(const string & line, const string & tag);

#endif
