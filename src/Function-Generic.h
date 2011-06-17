//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Function-Generic
// Generic methods
// Copyright 2006 Gabor T. Marth, Boston College
// All rights reserved
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef FUNCTION_GENERIC_H
#define FUNCTION_GENERIC_H

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
#include "Type-Hash.h"

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
// string to lower 
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

#endif
