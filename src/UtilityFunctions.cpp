
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//
//  UtilityFunctions.cpp
//  SpannerScan
//
//  Created by Chip Stewart on 5/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
#include "UtilityFunctions.h"


#ifndef UTILITYFUNCTIONS_CPP
#define UTILITYFUNCTIONS_CPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <sys/stat.h> 


using std::ios;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::istream;
using std::cin;
using std::cout;
using std::clog;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::multimap;

//#include "Type-Hash.h"

//------------------------------------------------------------------------------
// prints a vector
//------------------------------------------------------------------------------
template < class T >
void printVector(const vector< T > &vectorRef, const bool separator, const bool bracket, ostream &output2) {
    if ( vectorRef.empty() ) {
        if (bracket) {
            output2 << "()";
        }
    }
    else {
        if (bracket) {
            output2 << "(";
        }
        if (separator) {
            std::ostream_iterator< T > output( output2, " " );
            std::copy( vectorRef.begin(), vectorRef.end(), output );
        }
        else {
            std::ostream_iterator< T > output( output2 );
            std::copy( vectorRef.begin(), vectorRef.end(), output );
        }
        if (bracket) {
            output2 << ")";
        }
    }
}

//------------------------------------------------------------------------------
// converts a string class string to a short
//------------------------------------------------------------------------------
short string2Short (string sString) {
    
    // convert it to C-style string
    const char *cString = sString.c_str();
    
    // convert it to integer
    short iString = atoi(cString);
    
    // return
    return iString;
}

//------------------------------------------------------------------------------
// converts a string class string to an int
//------------------------------------------------------------------------------
int string2Int (string sString) {
    
    // convert it to C-style string
    const char *cString = sString.c_str();
    
    // convert it to integer
    int iString = atoi(cString);
    
    // return
    return iString;
}

//------------------------------------------------------------------------------
// converts a string class string to a long int
//------------------------------------------------------------------------------
long string2Long (string sString) {
    
    // convert it to C-style string
    const char *cString = sString.c_str();
    
    // convert it to integer
    long lString = atol(cString);
    
    // return
    return lString;
}

//------------------------------------------------------------------------------
// converts a string class string to a long long
//------------------------------------------------------------------------------
long long string2LongLong (string sString) {
    
    // convert it to C-style string
    const char *cString = sString.c_str();
    
    // convert it to integer
    long long llString = atol(cString);
    
    // return
    return llString;
}

//------------------------------------------------------------------------------
// converts a string class string to a double
//------------------------------------------------------------------------------
double string2Double (string sString) {
    
    // convert it to C-style string
    const char *cString = sString.c_str();
    
    // convert it to integer
    double fString = atof(cString);
    
    // return
    return fString;
}

//------------------------------------------------------------------------------
// converts a string class string to a float
//------------------------------------------------------------------------------
float string2Float (string sString) {
    
    // convert it to C-style string
    const char *cString = sString.c_str();
    
    // convert it to integer
    float fString = float(atof(cString));
    
    // return
    return fString;
}

//------------------------------------------------------------------------------
// converts a string class string to an double
//------------------------------------------------------------------------------
long double string2LongDouble (string sString) {
    
    // convert it to C-style string
    const char *cString = sString.c_str();
    
    // convert it to integer
    long double fString = atof(cString);
    
    // return
    return fString;
}



//------------------------------------------------------------------------------
// sorts keys of a hash in order of associated value
//------------------------------------------------------------------------------
template< typename keyType, typename valueType >
vector<keyType> sortKeysByValue(map<keyType, valueType, std::less<keyType> > hash, bool descend) {
    
    // instantiate inverse has as a multimap
    multimap<valueType, keyType, std::less<valueType> > inverseHash;
    
    // load elements of hash into inverseHash
    for (typename map<keyType, valueType, std::less<keyType> >::const_iterator iter = hash.begin(); iter != hash.end(); iter++) {
        keyType key = iter->first;
        valueType value = iter->second;
        inverseHash.insert(typename multimap<valueType, keyType, std::less<valueType> >::value_type(value, key));
    }
    
    // compose vector of original keys sorted by original values
    vector<keyType> sortedKeys;
    for(typename multimap<valueType, keyType, std::less<valueType> >::const_iterator iter = inverseHash.begin(); 
        iter != inverseHash.end(); iter++) {
        keyType key = iter->second;
        sortedKeys.push_back(key);
    }
    
    // reverse if descending order was required
    if (descend) {
        reverse(sortedKeys.begin(), sortedKeys.end());
    }
    
    // return
    return sortedKeys;
}

//=============================================================
// Split a string on a given character into a vector of strings
// The vector is passed by reference and cleared each time
// The number of strings split out is returned
// William S. Lear rael at see.sig
// Mon May 10 17:07:49 CEST 1999
//=============================================================
int  split(vector<string>& v, const string& str,
           const string& delimiters = " ")
{
    v.clear();
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
    
    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        v.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
    return int(v.size());
}
//=========================================================================
// looks at a path/file stub for a list of files with a given pattern
//=========================================================================
int selectfiles(vector<string>& s, const string & fstub)  {
    DIR *pdir;
    struct dirent *pent;
    size_t pos;
    // parse in path, name, ext  
    string path =extractpath(fstub);
    string name =extractfilename(fstub);
    // open path
    pdir=opendir(path.c_str()); //"." refers to the current dir
    if (!pdir){
        cerr << "opendir() failure in selectfiles: " << fstub << endl;
        exit(1);
    }
    // loop over files
    errno=0; 
    while ((pent=readdir(pdir))){
        string file1 = pent->d_name;
        pos = file1.find(name);
        if (string::npos != pos ) {
            // cout << file1 << endl;    
            string s1 = path+"/"+file1;
            s.push_back(s1);
        }
    }
    if (errno){
        cerr << "readdir() failure in selectfiles: " << fstub << endl;
        exit(1);
    }
    closedir(pdir);
    return int(s.size());
}

//=========================================================================
// extract path from full file stub - before last / 
//=========================================================================
string extractpath(const string & filename) {
    string path = ".";
    string::size_type pos = filename.find_last_of("/");
    if (string::npos != pos ) { 
        path = filename.substr(0,pos);
    }
    return path;
}

//=========================================================================
// extract file name from file stub - between last / and last .
//=========================================================================
string extractfilename(const string & filename) {
    string name="";
    string::size_type pos = filename.find_last_of("/");
    string::size_type pos1 = filename.find_last_of(".");
    if (string::npos == pos1 ) pos1=filename.size(); 
    if (string::npos == pos ) pos=-1;
    pos=pos+1;
    name = filename.substr(pos,pos1-pos);
    return name;
}
//=========================================================================
// extract extension from file stub - last string after .
//=========================================================================
string extractfileext(const string & filename) {
    string ext="";
    string::size_type pos = filename.find_last_of(".");
    if (string::npos != pos ) { 
        ext = filename.substr(pos,filename.size()-pos);
    }
    return ext;
}
//=========================================================================
// convert int to binary number string 
//=========================================================================
string int2binary(const int n1) {
    string out(32,'0');
    int n=n1;
    for (int i=31; i>=0; i--) {
        if ((n >> i) & 1) {
            out[i]='1';
        }
    }
    return out;
}


//------------------------------------------------------------------------------
// upperCase  -- converts string to upper case 
//------------------------------------------------------------------------------
string upperCase(const string s) {
    
    // declare rev comp
    string sUpperCase;
    
    // iterate through every char in dna string in reverse order
    for(string::const_iterator iter = s.begin(); iter != s.end(); iter++) {
        char base = *iter;
        
        char upperBase = toupper(base);
        
        // append lower case base to reverse
        sUpperCase += upperBase;
    }
    
    // return rev comp
    return sUpperCase;
}


//------------------------------------------------------------------------------
// lowerCase  -- converts string to lower case 
//------------------------------------------------------------------------------
string lowerCase(const string s) {
    
    // declare rev comp
    string sLowerCase;
    
    // iterate through every char in dna string in reverse order
    for(string::const_iterator iter = s.begin(); iter != s.end(); iter++) {
        char base = *iter;
        
        char lowerBase = tolower(base);
        
        // append lower case base to reverse
        sLowerCase += lowerBase;
    }
    
    // return rev comp
    return sLowerCase;
}

//------------------------------------------------------------------------------
// Titlecase  -- converts string to Title case 
//------------------------------------------------------------------------------
string Titlecase(const string s) {
    
    // declare rev comp
    string sTitlecase;
    
    
    // iterate through every char in dna string in reverse order
    for(string::const_iterator iter = s.begin(); iter != s.end(); iter++) {
        char base = *iter;
        
        char Titlecase = tolower(base);
        if (sTitlecase.size()==0) {
            Titlecase = toupper(base);
        }
        
        // append lower case base to reverse
        sTitlecase += Titlecase;
    }
    
    // return rev comp
    return sTitlecase;
}

// string trim utility functions
string trim_right (const string & s,  string t )
{ 
    string d (s); 
    string::size_type i (d.find_last_not_of (t));
    if (i == string::npos) {
        return "";
    } else {
        return d.erase (d.find_last_not_of (t) + 1) ;  
    }
}  // end of trim_right

string trim_left (const string & s, string  t) 
{ 
    string d (s); 
    return d.erase (0, s.find_first_not_of (t)) ; 
}  // end of trim_left

string trim (const string & s, string  t)
{ 
    string d (s); 
    return trim_left (trim_right (d, t), t) ; 
}  // end of trim


//=========================================================================
// extract text following Bam tag 
//=========================================================================
string extractBamTag(const string & line, const string & tag)
{
    string sout="";
    string s1(line);
    string::size_type pos = s1.find(tag);
    if (string::npos != pos ) { 
        //pos=pos+3;
        sout = s1.substr(pos);
        pos = sout.find("\n");
        if (string::npos != pos ) { 
            string s2=sout;
            sout = s2.substr(0,pos); 
        }        
    }
    return sout;
}

//=========================================================================
// extract text following Bam tag 
//=========================================================================
string extractBamSubTag(const string & line, const string & tag)
{
    string sout="";
    string s1(line);
    string::size_type pos = s1.find(tag);
    if (string::npos != pos ) { 
        pos=pos+tag.size();
        sout = s1.substr(pos);
        pos = sout.find("\t");
        if (string::npos != pos ) { 
            string s2=sout;
            sout = s2.substr(0,pos); 
        }        
    }
    return sout;
}


#endif
