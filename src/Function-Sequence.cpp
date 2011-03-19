//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Function-Sequence
// DNA sequence manipulation methods 
// Copyright 2006, 2007 Gabor T. Marth, Boston College
// All rights reserved
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef FUNCTION_SEQUENCE_CPP
#define FUNCTION_SEQUENCE_CPP

//#include <iostream>
//#include <fstream>
//#include <string>
//#include <vector>
//#include <map>
//#include <iterator>
//#include <algorithm>

// "boost" regular expression library
// #include <boost/regex.hpp>

#include "Function-Sequence.h"

/*
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
*/

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// static variable declarations and initializations
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// unpadDna -- unpads dna sequence 
//------------------------------------------------------------------------------
string unpadDna(const string dna) {
  
  // regular expression variables
  // boost::smatch what;
  // boost::regex pattern("[\\*|\\-]");
	string what;
	string pattern("[\\*|\\-]");
	string nullString = "";

  // unpadded dna output
  string dnaUnpadded = dna;
      
  // replace the occurence of every "*" or "-" with ""
  //	while (boost::regex_search(dnaUnpadded, what, pattern)) {
  //		string newDna = boost::regex_replace(dnaUnpadded, pattern, nullString);
   while (RE2::FullMatch(dnaUnpadded.c_str(),pattern.c_str(),&what)) {	   
	string newDna = dnaUnpadded;
	RE2::Replace(&newDna,pattern.c_str(), nullString.c_str());
    dnaUnpadded = newDna;
  }

  // return unpadded dna
  return dnaUnpadded;
}

//------------------------------------------------------------------------------
// makeBaseMap -- makes base map between padded and unpadded dna 
//------------------------------------------------------------------------------
vector< vector<int> > makeBaseMap(const string dna) {
  
  // make mapping between aligned (padded) dna and dna sequence
  vector<int> basemap_gapped_ungapped, basemap_ungapped_gapped;
  int gappedPos = 1;
  int ungappedPos = 1;
  int alPos = 1;
  int seqPos = 1;
  for(string::const_iterator iter = dna.begin(); iter != dna.end(); iter++) {
    char base = *iter;
    
    // if pad character
    if (base == '*' or base == '-') {
      basemap_gapped_ungapped.push_back(-1);
      gappedPos++;
      alPos++;
    }
    else {
      basemap_gapped_ungapped.push_back(ungappedPos);
      basemap_ungapped_gapped.push_back(gappedPos);
      gappedPos++;
      ungappedPos++;
      alPos++;
      seqPos++;
    }
  }

  // output vector
  vector< vector<int> > maps;
  
  // add maps to vector
  maps.push_back(basemap_gapped_ungapped);
  maps.push_back(basemap_ungapped_gapped);

  // return maps
  return maps;
}

//------------------------------------------------------------------------------
// padDna -- pads dna sequence according to base map
//------------------------------------------------------------------------------
string padDna(const string dna, const vector<int> basemapGappedUngapped) {

  // define output 
  string dnaPadded;
  
  // pad 
  for (int i=0; i<int(basemapGappedUngapped.size()); i++) {
    int pos = basemapGappedUngapped[i];

    if (pos == -1) {
      dnaPadded += "*";
    }
    else {
      dnaPadded += dna.substr(pos-1, 1);
    }
  }

  // return
  return dnaPadded;
}

//------------------------------------------------------------------------------
// padQual -- pads base quality sequence according to base map
//------------------------------------------------------------------------------
vector<short> padQual(const vector<short> qual, const vector<int> basemapGappedUngapped) {

  // define output 
  vector<short> baseQualPadded;
  
  // to each base assign the quality value of the closest unpadded base on the left
  // (if the base is not a pad its own quality value is assigned)
  int prevQ = 0;
  for (int i=0; i<int(basemapGappedUngapped.size()); i++) {
    int pos = basemapGappedUngapped[i];
    if (pos == -1) {
      baseQualPadded.push_back(prevQ);
    }
    else {
      short q = qual[pos-1];
      baseQualPadded.push_back(q);
      prevQ = q;
    }
  }

  // go back and take the minimum with the right neighbor's quality value
  prevQ = 0;
  for (int i=basemapGappedUngapped.size()-1; i>=0; i--) {
    int pos = basemapGappedUngapped[i];
    if (pos == -1) {
      if (prevQ < baseQualPadded[i]) {
	baseQualPadded[i] = prevQ;
      }
    }
    else {
      short q = qual[pos-1];
      prevQ = q;
    }
  }

  // return
  return baseQualPadded;
}

//------------------------------------------------------------------------------
// revCompDna -- reverse complements dna sequence 
//------------------------------------------------------------------------------
string revCompDna(const string dna) {
  
  // declare rev comp
  string dnaRevComp;

  // iterate through every char in dna string in reverse order
  for(string::const_reverse_iterator iter = dna.rbegin(); iter != dna.rend(); iter++) {
    char base = *iter;

    switch (base) {
    case 'a': 
      base = 't';
      break;
    case 'A': 
      base = 'T';
      break;
    case 'c': 
      base = 'g';
      break;
    case 'C': 
      base = 'G';
      break;
    case 'g': 
      base = 'c';
      break;
    case 'G': 
      base = 'C';
      break;
    case 't': 
      base = 'a';
      break;
    case 'T': 
      base = 'A';
      break;
    }

    // append complemented base to reverse
    dnaRevComp = dnaRevComp + base;
  }

  // return rev comp
  return dnaRevComp;
}

//------------------------------------------------------------------------------
// lowerCaseDna -- turns dna sequence into lower case 
//------------------------------------------------------------------------------
string lowerCaseDna(const string dna) {
  
  // declare rev comp
  string dnaLowerCase;

  // iterate through every char in dna string in reverse order
  for(string::const_iterator iter = dna.begin(); iter != dna.end(); iter++) {
    char base = *iter;

    char lowerBase = tolower(base);

    // append lower case base to reverse
    dnaLowerCase += lowerBase;
  }

  // return rev comp
  return dnaLowerCase;
}

//------------------------------------------------------------------------------
// upperCaseDna -- turns dna sequence into upper case 
//------------------------------------------------------------------------------
string upperCaseDna(const string dna) {
  
  // declare rev comp
  string dnaUpperCase;

  // iterate through every char in dna string in reverse order
  for(string::const_iterator iter = dna.begin(); iter != dna.end(); iter++) {
    char base = *iter;

    char upperBase = toupper(base);

    // append lower case base to reverse
    dnaUpperCase += upperBase;
  }

  // return rev comp
  return dnaUpperCase;
}

//------------------------------------------------------------------------------
// gapMapFromDna -- makes a map of gaps in padded dna
//------------------------------------------------------------------------------
map<int, int, std::less<int> > gapMapFromDna(const string dna) {
  
  // output gap map
  map<int, int, std::less<int> > gapMap;
  
  // variables
  int p = 0;
  int gapStart = 0;
  int gapLength = 0;
  bool inGap = false;

  // extract gaps
  for(string::const_iterator iter = dna.begin(); iter != dna.end(); iter++) {
    char b = *iter;
    
    // if pad character
    if (b == '*' or b == '-') {

      // if already in gap, increment gapLength
      if (inGap) {
	gapLength++;
      }

      // otherwise, start gap
      else {
	gapStart = p;
	inGap = true;
	gapLength = 1;
      }
    }

    // if not pad
    else {

      // if gap just ended, register it
      if (inGap) {
	inGap = false;
	gapMap[gapStart] = gapLength;
      }
      p++;
    }
  }

  // return maps
  return gapMap;
}

//------------------------------------------------------------------------------
// unpaddedPosMap -- makes a map from padded to unpadded dna position:
//                   positions of pads are assigned the value 0
//------------------------------------------------------------------------------
vector<int> unpaddedPosMap(const string dna) {
  
  // output map
  vector<int> posMap;

  // variables
  int zero = 0;
  int up = 0;
  for(string::const_iterator iter = dna.begin(); iter != dna.end(); iter++) {
    char b = *iter;
    
    // if pad character
    if (b == '*' or b == '-') {
      posMap.push_back(zero);
    }
    else {
      up++;
      posMap.push_back(up);
    }
  }

  // return 
  return posMap;
}

//------------------------------------------------------------------------------
// unpaddedPosMapBegin -- makes a map from padded to unpadded dna position:
//                        positions of pads are assigned the base towards the beginning
//------------------------------------------------------------------------------
vector<int> unpaddedPosMapBegin(const string dna) {
  
  // output map
  vector<int> posMap;

  // variables
  int up = 0;
  for(string::const_iterator iter = dna.begin(); iter != dna.end(); iter++) {
    char b = *iter;
    
    // if pad character
    if (b == '*' or b == '-') {
      posMap.push_back(up);
    }
    else {
      up++;
      posMap.push_back(up);
    }
  }

  // return 
  return posMap;
}

//------------------------------------------------------------------------------
// unpaddedPosMapEnd -- makes a map from padded to unpadded dna position:
//                      positons of pads are assigned the base towards the end
//------------------------------------------------------------------------------
vector<int> unpaddedPosMapEnd(const string dna) {
  
  // output map
  vector<int> posMap;

  // unpad dna
  string dnaUnpadded = unpadDna(dna);

  // variables
  int up = dnaUnpadded.size()+1;
  for(string::const_reverse_iterator iter = dna.rbegin(); iter != dna.rend(); iter++) {
    char b = *iter;
    
    // if pad character
    if (b == '*' or b == '-') {
      posMap.push_back(up);
    }
    else {
      up--;
      posMap.push_back(up);
    }
  }
  // reverse order
  std::reverse(posMap.begin(), posMap.end());

  // return 
  return posMap;
}

//------------------------------------------------------------------------------
// paddedPosMap -- makes a map from unpadded to padded dna
//------------------------------------------------------------------------------
vector<int> paddedPosMap(const int length, map<int, int, std::less<int> > gapMap) {
  
  // output map
  vector<int> posMap;

  // variables
  int pp = 1;
  for (int up=1; up<=length; up++) {
    posMap.push_back(pp);
    pp++;

    if (gapMap.count(up) > 0) {
      pp += gapMap[up];
    }
  }

  // return 
  return posMap;
}

//------------------------------------------------------------------------------
// padDnaGapMap -- pad dna sequence
//------------------------------------------------------------------------------
string padDnaGapMap(const string dna, map<int, int, std::less<int> > gapMap) {
  
  // output dna
  string paddedDna("");

  // variables
  int up = 0;
  for(string::const_iterator iter = dna.begin(); iter != dna.end(); iter++) {
    char b = *iter;

    up++;
    paddedDna += b;
    if (gapMap.count(up) > 0) {
      for (int i=0; i<gapMap[up]; i++) {
	paddedDna += "*";
      }
    }
  }

  // return 
  return paddedDna;
}

//------------------------------------------------------------------------------
// padQualGapMap -- pad qual sequence
//------------------------------------------------------------------------------
vector<short> padQualGapMap(const vector<short> qual, map<int, int, std::less<int> > gapMap) {
  
  // output dna
  vector<short> paddedQual;

  // variables
  int up = 0;
  short prevQ = 0;
  for(vector<short>::const_iterator iter = qual.begin(); iter != qual.end(); iter++) {
    short q = *iter;

    // increment unpadded pos
    up++;
    
    // register qual value at this position
    paddedQual.push_back(q);

    // if there are gaps after this position, register average between last and next 
    // base's qual value
    if (gapMap.count(up) > 0) {
      prevQ = q;
      short nextQ = qual[up];
      short gapQ = prevQ;
      if (nextQ < prevQ) {
        gapQ = nextQ;
      }
      for (int i=0; i<gapMap[up]; i++) {
        paddedQual.push_back(nextQ);
      }
    }

    // update previous qual value
    prevQ = q;
  }

  // return 
  return paddedQual;
}

//------------------------------------------------------------------------------
// dnaGenotypes -- returns a list of dna genotypes (haploid or diploid)
//------------------------------------------------------------------------------
vector<string> dnaGenotypes(const int mult, const bool includePad) {
  
  // complain in bomb if multiplicity is not 1 or 2
  if (mult != 1 and mult !=2) {
    clog << "genotypes of nonsensical multiplicity requested. terminating." << endl;
    exit(1);
  }

  // define genotype vector
  vector<string> Genotypes;

  // haploid alleles
  if (mult == 1) {
    Genotypes.push_back("a");
    Genotypes.push_back("c");
    Genotypes.push_back("g");
    Genotypes.push_back("t");
    if (includePad) {
      Genotypes.push_back("*");
    }
  }
  
  // diploid alleles
  else if (mult == 2) {
    Genotypes.push_back("aa");
    Genotypes.push_back("ac");
    Genotypes.push_back("ag");
    Genotypes.push_back("at");
    Genotypes.push_back("cc");
    Genotypes.push_back("cg");
    Genotypes.push_back("ct");
    Genotypes.push_back("gg");
    Genotypes.push_back("gt");
    Genotypes.push_back("tt");
    if (includePad) {
      Genotypes.push_back("a*");
      Genotypes.push_back("c*");
      Genotypes.push_back("g*");
      Genotypes.push_back("t*");
      Genotypes.push_back("**");
    }
  }
  
  // return genotype vector
  return Genotypes;
}

//------------------------------------------------------------------------------
// dnaGenotypeIsValid -- returns bool according to whether or not input genotype
//                       is valid
//------------------------------------------------------------------------------
bool dnaGenotypeIsValid(const int mult, const bool includePad, const string genotype) {
  
  // complain in bomb if multiplicity is not 1 or 2
  if (mult != 1 and mult !=2) {
    clog << "genotypes of nonsensical multiplicity requested. terminating." << endl;
    exit(1);
  }

  // valid flag
  bool valid = false;

  // haploid alleles
  if (mult == 1) {
    if (genotype == "a" or genotype == "c" or genotype == "g" or genotype == "t" or (includePad and (genotype == "*"))) {
      valid = true;
    }
  }
  
  // diploid alleles
  else if (mult == 2) {
    if (
	genotype == "aa" or
	genotype == "ac" or
	genotype == "ag" or
	genotype == "at" or
	genotype == "cc" or
	genotype == "cg" or
	genotype == "ct" or
	genotype == "gg" or
	genotype == "gt" or
	genotype == "tt" or
	(includePad and (
			 genotype == "a*" or
			 genotype == "c*" or
			 genotype == "g*" or
			 genotype == "t*" or
			 genotype == "**")
	 )
	) {
      valid = true;
    }
  }
  
  // return validity
  return valid;
}

//------------------------------------------------------------------------------
// getDnaGenotypeIndex -- returns the int index of the input genotype
//                     -- return the int index of the complement if complement is true
//------------------------------------------------------------------------------
int getDnaGenotypeIndex(const int mult, string genotype, bool complement) {
  
  // complain in bomb if multiplicity is not 1 or 2
  if (mult != 1 and mult !=2) {
    clog << "genotypes of nonsensical multiplicity requested. terminating." << endl;
    exit(1);
  }

  // define output index
  int index = -1;

  // if non-complemented
  if (! complement) {
    if (mult == 1) {
      if (genotype == "a") {index = 0;}
      else if (genotype == "c") {index = 1;}
      else if (genotype == "g") {index = 2;}
      else if (genotype == "t") {index = 3;}
    }
    else if (mult == 2) {
      if (genotype == "aa") {index = 0;}
      else if (genotype == "ac") {index = 1;}
      else if (genotype == "ag") {index = 2;}
      else if (genotype == "at") {index = 3;}
      else if (genotype == "cc") {index = 4;}
      else if (genotype == "cg") {index = 5;}
      else if (genotype == "ct") {index = 6;}
      else if (genotype == "gg") {index = 7;}
      else if (genotype == "gt") {index = 8;}
      else if (genotype == "tt") {index = 9;}
    }
  }
  else {
    if (mult == 1) {
      if (genotype == "a") {index = 3;}
      else if (genotype == "c") {index = 2;}
      else if (genotype == "g") {index = 1;}
      else if (genotype == "t") {index = 0;}
    }
    else if (mult == 2) {
      if (genotype == "aa") {index = 9;}
      else if (genotype == "ac") {index = 8;}
      else if (genotype == "ag") {index = 6;}
      else if (genotype == "at") {index = 3;}
      else if (genotype == "cc") {index = 7;}
      else if (genotype == "cg") {index = 5;}
      else if (genotype == "ct") {index = 2;}
      else if (genotype == "gg") {index = 4;}
      else if (genotype == "gt") {index = 1;}
      else if (genotype == "tt") {index = 0;}
    }
  }
  return(1);
}

//------------------------------------------------------------------------------
// extractTemplateName -- extracts template name from sequence name
//------------------------------------------------------------------------------
string extractTemplateName (const string sequenceName, const string templateSource) {
  
  // define templateName as sequenceName
  string templateName = sequenceName;

  if (templateSource == "baylor") {

    // whatever was matched
    string match;
    
    // Baylor read naming convention
    string patternBaylor("^(\\S+)\\_(\\S+)\\_(\\S+)\\.scf");
    
    // check pattern for Baylor template
	if (RE2::FullMatch(sequenceName.c_str(),patternBaylor.c_str(),&match) ) {
      templateName = match;
    }
  }
  else if (templateSource == "elegans") {

    // whatever was matched
    string match;
    
    // Solexa read naming convention for elegans strains
    string patternElegans("(\\S+)\\_\\S+\\_\\S+\\_\\S+\\_\\S+\\_\\S+");

    // check pattern for elegans template
	if (RE2::FullMatch(sequenceName.c_str(),patternElegans.c_str(),&match) ) {
      templateName = match;
    }
  }

  else if (templateSource == "drosophila") {

    // whatever was matched
    string match;
    
    // Capillary read naming convention for drosophila validation reads
    string patternDroVal("^(DLAH)");

    // Solexa read naming convention for drosophila strains
    string patternDro454("^(\\S+)\\_");

    // check pattern for drosophila template
	if (RE2::FullMatch(sequenceName.c_str(),patternDroVal.c_str(),&match) ) {
      templateName = match;
    }
    else if (RE2::FullMatch(sequenceName.c_str(),patternDro454.c_str(),&match) ) {
	
      templateName = match;
    }
  }

  else if (templateSource == "multiple") {

    // whatever was matched
    string match;
    
    // Solexa read naming convention for drosophila strains
	string patternMultiple("^(\\S+)\\-");

    // check pattern for template
	if (RE2::FullMatch(sequenceName.c_str(),patternMultiple.c_str(),&match) ) {
      templateName = match;
    }
  }

  else if (templateSource == "single") {
    templateName = "single";
  }

  else if (templateSource == "unknown") {
    templateName = sequenceName;
  }

  // return
  return templateName;
}

//------------------------------------------------------------------------------
// printDna -- prints dna sequence to ostream
//------------------------------------------------------------------------------
void printDna(ostream &out, const string dna, const int L) {
  for (int lineCount = 0; lineCount*L < int(dna.length()); lineCount++) {
    string dnaLine = dna.substr(lineCount*L, L);
    out << dnaLine << endl;
  }
}

//------------------------------------------------------------------------------
// printDnaFasta -- prints dna sequence to ostream in FASTA format
//------------------------------------------------------------------------------
void printDnaFasta(ostream &out, const string header, const string dna, const int L) {
  out << ">" << header << endl;
  for (int lineCount = 0; lineCount*L < int(dna.length()); lineCount++) {
    string dnaLine = dna.substr(lineCount*L, L);
    out << dnaLine << endl;
  }
}

//------------------------------------------------------------------------------
// printQual -- prints qual sequence to ostream
//------------------------------------------------------------------------------
void printQual(ostream &out, const vector<short> qual, const int L) {
  int numberQual = qual.size();
  
  bool first = true;
  int l = 0;
  for (int i=0; i<numberQual; i++) {
    l++;
    if (!first) {
      out << " ";
    }
    first = false;
    out << qual[i];
    if (l>=L and i<numberQual-1) {
      out << endl;
      l=0;
      first = true;
    }
  }
  out << endl;
}

//------------------------------------------------------------------------------
// printQualFasta -- prints qual sequence to ostream in FASTA format
//------------------------------------------------------------------------------
void printQualFasta(ostream &out, const string header, const vector<short> qual, const int L) {
  int numberQual = qual.size();
  
  out << ">" << header << endl;
  bool first = true;
  int l = 0;
  for (int i=0; i<numberQual; i++) {
    l++;
    if (!first) {
      out << " ";
    }
    first = false;
    out << qual[i];
    if (l>=L and i<numberQual-1) {
      out << endl;
      l=0;
      first = true;
    }
  }
  out << endl;
}

//------------------------------------------------------------------------------
// printBpos -- prints bpos sequence to ostream
//------------------------------------------------------------------------------
void printBpos(ostream &out, const vector<int> bpos, const int L) {
  int numberBpos = bpos.size();
  
  bool first = true;
  int l = 0;
  for (int i=0; i<numberBpos; i++) {
    l++;
    if (!first) {
      out << " ";
    }
    first = false;
    out << bpos[i];
    if (l>=L and i<numberBpos-1) {
      out << endl;
      l=0;
      first = true;
    }
  }
  out << endl;
}

//------------------------------------------------------------------------------
// printBposFasta -- prints bpos sequence to ostream in FASTA format
//------------------------------------------------------------------------------
void printBposFasta(ostream &out, const string header, const vector<int> bpos, const int L) {
  int numberBpos = bpos.size();
  
  out << ">" << header << endl;
  bool first = true;
  int l = 0;
  for (int i=0; i<numberBpos; i++) {
    l++;
    if (!first) {
      out << " ";
    }
    first = false;
    out << bpos[i];
    if (l>=L and i<numberBpos-1) {
      out << endl;
      l=0;
      first = true;
    }
  }
  out << endl;
}

//------------------------------------------------------------------------------
// getCigarLength -- returns the total sequence length of a cigar string
//   opt=0: query length
//   opt=1: reference length
//------------------------------------------------------------------------------
int getCigarLength(string cigar, int opt=0) 
{	
	RE2 patternCigar("(\\d+)(\\S)");

	StringPiece sp(cigar);    // Wrap a StringPiece around it

	int len1;
	string op1;
	
	int len=0;
	while (RE2::Consume(&sp,patternCigar, &len1, &op1)) {
		if (opt==0) {
			// query length does not count deletions in reference
			size_t found=op1.find("D");
			if (found==string::npos) 	// skip D's (insertions in query)
				len+=len1;
		} else if (opt==1) {
			// reference length does not count insertions  in reference
			size_t found=op1.find("I");
			if (found==string::npos) 	// skip D's (insertions in query)
				len+=len1;
		}			
	}
    return len;
}

//------------------------------------------------------------------------------
// getCigarLengths -- returns subtotals sequence length within a cigar string
//   LQ: query length
//   LR: reference length
//   MM: mismatch length
//------------------------------------------------------------------------------
bool getCigarLengths(string cigar, int &LQ, int &LR, int &MM) 
{	
	RE2 patternCigar("(\\d+)(\\S)");
	
	StringPiece sp(cigar);    // Wrap a StringPiece around it
	
	int len1;
	string op1;
	
	LQ=0;
	LR=0;
	MM=0;
	
	while (RE2::Consume(&sp,patternCigar, &len1, &op1)) {
		size_t found=op1.find("S");
		if (found!=string::npos) 	// skip S's (clips in query)
			continue;
		// query length does not count deletions in reference
		found=op1.find("D");
		if (found==string::npos) 	// skip D's (insertions in query)
			LQ+=len1;
		// reference length does not count insertions  in reference
		found=op1.find("I");
		if (found==string::npos) 	// skip D's (insertions in reference)
			LR+=len1;
		found=op1.find("M");
		if (found==string::npos) 	// skip M's (matched - substitutions?)
			MM+=len1;
	}
	return (LQ+LR+MM)>0;
}

//------------------------------------------------------------------------------
// getCigarMismatchCount -- returns the total number of mismatches in cigar string
//------------------------------------------------------------------------------
int getCigarMismatchCount(string cigar) 
{
		
	RE2 patternCigar("(\\d+)(\\S)");
	
	StringPiece sp(cigar);    // Wrap a StringPiece around it
	
	int len1;
	string op1;
	
	int mm=0;
	//while (RE2::Consume(&sp,patternCigar.c_str(), &len1, &op1)) {
	while (RE2::Consume(&sp,patternCigar, &len1, &op1)) {
		
		size_t found=op1.find("S");
		if (found!=string::npos) 	// skip S's (clips in query)
			continue;
		found=op1.find("M");
		if (found==string::npos) 		
			mm+=len1;
	}
	return mm;
}


//------------------------------------------------------------------------------
// getMDMismatchCount -- returns the total number of mismatches in MD string
//------------------------------------------------------------------------------
int getMDMismatchCount(string MD) 
{
		
	RE2 patternMD("(\\D+)"); // anything not a digit
	
	StringPiece sp(MD);    // Wrap a StringPiece around it
	
	string se;
	
	int mm=0;
	while (RE2::Consume(&sp,patternMD, &se)) {
		
		size_t len1=se.size();
		mm+=len1;
	}
	return mm;
}





#endif
