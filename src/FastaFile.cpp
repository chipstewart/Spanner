/*
 *  FastaFile.cpp
 *
 *  Created by Chip Stewart on 11/5/07.
 *  Copyright 2007 Boston College. All rights reserved.
 *
 */
#include "FastaFile.h"

struct ToLower
{
  char operator() (char c) const  { return std::tolower(c); }
};

struct ToUpper
{
  char operator() (char c) const  { return std::toupper(c); }
};
  
//---------------------------------------------------------------------------------
// FastaObj is the container class for a Fasta file sequence 
//---------------------------------------------------------------------------------
// constructor
FastaObj::FastaObj(const string & fn, const string & sn)
{
	this->FastaFile=fn;
  this->numberSeq=0;                                 //  number of sequences
  if (fn.length()==0) {
     return;
  }
  //-------------------------------------------------------------------------------
  // if first argument is "ace", then forgo loading ace file here
  //-------------------------------------------------------------------------------
  if (fn=="ace") {
     return;
  }
  // boost regex matches
  //boost::smatch match,match2;
  string match, match2;

  // Patterns to match
  string patternFastaHeader("^>(\\s*\\S+\\s*.*)$");
  string patternFastaName("^\\s*(\\S+)");
  //boost::regex patternFastaHeader("^>(\\s*\\S+\\s*.*)$");
  //boost::regex patternFastaName("^\\s*(\\S+)");
  
  // unregistered flag to indicate if dna still has to be registered with sequence
  bool registered = true;

  // select seqName flag to indicate if this sequence is to be registered with sequence
  bool selected = false;

  // sequence counter
  int seqCount = 0;

  // sequence name
  string seqName;

  // sequence dna
  string seqDna;

  // sequence header
  string seqHeader;

  // input line
  string line;


  //----------------------------------------------------------------------------
  // input FASTA DNA file
  //----------------------------------------------------------------------------

  // open
  ifstream dnaIn(this->FastaFile.c_str(), ios::in);
  
  if (!dnaIn) {
    cerr << "Unable to open file: " << this->FastaFile << endl;
    exit(1);
  }
  
  while (getline(dnaIn, line)) {
    
    // header line (long format): register previous sequence and start new
    // if (boost::regex_search(line, match, patternFastaHeader)) {
	if (RE2::FullMatch(line.c_str(),patternFastaName.c_str(),&match) ) {
      
      // if previous sequence info not yet registred, register it
      if ( (!registered)&(selected) ) {
	
        // add sequence
        string sn1=seqName;
        transform(sn1.begin(), sn1.end(), sn1.begin(), ToUpper());
        this->seq[sn1]=seqDna;	
        this->seqNames.push_back(sn1);
      }
      
      // increment read count
      seqCount++;
    
      // reset seqDna
      seqDna = "";
      
      // retreive info for new sequence
      seqHeader = match;
      
      // parse out sequence name
      seqName = seqHeader;
		
	  //if (boost::regex_search(seqHeader, match2, patternFastaName)) {     
	  if (RE2::FullMatch(seqHeader.c_str(),patternFastaName.c_str(),&match2) ) {
          seqName = match2;
      }
      
      // select sequence by name
      string sn1=seqName;
      string sn2=sn;      
      transform(sn1.begin(), sn1.end(), sn1.begin(), ToUpper());
      transform(sn2.begin(), sn2.end(), sn2.begin(), ToUpper());
      selected = true;
      if (sn.length()>0) {
        selected = (sn1==sn2);
      }
      // set registered read info flag
      registered = false;
    }
    
    else if (selected) {
      seqDna += line;
    }
  }
  
  // register last DNA entry if necessary
   if ( (!registered)&(selected) ) {
    
    // add sequence
    string sn1=seqName;
    std::transform(sn1.begin(), sn1.end(), sn1.begin(), ToUpper());
    this->seq[sn1]=seqDna;	
    this->seqNames.push_back(sn1);
    //this->seq[seqName]=seqDna;
    //this->seqNames.push_back(seqName);
   
  }
  this->numberSeq = this->seq.size();  //  number of sequences
  dnaIn.close();
}

void FastaObj::addSeq(const string & contigName1, const string & seqDna)
{
  if (seqDna.length()==0) {
     return;
  }
  this->seq[contigName1]=seqDna;	
  this->seqNames.push_back(contigName1);
  this->numberSeq = this->seq.size();  //  number of sequences
}

//========================================================================
// fetch a sequence of name s
//========================================================================
string FastaObj::getSeq(const string & s) 
{
  return this->seq[s];
}
//========================================================================
// fetch a sequence of name s at p of length n
//========================================================================
string FastaObj::getSubSeq(const string & s, size_t p, size_t n) 
{
  return this->seq[s].substr(p,n);
}
//========================================================================
// fetch Number of sequences 
//========================================================================
int FastaObj::getNumberSeq()
{
   return this->numberSeq;
}
size_t FastaObj::getSeqLength(const string & s) 
{
  return this->seq[s].size();
}   
