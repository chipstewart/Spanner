/*
 *  RunControlParameterFile.h
 *  cnv
 *
 *  Created by Chip Stewart on 10/24/07.
 *  Copyright 2007 BC. All rights reserved.
 *
 */
#ifndef RUNCONROLPARAMETERS_H
#define RUNCONROLPARAMETERS_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <limits.h>
// "boost" regular expression library
// #include <boost/regex.hpp>
// 
//"RE2" regular expression library  (replaces boost...)
#include <re2/re2.h>
using namespace re2;

#include "Function-Generic.h"
#include "Histo.h"

using std::string;

class RunControlParameters
{
	friend ostream &operator<<(ostream &, const RunControlParameters &);
public:
	RunControlParameters();	// default constructor
	RunControlParameters( const string &);	// constructor
	RunControlParameters&operator=(const RunControlParameters &rhs);
	string getParameterFile() const;       // return file name
	double getFragmentLengthWindow() const;// return fragment window in %
	double getFragmentLengthLo() const;    // return low fragment length
	double getFragmentLengthHi() const;    // return high fragment length
	double getFragmentLength() const;      // nominal  fragment length
	int getMinClustered() const;           // minumum reads in cluster   
	int getMinClusteredRetro() const;      // minumum reads in cluster for insertions   
	int getClusteringLength() const;       // cluster length parameter    
	int getClusteringLengthRetro() const;  // cluster length parameter for insertions    
	double getDepthLo() const;             // return low depth of coverage
	double getDepthHi() const;             // return high depth of coverage
	int getMaxFragments() const;           // return max Fragment count
	string getAllowContigsRegex() const;   // return regexp for allowed contigs 
	
	int getDbg() const;                    // get debug level
	void setParameterFile(const string &); // set file name
	void setFragmentLengthWindow(const double);  // set fragment window in %
	void setFragmentLengthLo(const double);// set low fragment length
	void setFragmentLengthHi(const double);// set high fragment length
	void setFragmentLength(const double);  // nominal fragment length
	void setMinClustered(const int);       // set minumum reads in cluster
	void setMinClusteredRetro(const int);  // set minumum reads in cluster for insertions
	void setClusteringLength(const int);	   // set cluster length parameter   
	void setClusteringLengthRetro(const int);// set cluster length parameter  for insertions   
	void setDepthLo(const double);         // set low depth of coverage
	void setDepthHi(const double);         // set high depth of coverage
	void setMaxFragments(const int);       // set max fragments to load 
	void setAllowContigsRegex(const string &); // set allowed list of contigs with names / regex
	
	void setDbg(const int);                // set debug level
	void print() const;                    // print RCP object
	string getReadSetRegex() const;        // Pair Set Regex string
	void setReadSetRegex(const string &);	
	string getReadIDRegex() const;         // Pair ID Regex string
	void setReadIDRegex(const string &);	
	string getReadNumberRegex() const;     // Pair Number Regex string
	void setReadNumberRegex(const string &);	
	string getReadEnd1() const;            // Pair Number 1 swap string
	string getReadEnd2() const;            // Pair Number 2 swap string
	void setReadEnds(const string &);	
	string getReferenceFastaFile() const;         // return file name
	void  setReferenceFastaFile(const string &);  // set ReferenceFastaFile Fastafile name
	int getReadPairSenseConfig() const;     // True if DNA sense if flipped within a paired-end
	void setReadPairSenseConfig(const int); 
	bool getPrintTextOut() const;          // print text output files
	void setPrintTextOut(const bool); 
	bool getDoInsertions() const;          // run Insertion detector
	void setDoInsertions(const bool); 
	bool getDoScanOnly() const;            // limit processing to scan only
	void setDoScanOnly(const bool); 
	string getOutputDir() const;           // output dir name
	void  setOutputDir(const string &);  
	string getPrefix() const;              // output file name prefix
	void  setPrefix(const string &);  
	C_HistoGroups getHistoGroups() const;     
	void setHistoGroups(const C_HistoGroups &); 
	//HistObj getFragHist() const;               
	//void setFragHist(const HistObj &); 
	//HistObj getReadLengthHist() const;          
	//void setReadLengthHist(const HistObj &); 
	//void setFragmentLengthLimits();
	
	string getAnchorfile() const;             // output file name prefix
	void  setAnchorfile(const string &);  
	bool getWriteDepth() const;               // write depth of coverage files
	void setWriteDepth(const bool); 
	void setMinLength(const int);             // set minumum event length
	int getMinLength() const;                 // set minumum event length
	int getDepthBinSize() const;              // Bin size for Depth coverage read counts
	void setDepthBinSize(const int);          // Bin size for Depth coverage read counts
	int getMaxMisMatches() const;                // maximum number of mismatches 
	void setMaxMisMatches(const int) ;      
	float getFractionMaxMisMatches() const;      // fractional maximum number of mismatches 
	void setFractionMaxMisMatches(const float) ;      
	int getMinReadLength() const;                // minimum read length 
	void setMinReadLength(const int) ;      
	string getRDreferenceFile() const;           // minimum read length 
	void setRDreferenceFile(const string &) ;      
	int getRDminMedian() const;                  // minimum median read counts in a bin 
	void setRDminMedian(const int) ;      
	float getRDminExp() const;                   // minimum expected read counts in a bin 
	void setRDminExp(const float) ;       
	bool getDoMasking() const;               // repeat masking 
	void setDoMasking(const bool) ;         
	int getRDmaxGap() const;                 // Maximum gap in contiguous bins for 1 CNV 
	void setRDmaxGap(const int) ;      
	int getRDallow() const;                  // Number of allowed CNV per chromsome w/o chisq penalty
	void setRDallow(const int) ;      
	int getRDminbin() const;                 // Minium number of coverage bins in a CNV
	void setRDminbin(const int) ;
	int getQmin() const;                     // Minium number of coverage bins in a CNV
	void setQmin(const int) ;
	int getCNslosh() const;                  // tight/loose PE CNV selection cut
	void setCNslosh(const int) ;
	int getDupRemove() const;                // duplicate fragment removal 
	void setDupRemove(const int) ;
	//bool getDoRepeatCheck() const;           // keep track of repeat locations  
	//void setDoRepeatCheck(const bool) ;
	vector<string> getMobileElements() const;// mobile elements
	void setMobileElements(const vector<string> &) ;      
	string getMobiMaskFile() const;          // mask files
	void setMobiMaskFile(const string &) ;      
	int getBamZA() const;											// require ZA info from bam file 
	void setBamZA(const int) ;
	int getSpannerMode() const;               // Spanner processing mode (scan, build, detect)
	void setSpannerMode(const int); 
	
  // Regex Fragment Length Window 
  string spatternFLWIN;
  // Regex Fragment Length Low 
  string spatternFLLO;
  // Regex Fragment Length High 
  string spatternFLHI;
  // Regex Depth Low 
  string spatternDPLO;
  // Regex Fragment Length High 
  string spatternDPHI;
  // Regex Fragment Length Nominal
  string spatternFL;
  // Regex Min Clustered size
  string spatternMinClus;
  // Regex Min Clustered size for insertions
  string spatternMinClusIns;
  // Regex Clustering Length
  string spatternClusLen;
  // Regex Clustering Length for Insertions
  string spatternClusLenIns;
  // Regex Read Set Label in Read Name
  string spatternReadSet; 
  // Regex Read ID in Read Name
  string spatternReadID;
  // Regex Read Number in Read Name
  string spatternReadNumber;
  // Regex Read End swap string
  string spatternReadEnds;
  // Regex Paired Read default orientiation (short solexa false, long true) 
  string spatternReadPairSenseConfig;
  // Regex do insertion detection
  string spatternDoInsertion;
  // ReferenceFastaFile fasta file full path name
  string spatternReferenceFastaFile;
  // Output area
  string spatternOutputDir;
  // Output file prefix
  string spatternPrefix;
  // specify regexp for allowed contig names
  string spatternAllowContigsRegex;
  // specify maximum number of fragments to process
  string spatternMaxFragments;
  // Regex do scan only
  string spatternDoScanOnly;
  // Regex print text output
  string spatternPrintTextOut;
  // Regex write depth
  string spatternWriteDepth;
  // anchor info file name
  string spatternAnchorfile;                   
  // Regex Min Event Length
  string spatternMinLength;  
  // Regex maximum number of mismatches 
  string spatternMaxMisMatches;  
  // Regex fractional maximum number of mismatches 
  string spatternFractionMaxMisMatches;  
  // Regex minimum read length 
  string spatternMinReadLength;  
  // Regex minimum read length 
  string spatternRDreferenceFile;  
  // Regex minimum median read counts per bin 
  string spatternRDminMedian;  
  // Regex minimum expected read counts per bin 
  string spatternRDminExp;  
  // Regex Depth bin size 
  string spatternDepthBinSize;
  // Regex Masking
  string spatternDoMasking;
  // Regex Maximum gap in contiguous bins for 1 CNV
  string spatternRDmaxGap;                       
  // Regex Number of allowed CNV per chromsome w/o chisq penalty
  string spatternRDallow;                         
  // Regex Minium number of coverage bins in a CNV
  string spatternRDminbin;
  // Regex Minium map Q
  string spatternQmin;
  // Regex PE CNV slosh cut
  string spatternCNslosh;
  // Regex duplicate frag removal
  string spatternDupRemove;
  // Regex flag for do repeat locations 
  string spatternDoRepeatCheck;   
  // Regex mobile elements
  string spatternMobileElements;  
  // Regex mobile mask bed file 
  string spatternMobiMaskFile;  
	// Regex histoGroups 
	string spatternHistoGroups;  	
	// Regex BamZA 
	string spatternBamZA;  
	// Regex SpannerMode 
	string spatternSpannerMode;  

	
private:
  string ParameterFile;                // file name
  double FragmentLengthWindow;         // fragment length window in % 
  double FragmentLengthLo;             // low fragment length  
  double FragmentLengthHi;             // high fragment length
  double FragmentLength;               // nominal fragment length
  int MinClustered;                    // minumum reads in cluster
  int MinClusteredRetro;               // minumum reads in cluster for insertions
  int ClusteringLength;                // return cluster length parameter
  int ClusteringLengthRetro;           // return cluster length parameter for insertions
  double DepthLo;                      // threshold low depth of coverage to consider for CNV
  double DepthHi;                      // threshold high depth of coverage to consider for CNV
  int dbg;                             // debug level 
  int MaxFragments;                    // stop loading alignments after this many reads from each end
  string ReadSetRegex;                 // Pair Set Name Regex string
  string ReadIDRegex;                  // Pair ID Regex string
  string ReadNumberRegex;              // Pair Number Regex string 
  string ReadEnd1;                     // Pair End 1 string 
  string ReadEnd2;                     // Pair End 2 string 
  string AllowContigsRegex;            // only process contigs with names consistent with this regex
  int ReadPairSenseConfig;              // True if DNA sense if flipped within a paired-end
  string ReferenceFastaFile;           // Fasta file with masked bases
  bool doInsertions;                   // switch for Insertion detector
  bool doScanOnly;                     // switch for scanning input files for fragment length only
  bool doPrintTextOut;                 // switch for print text output
  //bool doBuild;                        // build span files from parsed alignment files
  bool writeDepth;                     // write depth of coverage files
  string OutputDir;                    // name of output directory
  string Prefix;                       // prefix in name of output files   
  C_HistoGroups HistoGroups;           // template fragment length distribution
  int MinLength;                       // minumum event length 
  string Anchorfile;                   // anchor info file name  
  int MaxMisMatches;                   // maximum number of mismatches 
  float FractionMaxMisMatches;         // fractional maximum number of mismatches 
  int MinReadLength;                   // minimum read length 
  string RDreferenceFile;              // file name of alignabilty depth info
  int DepthBinSize;                    // Depth bin size
  int RDminMedian;                     // Minium nonzero median read counts per bin (RD rebinning parameter)
  float RDminExp;                      // Minium expected read counts per bin (RD clean parameter)
  bool doMasking;                      // switch for Reference Repeat Masking
  int RDmaxGap;                        // Maximum gap in contiguous bins for 1 CNV
  int RDallow;                         // Number of allowed CNV per chromsome w/o chisq penalty
  int RDminbin;                        // Minium number of coverage bins in a CNV
  int Qmin;                            // Minium map Q for unique mapped read
  int CNslosh;                         // CN buffer added to PE selection cut (del:CN < 2+slosh, dup: CN > 2-slosh)   
  int DupRemove;                       // remove duplicate fragments (0=no, 1=yes, 2=really yes...)  
  bool doRepeatCheck;                  // keep track of repeat locations 
  vector<string> MobileElements;       // list of mobile elements to detect insertions
  string MobiMaskFile;                 // Mask file template (*) for element name 
  int BamZA;                           // Require ZA tag in bam file
	int SpannerMode;                     // SpannerMode (0=scan, 1=build...)
	
  // parameter file strings
  // list of stuff to trim at ends of parameter strings 
  string SPACES;   
  // Regex split by "/"
  string patternSplit;
  // whatever was matched
  string match;
}; // end class 

#endif


