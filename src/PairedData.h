/*
 *  PairedReads.h
 *  cnv
 *
 *  Created by Chip Stewart on 10/28/07.
 *  Copyright 2007 Boston College. All rights reserved.
 *
 */
#ifndef PAIREDREADS_H
#define PAIREDREADS_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <list>
#include <math.h>
#include <dirent.h> 
#include <stdio.h> 
#include <time.h>
// private
#include "Spanner.h"
#include "Type-Hash.h"
// "boost" regular expression library
// #include <boost/regex.hpp>
// local includes
#include "RunControlParameterFile.h"
#include "Histo.h"
#include "Function-Generic.h"
#include "Function-Sequence.h"
#include "MosaikAlignment.h"
#include "headerSpan.h"
#include "api/BamMultiReader.h"
#include "SHA1.h"

using namespace std;
using namespace BamTools;
//using namespace MosaikReadFormat;
//namespace fs = boost::filesystem;



//------------------------------------------------------------------------------
// Basic single read-map properties
//------------------------------------------------------------------------------
class C_readmap {
  friend ostream &operator<<(ostream &, const C_readmap &);
  public:
    unsigned int pos;           // position in contig (unpadded)
    unsigned short len;         // length of this read aligment in contig coordinates
    unsigned short anchor;      // anchor index  
    char sense;                 // forward ('F') or reverse complement ('R')
		char q;                     // mapping quality
  	char q2;                    // mapping quality of next best alignment
	  unsigned short nmap;        // number of mappings for this read 
	  unsigned short mm;          // number of mismatches     
  	string mob;                 // two char flag for special contig hit
	  C_readmap();  
    C_readmap(const C_readmap &);  
    C_readmap(unsigned int, unsigned short, unsigned short, char, char,char);
    ~C_readmap(){};  
    C_readmap&operator=(const C_readmap &rhs);
    int operator==(const C_readmap &rhs) const;
    int operator<(const C_readmap &rhs) const;  
}; 

//------------------------------------------------------------------------------
// multiply-mapped read class 
//------------------------------------------------------------------------------
class C_readmaps {
  friend ostream &operator<<(ostream &, const C_readmaps &);
  public:
    C_readmaps();                   // constructor
    C_readmaps(const C_readmaps &);  
	  vector<C_readmap> align;         // vector of all alignments for this read
    ~C_readmaps(){};  
    C_readmaps&operator=(const C_readmaps &rhs);
    char element;
    unsigned int Nalign;
};

//------------------------------------------------------------------------------
// Paired-read class  
//------------------------------------------------------------------------------
class C_pairedread {
 friend ostream &operator<<(ostream &, const C_pairedread &); 
  public:
    C_pairedread();              // constructor
    C_readmaps read[2];          // read[0] or read[1]
    ~C_pairedread(){};  
    unsigned int ReadGroupCode;
		string Name;
  private:
    unsigned int Nend;

};

//------------------------------------------------------------------------------
// anchor info 
//------------------------------------------------------------------------------
class C_anchorinfo {
  friend ostream &operator<<(ostream &, C_anchorinfo &);
  public:
    C_anchorinfo();               // constructor
    C_anchorinfo(const C_anchorinfo &);           // copy constructor
    C_anchorinfo(char);           // AB constructor
    C_anchorinfo(string &);       // constructor - load from anchor info file
    C_anchorinfo(Mosaik::CAlignmentReader &);    // constructor - load from Mosaik file
    //C_anchorinfo(BamReader &);    // constructor - load from Bam  file
    C_anchorinfo(BamMultiReader &);    // constructor - load from Bam  file
    ~C_anchorinfo(){};            // destructor
    int operator==(const C_anchorinfo &rhs) const;
    void anchorlimit(string &);   // use regex to limit anchors to be processed (memory limits)
    bool anchorElements(vector<string> & );  // set anchor elements (for umulti-fragment elements)
    bool anchorElements();        // set default anchor elements 
    unsigned int anchorMinElement();
    void printAnchorInfo(string &);
    char anchorIndex(string &);   // return anchor index given anchor name 
	map<string, unsigned int, less<string> > L;
    vector<char> use; 
    vector<char> element; 
    vector <string> names;
    string source;
    void push_anchor (string &, unsigned int);
};

    
//------------------------------------------------------------------------------
// library info 
//------------------------------------------------------------------------------
class C_libraryinfo {
  friend ostream &operator<<(ostream &, C_libraryinfo &);
  public:
    C_libraryinfo();                          // constructor
    C_libraryinfo(const C_libraryinfo &);     // copy constructor
    Mosaik::ReadGroup Info; 
    /* BAM/Mosiak header info:
    unsigned int MedianFragmentLength;
		unsigned int ReadGroupCode;
		SequencingTechnologies SequencingTechnology;
		string CenterName;
		string Description;
		string LibraryName;
		string PlatformUnit;
		string ReadGroupID;
		string SampleName;
    */   
    // Fragment lengths 
    HistObj fragHist;                    // fragment length distribution
    HistObj readLengthHist;              // read length distribution
    int LM;
    double tailcut;
    int LMlow;
    int LMhigh;    
    // read lengths
    double LR;
    int LRmin;
    int LRmax;
    // redundant fragment fractions
    int NPair;
    int NSingle;
    int NPairRedundant;
    int NSingleRedundant; 

};

//------------------------------------------------------------------------------
// multiple libraries -> C_libraries libmap. 
//------------------------------------------------------------------------------
typedef std::map<unsigned int, C_libraryinfo, std::less<unsigned int> >  C_librarymap;

class C_libraries {
  friend ostream &operator<<(ostream &, C_libraries &);
public:
	C_libraries();                            // constructor
	C_libraries(const C_libraries &);         // copy constructor
	C_libraries(Mosaik::CAlignmentReader &);  // constructor - load from Mosaik file    
	//C_libraries(BamReader  & );               // constructor - load from BAM file    
	C_libraries(BamMultiReader  & );               // constructor - load from BAM file    
	C_libraries(string  &);                   // load constructor from library.span file  (& AB)
	void printLibraryInfo(string &);          // text dump
	void writeLibraryInfo(string &, string &);// binary file output
	int maxLF();                              // maximum library fragment size
	double getTailcut();                      // return tailcut from first library
	void resetFragLimits(double);             // set frag limits LMlow, LMhigh by new tailcut 
	vector<Mosaik::ReadGroup> GetReadGroups(string &); // parse RG info from sam/bam header  
	vector<char> getSequencingTechnology();   // list of platform techologies
	C_librarymap  libmap;                     // map of library info records indexed by unit code
	C_anchorinfo anchorinfo;                  // only one anchor info needed for all libs
	map<string, double, less<string> > readFractionSamples();  // read fractions for each sample
	map<string, unsigned int , less<string> > ReadGroupID2Code;  // read fractions for each sample
	
};  




//  local read pair structure  (mapped fragment contained on one anchor reference)
class C_localpair {
  friend ostream &operator<<(ostream &, const C_localpair &);
  public:
   unsigned int pos;            // position in contig (unpadded)
   int lm;                      // length from F to R ends
   unsigned short anchor;       // anchor index  
   unsigned short len1;         // length of F read aligment 
   unsigned short len2;         // length of R read aligment 
   char orient;                 // orient FR ('') or both forward ('F') or both reverse complement ('R')
   char q1;                     // F mapping quality
   char q2;                     // R mapping quality
   char mm1;                    // F mismatches 
   char mm2;                    // R mismatches
   char constrain;              // marks non-unique ends 0-both u, 1-FU,2-RU,3-both NU
   unsigned int ReadGroupCode;  // library info index
   C_localpair();  
   C_localpair(const C_localpair &);  
   C_localpair(unsigned int, int, unsigned short,unsigned short 
       ,unsigned short, char, char, char,char, char, char,unsigned int);
   C_localpair(C_pairedread &,char);
   ~C_localpair(){};  
   C_localpair&operator=(const C_localpair &rhs);
   int operator==(const C_localpair &rhs) const;
   int operator<(const C_localpair &rhs) const;  
}; 


// cross read pair structure ( mapped fragment across two anchor references)
class C_crosspair {
  friend ostream &operator<<(ostream &, const C_crosspair &);
  public:
   C_readmap read[2];          // position in contig (unpadded)
   unsigned int ReadGroupCode; // library info index
   C_crosspair();  
   C_crosspair(const C_crosspair &);  
   C_crosspair(const C_readmap  &, const C_readmap  &, const unsigned int);
   ~C_crosspair(){};  
   C_crosspair&operator=(const C_crosspair &rhs);
   int operator==(const C_crosspair &rhs) ; //const;
   int operator<(const C_crosspair &rhs) const;  
}; 

//------------------------------------------------------------------------------
// u-multi read pair structure ( fragment with one unique and one multiple map)
//------------------------------------------------------------------------------
class C_umpair {
  friend ostream &operator<<(ostream &, const C_umpair &);
  public:
   C_readmap read[2];          // position in contig (unpadded)
   int nmap;                   // Number of different element subclasses (moblist hits)
   int nmapA;                  // Total number of moblist hits
   int elements;
   unsigned int ReadGroupCode;  // library info index
   C_umpair();  
   C_umpair(const C_umpair &); 
   C_umpair(const C_readmap  &, const C_readmap  &, int, int, unsigned int);
   C_umpair(const C_pairedread  &, const C_anchorinfo & );
   ~C_umpair(){};  
   C_umpair&operator=(const C_umpair &rhs);
   bool constrain(int LMlow, int LMhigh);
   int operator==(const C_umpair &rhs) ; //const;
   int operator<(const C_umpair &rhs) const;  
}; 

//------------------------------------------------------------------------------
// C_read class for single end reads - adds ReadGroupCode to readmap class
//------------------------------------------------------------------------------
class C_singleEnd : public C_readmap {
  public:
    C_singleEnd();  
    C_singleEnd(const C_singleEnd &); 
    C_singleEnd(unsigned int, unsigned short,unsigned short, char, char, char, unsigned int);
    ~C_singleEnd(){};  
    C_singleEnd&operator=(const C_singleEnd &rhs);
    unsigned int ReadGroupCode; 
};

//------------------------------------------------------------------------------
// depth of coverage class
//------------------------------------------------------------------------------
class C_depth {
  friend ostream &operator<<(ostream &, const C_depth &);
  public:
    C_depth();              // constructor
    C_depth(int);           // constructor
    ~C_depth() {};          // destructor
    int pos0;
    int pos1;
    long long light;
    vector<float> n;        // depth of coverage
    int nzbins;
    double nzMedian;    
    StatObj Stats;          // summary stats for depth
    string name;
    unsigned int ReadGroupCode;  // library info index
    int addMarks(const string &,const string &, char); 
    int GCcontent(const string &,const string &); 
    void calcStats();              // default stat calculator
    void calcStats(int,float,float);     // specify bins
 
};

// repeat marker class
class C_marker {
  friend ostream &operator<<(ostream &, const C_marker &);
  public:
    C_marker() {};          // constructor
    C_marker(int);          // constructor
    ~C_marker() {};         // destructor
    vector<bool> x;         // mark bases 
    StatObj Stats;          // summary stats for marker
    string name;
    unsigned int ReadGroupCode;  // library info index
    int addMarks(const string &,const string &, char); 
};


// pair map type
typedef std::map<string, C_headerSpan, std::less<string> >  C_headers;

//contig class
class C_contig {
  friend ostream &operator<<(ostream &, const C_contig &);
  public:
    C_contig() {};                               // default constructor
    C_contig(string &, int, bool);               // constructor
    void calcStats();                            // calculate stats
    void calcLengths();                          // calculate average read length
    void calcDepth(int);                         // calculate read depth of coverage 
    void calcDepth(int,int,int);                 // calculate read depth in region
    void calcFragDepth(int,int);                 // calculate fragment depth of coverage
    void calcFragDepth(int,int,int,int);         // calculate fragment coverage in region
    void calcStarts(int);                        // calculate read starts
    long long countReads(int,int);               // count read start in range 
    int howFar(int,int);                         // end position of region with set number of non-repeat bases
    int setLengthFromAnchor();                   // set contig Length from anchor[contigName].L
    list<C_localpair> localpairs;                // local unique pair  map
    list<C_crosspair> crosspairs;                // cross contig unique pair map
    list<C_umpair> umpairs;                      // unique-multiple pair (more info than uniqueEnd)
    list<C_singleEnd> dangle;                    // dangling partner of missing end of pair
    list<C_singleEnd> singleton;                 // one end of pair unique - other end ambiguous 
    string getContigName() const;                // contig name
    void setContigName(string &) ;               // contig name
    C_depth  read_depth;                         // depth of coverage 
    C_depth  frag_depth;
    C_depth  read_start;                         // depth of read start positions (1/read)
    C_depth  repeat;                             // repeat map based on multiply aligned reads
    C_depth  read_counts;                        // count of reads per bin
    StatObj pairStats;                           // summary stats for local fragment lenths 
    StatObj crossStats;                          // summary stats for cross fragments 
    StatObj dangleStats;                         // summary stats for cross fragments 
    int Length;                                  // contig length;
    double aLR;                                  // mean read length
    double sLR;                                  // stdev read length
    double totalUniqueReads;                     // total count of unique reads 
    double totalRepeatBases;                     // total count of repeat bases 
    double totalNoCovBases;                      // total count of leading & trailing unaccessable bases  
    int  uniquified;                             // uniquify flag (0=not yet, 1=sorted already, 2 unique already ,...)
    C_anchorinfo anchors;
    unsigned short getAnchorIndex(); 
    string setName;             // set name
    void writePairs(string & ); // const // odd that const kills the list writing...;
    void printPairs(string & ); // const // odd that const kills the list writing...;
    void writeEnd(string &, list<C_singleEnd> & ); // const // odd that const kills the list writing...;
    void printEnd(string &,  list<C_singleEnd> & ); 
    void writeCross(string & );
    void printCross(string & );
    void writeMulti(string &, list<C_umpair> &);
    void printMulti(string &, list<C_umpair> &);
    void writeDepth(string &, C_depth &);
    void writeDepth(string &, C_depth &, int binsize, int totbin);
    void writeMarker(string &, C_marker &);
    void writeStats(string & );
    void printStats(string & );
    C_headerSpan loadHeader(fstream &);
    C_depth loadDepth(string & ); 
    C_depth countReads(int binsize);
    C_marker loadMarker(string & ); 
    list<C_singleEnd>   loadEnd(string & );                    
    void loadPairs(string & ); 
    void loadCross(string & );  
    list<C_umpair> loadMulti(string &);  
    //void loadRetroStart(string &);  
    void loadMultipairs(string & );  
    void loadRepeat(string & );  
    //void loadUniqueEnd(string & );  
    void loadDangle(string & );  
    //void loadDangleEnd(string & );  
    void loadSingleton(string & );      
    void sort();
    void uniquify();
    //void uniquifyBam();
    static bool isRedundantPair(const C_localpair &p1, const C_localpair &p2);
    static bool isRedundantCross(const C_crosspair &p1, const C_crosspair &p2);
    static bool isRedundantRead(const C_singleEnd &r1, const C_singleEnd &r2);
    static bool isRedundantMulti(const C_umpair &p1, const C_umpair &p2);
    //static bool isRedundantPairBam(C_localpair &p1, C_localpair &p2);
   private:
    string contigName;                              // contig name
}; // end class 

// pair map type
typedef std::map<string, C_pairedread, std::less<string> >  C_pairedreads;

// contigs type
typedef std::map<string, C_contig, std::less<string> >  C_contigs;

/*
class C_readGroupTags {
	friend ostream &operator<<(ostream &, const C_readGroupTags &);
	public:	
		string sample;
		string library;
		string description;
		string platformUnit;
		string predictedInsertSize;
		string platformTech;
		string daterun;
		string seqcenter;
};
*/

// pair map container class for each set 
class C_set {
  friend ostream &operator<<(ostream &, const C_set &);
public:
	C_set(string &, RunControlParameters &);               // constructor
	C_contigs contig;
	string getSetName() const;
	void  setSetName(string &);
	string getFileName() const;
	void  setFileName(string &);
	StatObj repeatStats;                            // summary stats for repeat multi-mapped reads 
	StatObj lengthStats;                            // summary stats for read lenths 
	StatObj fragStats;                              // summary stats for fragment lenths 
	StatObj qStats;                                 // summary stats for mapping quality 
	StatObj pairCountStats;                         // summary stats for pair counting 0 1 N     
	StatObj pairModelStats;				     		// summary stats for pair orientation models 1-8
	StatObj spanStats;                              // summary stats for span lenths 
	StatObj refStats;                               // summary stats for reference ID 
	
	C_anchorinfo anchors;
	void initContigs();                             // loop over anchor info to init contigs map
	void processRepeat(C_readmaps &, int);
	void processSingleton(C_readmaps &, unsigned int);
	void processPair(C_pairedread &,char);
	int  Fraglength(C_pairedread & );
	int  spanLength(C_pairedread & );
	int  pairModel(C_pairedread &); 	
	C_pairedread resolvePairConstraint(C_pairedread &); 
	void calcStats();
	void sort();
	void uniquify();
	void write();                                   
	void printOut();  
	void calcDepth();                               // calculate depth of coverage
	void writeDepth();                          
	void countReads(int);
	void writeCountReads();  
	unsigned int Nfrag;
	unsigned int Npair;
	RunControlParameters pars;
	unsigned int Npair_11o;
	unsigned char elementMinAnchor;                 // start of element anchor indices...
	C_libraries libraries;
	void calcLibraryRedundancy(C_libraryinfo &);    // fill in redundant read part of library info
	//void fixBamRedundancy(C_libraryinfo & lib1);    // fix Bam redundancy & fill lib info 
	// SAM @RG header info
	//C_readGroupTags RG;
private:
	string setName;
	string fileName;
	unsigned int Npair_00;    
	unsigned int Npair_01;
	unsigned int Npair_0N;
	unsigned int Npair_11;
	unsigned int Npair_1N;
	unsigned int Npair_NN;    
	// debug
	/*
	 unsigned int Npair_ALU;    
	 unsigned int Npair_ALUC;   
	 */
	// 
	C_pairedreads pairs;
}; // end class 


// ReadGroupID to ReadGroupCode
typedef std::map<string, unsigned long int, std::less<string> > C_ReadGroupID2Code;  

// ReadGroupID to set index
typedef std::map<string, int, std::less<string> >  C_ReadGroupID2set;

// ReadGroupCode to set index
typedef std::map<unsigned long int, int, std::less<unsigned long int> >  C_ReadGroupCode2set;

// Paired-read map file container class
class C_pairedfiles {
  friend ostream &operator<<(ostream &, const C_pairedfiles &);
public:
	C_pairedfiles(string &, RunControlParameters &);  // constructor
	void summaryLog(); 
	void printScanResults(); 
	int Nset;
	vector <C_set> set;
	char inputType;
	char checkSetName(string &);
	string infile[2];
	vector <string> setNames;
	C_ReadGroupID2set ReadGroupID2set;
	C_ReadGroupCode2set ReadGroupCode2set;
	C_ReadGroupID2Code ReadGroupID2Code;
private:
	void loadBam( string  &, RunControlParameters &);
	void testMultiMapBam( string  &, RunControlParameters &);
	//void loadBamSortedByName( string  &, RunControlParameters &);
	//void loadMosaik( string  &, RunControlParameters &);
	void loadSpanner( string  & , RunControlParameters &);
	void loadSpannerDirectory( string  & , RunControlParameters &);
	void loadBamDirectory( string  & ,  RunControlParameters &);
	bool checkSpannerFiles(string &);
	bool checkSpannerDirectory(string &);
	bool checkBamDirectory(string &);
	bool checkBamFile(string &);
	int makeSetsFromLibs(string &, RunControlParameters &, C_anchorinfo &);
	int makeOneSetFromLibs(string &, RunControlParameters &, C_anchorinfo &);
	int  strnum_cmp(const string & a0, const string & b0);
	bool  nextBamAlignmentPair( BamReader & ar1, BamReader & ar2, C_pairedread & p1); 
	bool  nextBamAlignmentZA( BamReader & ar1, C_pairedread & p1); 
	bool  nextBamAlignmentPairSortedByName( BamReader & ar1, C_pairedread & p1);
	bool  nextBamAlignmentJump( BamReader & ar1, BamReader & ar2, C_pairedread & pr1); 
	bool  nextBamAlignmentPairSpecial( BamReader & ar1, C_pairedread & p1);
	int  BamZA2PairedRead(BamAlignment & ba1, C_pairedread & pr1); 
	int  BamBam2PairedRead(BamAlignment & ba1, BamAlignment & ba2, C_pairedread & pr1); 
	int  BamSpecial2PairedRead(BamAlignment & ba1, BamAlignment & ba2, C_pairedread & pr1); 
	//C_pairedread   Mosaik2pair(  Mosaik::AlignedRead& mr);
	//bool  nextBamAlignmentPair( BamMultiReader & ar1, BamMultiReader & ar2, Mosaik::AlignedRead & mr1); 
	//bool  nextBamAlignmentPair( BamReader & ar1, BamReader & ar2, Mosaik::AlignedRead & mr1); 
	//bool  nextBamAlignmentScan( BamMultiReader & ar1, Mosaik::AlignedRead & mr1); 
	//bool  nextBamAlignmentScan( BamReader & ar1, Mosaik::AlignedRead & mr1); 
	//bool  nextBamAlignmentZA( BamReader & ar1, Mosaik::AlignedRead & mr1); 
	//bool  nextBamAlignmentRead( BamMultiReader & ar1, Mosaik::AlignedRead & mr1, bool checkName); 
	//bool  nextBamAlignmentRead( BamReader & ar1, Mosaik::AlignedRead & mr1, bool checkName); 
	//bool  nextBamAlignmentPairSortedByName( BamMultiReader & ar1, Mosaik::AlignedRead & mr1);
	//bool  nextBamAlignmentPairSortedByName( BamReader & ar1, Mosaik::AlignedRead & mr1);
	//int  Bam2MosaikPair( BamAlignment & ba1,BamAlignment & ba2, Mosaik::Alignment & ma1, Mosaik::Alignment & ma2); 
	//int  Bam2MosaikScan( BamAlignment & ba1, Mosaik::Alignment & ma1, Mosaik::Alignment & ma2); 
	//int  Bam2MosaikRead( BamAlignment & ba1, Mosaik::Alignment & ma1); 
	//int  Bam2Mosaik( BamAlignment & ba1, Mosaik::Alignment & ma1, Mosaik::Alignment & ma2); 
	//int  BamZA2Mosaik(BamAlignment & ba1, Mosaik::Alignment & ma1, Mosaik::Alignment & ma2); 
  bool  parseBamReadName(string &, int &, string &);
	int  parseBamTagData(string & tagData, string & tag);
	string  parseBamTagDataString(string & tagData, string & tag);
	string  parseBamHeader(string & samheader, const string & tag) const;
	void inputcheck(string &,  RunControlParameters &);
	int BamCigarData2Len(vector<CigarOp>, bool);
	int BamCigarData2mm(vector<CigarOp>);
	C_anchorinfo anchors;
	C_libraries libraries;
	C_headers headers;
	vector <string> SpannerFileNames;
	vector <string> BamFileNames;
	vector <short int> Shots2Mate;
	int Nbad;
	int NbadPos;
	time_t tprev, tnow;
	float elapsedtime();
	int  MateMode;  // 0: FR (illumina RP short); 1: R1R2/F2F2 (454);  2: F1F2/R2R1 (SOLiD)
	char SpannerMode;
	int BamZ;
	int Qmin;
	// Spanner modes 
	//bool scan, build, detect, genotype, mei, multi;
	int MaxFragments;
	map<string, bool, less<string> > doneFrag;  // keep names of reads already parsed to prevent double counting bam records
	
}; // end class 



// chromosome length and name arrays for AB files, which do not include 
// information about the reference anchor 
const unsigned long ABCHROMOSOMELENGTH[25]={ 
 247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 
 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915,  88827254, 
  78774742,  76117153,  63811651,  62435964,  46944323,  49691432, 154913754,  57772954, 16571};
//---------------------------------------------------------------------------  
//  updated after email from Yongmin - AB used build36.1
//---------------------------------------------------------------------------  
/*
  250781858,246421880,202351854,194005536,183441550,173341421,161090302
 ,148364467,142277156,137308662,136373133,134240242,115773594,107888137
 ,101772329, 90096215, 79900096, 77204541, 64723247, 63327907, 47614957
 , 50401310,157126808, 58598282,    16808};
*/
const string ABCHROMOSOMENAME[25]={"1","2","3","4","5","6","7","8","9","10",
  "11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"};
  
// add a prototype for the Ariya Hidayat's FastLZ library
extern "C" {
  int fastlz_compress(const void* input, int length, void* output);	
  int fastlz_decompress(const void* input, int length, void* output, int maxout);
}

const int DPCOMPRESS=1000000;

union light_t {
    long long i64;
    int i32[2];
    short i16[4];
    char i8[8];
  };

#endif
