/*
 *  RunControlParameterFile.cpp
 *  cnv
 *
 *  Created by Chip Stewart on 10/24/07.
 *  Copyright 2007 Boston College. All rights reserved.
 *
 */
#include "RunControlParameterFile.h"

// list of stuff to trim at ends of parameter strings 
string SPACES=" \t\r\n\"";

// string trim utility functions
inline string trim_right (const string & s, const string & t = SPACES)
{ 
  string d (s); 
  string::size_type i (d.find_last_not_of (t));
  if (i == string::npos) {
    return "";
  } else {
    return d.erase (d.find_last_not_of (t) + 1) ;  
  }
}  // end of trim_right

inline string trim_left (const string & s, const string & t = SPACES) 
{ 
  string d (s); 
  return d.erase (0, s.find_first_not_of (t)) ; 
}  // end of trim_left

inline string trim (const string & s, const string & t = SPACES)
{ 
  string d (s); 
  return trim_left (trim_right (d, t), t) ; 
}  // end of trim

// constructor with default parameters
RunControlParameters::RunControlParameters() {
   RunControlParameters p1("none");
  *this = p1;
}

// constructor loading from file
RunControlParameters::RunControlParameters( const string & filename)
{
  // defaults
  // Regex split by "/"
  string patternSplit1("^(\\S+)/(\\S+)");
  patternSplit=patternSplit1;
  // defaults
  setFragmentLengthWindow(99.9);   // anything outside of window considered for SV
  setFragmentLengthLo(10.0);       // anything of LF less than this is odd
  setFragmentLengthHi(10000.0);    // anything of LF greater than this is also odd
  setFragmentLength(200.0);        // nominal LF - for breakpoint determination 
  setMinClustered(2);          // min number of confirming fragments for SV candidate
  setMinClusteredRetro(2);     // min number of confirming fragments for Mobi candidate
  setClusteringLength(100);    // spatial scale (bp units) for clustering
  setClusteringLengthRetro(60);// spatial scale (bp units) for clustering un-paired reads for insertions
  setDepthLo(0.5);             // region with DofC/nominal less than this considered as deletion
  setDepthHi(1.5);             // region with DofC/nominal greater than this considered as duplication
  setMinLength(1);             // mimum event length 
  setParameterFile("none");
  setReadPairSenseConfig(0);
  setDoInsertions(false);      // turn off insertion detector until it works better
  setPrintTextOut(false);      // Print Text Outout 
  setDoScanOnly(false);        // limit processing to scan -only
  setSpannerMode(1);           // limit processing parse aligment files to build span files
  setReadSetRegex(".");
  setReadIDRegex("^(\\S+)");
  setReadNumberRegex("./(\\d)$");
  setReadEnds("./.");
  setReferenceFastaFile("");
  setAllowContigsRegex(".");
  setMaxFragments(INT_MAX);
  setOutputDir(".");
  setPrefix("");
  setAnchorfile("");
  setWriteDepth(false);
  setMaxMisMatches(1000);      
  setFractionMaxMisMatches(float(1.0));      
  setMinReadLength(1);   
  setRDreferenceFile("");
  setDepthBinSize(1000);  
  setRDminMedian(40);
  setRDminExp(0.8);
  setRDmaxGap(100);
  setRDallow(5);
  setRDminbin(5);
  setDoMasking(false);        
  setQmin(20);
  setCNslosh(0);
  setDupRemove(1);
  //setDoRepeatCheck(true);
  string mobs[] = {"moblist_ALU","moblist_L1","moblist_SVA","moblist_ERV"};
  vector<string> mobv(mobs, mobs + 4);
  setMobileElements(mobv);
  setMobiMaskFile("");
	setBamZA(0);
  
  // Regex Fragment Length Window 
  spatternFLWIN="FragmentLengthWindow";
  // Regex Fragment Length Low 
  spatternFLLO="FragmentLengthLo";
  // Regex Fragment Length High 
  spatternFLHI="FragmentLengthHi";
  // Regex Depth Low 
  spatternDPLO="DepthCoverageRatioLo";
  // Regex Fragment Length High 
  spatternDPHI="DepthCoverageRatioHi";
  // Regex Fragment Length Nominal
  spatternFL="FragmentLengthNominal";
  // Regex Min Clustered size
  spatternMinClus="MinClustered";
  // Regex Min Clustered size for insertions
  spatternMinClusIns="MinClusteredRetroertions";
  // Regex Clustering Length
  spatternClusLen="ClusteringLength";
  // Regex Clustering Length for Insertions
  spatternClusLenIns="ClusteringLengthRetroertsions";
  // Regex Minimum event length 
  spatternMinLength="MinEventLength";
  // Regex Read Set Label in Read Name
  spatternReadSet="ReadSetRegex"; 
  // Regex Read ID in Read Name
  spatternReadID="ReadIDRegex";
  // Regex Read Number in Read Name
  spatternReadNumber="ReadNumberRegex";
  // Regex Read End swap string
  spatternReadEnds="ReadEndSwap";
  // Regex Paired Read default orientiation (short solexa false, long true) 
  spatternReadPairSenseConfig="ReadPairSenseConfig";
  // Regex do insertion detection
  spatternDoInsertion="DoInsertionDetection";
  // ReferenceFastaFile fasta file full path name
  spatternReferenceFastaFile="ReferenceFastaFile";
  // Output area
  spatternOutputDir="OutputDir";
  // Output file prefix
  spatternPrefix="Prefix";
  // specify regexp for allowed contig names
  spatternAllowContigsRegex="AllowContigsRegex";
  // specify maximum number of fragments to process
  spatternMaxFragments="MaxFragments";
  // Regex do PrintTextOut
  spatternPrintTextOut="DoPrintTextOut";
  // Regex do scan only
  spatternDoScanOnly="DoScanOnly";
  // Regex anchorfile
  spatternAnchorfile="Anchorfile";
  // Regex do scan only
  spatternWriteDepth="WriteDepthOfCoverageFiles";
  // Regex maximum number of mismatches 
  spatternMaxMisMatches="MaximumNumberOfMismatches";  
  // Regex fractional maximum number of mismatches 
  spatternFractionMaxMisMatches="FractionalMaximumNumberOfMismatches";  
  // Regex minimum read length 
  spatternMinReadLength="MinimumReadLength";  
  // Regex RDreference filename
  spatternRDreferenceFile="RDreferenceFile"; 
  // Regex Depth bin size
  spatternDepthBinSize="DepthBinSize";    
  // Regex RDminMedian  minimum median read count per bin 
  spatternRDminMedian="RDminMedian";    
  // Regex RDminExp  minimum expected read count per bin 
  spatternRDminExp="RDminExp";    
  // Regex DoMasking
  spatternDoMasking="DoMasking";    
  // Regex Maximum gap in contiguous bins for 1 CNV
  spatternRDmaxGap="RDmaxGap";                       
  // Regex Number of allowed CNV per chromsome w/o chisq penalty
  spatternRDallow="RDallow";                         
  // Regex Minium number of coverage bins in a CNV
  spatternRDminbin="RDminbin";
  // Regex Minium Q map value
  spatternQmin="Qmin";
  // Regex PE CNV slosh 
  spatternCNslosh="CNslosh";
  // Regex duplicate removal
  spatternDupRemove="DupRemove";
  // Regex do repeat 
  spatternDoRepeatCheck="DoRepeatCheck";
  // Regex Mobile elements
  spatternDoRepeatCheck="MobileElements";
  // Regex Mobile Mask
  spatternMobiMaskFile="MobiMaskFile";
  // Regex Bam ZA 
	spatternBamZA="BamZA";

  // list of stuff to trim at ends of parameter strings 
  SPACES=" \t\r\n\"";  
  // regex stuff to parse parameter file 
  string patternFLWIN("^"+spatternFLWIN+"=(\\S+)");
  string patternFLLO("^"+spatternFLLO+"=(\\S+)");
  string patternFLHI("^"+spatternFLHI+"=(\\S+)");
  string patternDPLO("^"+spatternDPLO+"=(\\S+)");
  string patternDPHI("^"+spatternDPHI+"=(\\S+)");
  string patternFL("^"+spatternFL+"=(\\S+)");
  string patternMinClus("^"+spatternMinClus+"=(\\d+)");
  string patternMinClusIns("^"+spatternMinClusIns+"=(\\d+)");
  string patternClusLen("^"+spatternClusLen +"=(\\d+)");
  string patternClusLenIns("^"+spatternClusLenIns+"=(\\d+)");
  string patternReadSet("^"+spatternReadSet+"=""(\\S+)"""); 
  string patternReadID("^"+spatternReadID+"=""(\\S+)""");
  string patternReadNumber("^"+spatternReadNumber+"=""(\\S+)""");
  string patternReadEnds("^"+spatternReadEnds+"=""(\\S+)""");
  string patternReadPairSenseConfig("^"+spatternReadPairSenseConfig+"=""(\\d+)""");
  string patternDoInsertion("^"+spatternDoInsertion+"=""(\\S+)""");
  string patternReferenceFastaFile("^"+spatternReferenceFastaFile+"=""(\\S+)""");
  string patternOutputDir("^"+spatternOutputDir+"=""(\\S+)""");
  string patternPrefix("^"+spatternPrefix+"=""(\\S+)""");
  string patternAllowContigsRegex("^"+spatternAllowContigsRegex+"=(\\S+)");
  string patternMaxFragments("^"+spatternMaxFragments+"=(\\S+)");
  string patternPrintTextOut("^"+spatternPrintTextOut+"=(\\S+)");
  string patternDoScanOnly("^"+spatternDoScanOnly+"=(\\S+)");
  string patternAnchorfile("^"+spatternAnchorfile+"=""(\\S+)""");
  string patternWriteDepth("^"+spatternWriteDepth+"=(\\S+)");
  string patternMinLength("^"+spatternMinLength+"=(\\d+)");
  string patternMaxMisMatches("^"+spatternMaxMisMatches+"=(\\d+)");
  string patternFractionMaxMisMatches("^"+spatternFractionMaxMisMatches+"=([-+]?[0-9]*\\.?[0-9]*)");
  string patternMinReadLength("^"+spatternMinReadLength+"=(\\d+)");
  string patternRDreferenceFile("^"+spatternRDreferenceFile+"=(\\S+)");
  string patternDepthBinSize("^"+spatternDepthBinSize+"=(\\d+)");
  string patternRDminMedian("^"+spatternRDminMedian+"=(\\d+)");
  string patternRDminExp("^"+spatternRDminExp+"=(\\d+)");
  string patternDoMasking("^"+spatternDoMasking+"=(\\d+)");
  string patternRDmaxGap("^"+spatternRDmaxGap+"=(\\d+)");
  string patternRDallow("^"+spatternRDallow+"=(\\d+)");
  string patternRDminbin("^"+spatternRDminbin+"=(\\d+)");
  string patternQmin("^"+spatternQmin+"=(\\d+)");
  string patternCNslosh("^"+spatternCNslosh+"=(\\d+)");
  string patternDupRemove("^"+spatternDupRemove+"=(\\d+)");
  string patternDoRepeatCheck("^"+spatternDoRepeatCheck+"=(\\S+)");
  string patternMobileElements("^"+spatternMobileElements+"=(\\S+)");
  string patternMobiMaskFile("^"+spatternMobiMaskFile+"=(\\S+)");
  string patternBamZA("^"+spatternBamZA+"=(\\d+)");

  //
  if (filename=="none") {
    return;
  }
  // open input parameter file
  ifstream par1(this->ParameterFile.c_str(), ios::in);
  
  // check
  if (!par1) {
    cerr << "Unable to open RCP file: " << filename << endl;
    exit(131);
  }
  
  // define RCP input line
  string line;
  
  // read through RCP file
  while (getline(par1, line)) {
    
    //--------------------------------------------------------------------------
    // load parameters from specified rcp file
    //--------------------------------------------------------------------------
    if (RE2::FullMatch(line.c_str(),patternFLWIN.c_str(),&match)) {
      setFragmentLengthWindow(string2Double(match));
    } else if (RE2::FullMatch(line.c_str(),patternFLLO.c_str(),&match) ) {
      setFragmentLengthLo(string2Double(match));
    } else if (RE2::FullMatch(line.c_str(),patternFLHI.c_str(),&match) ) {
      setFragmentLengthHi(string2Double(match));
    } else if (RE2::FullMatch(line.c_str(),patternDPLO.c_str(),&match) ) {
      setDepthLo(string2Double(match));
    } else if (RE2::FullMatch(line.c_str(),patternDPHI.c_str(),&match) ) {
      setDepthHi(string2Double(match));
    } else if (RE2::FullMatch(line.c_str(),patternFL.c_str(),&match) ) {
      setFragmentLength(string2Double(match));
    } else if (RE2::FullMatch(line.c_str(),patternMinClus.c_str(),&match) ) {
      setMinClustered(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternMinClusIns.c_str(),&match) ) {
      setMinClusteredRetro(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternClusLen.c_str(),&match) ) {
      setClusteringLength(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternClusLenIns.c_str(),&match) ) {
      setClusteringLengthRetro(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternReadSet.c_str(),&match) ) {
      string s = trim(match);
      setReadSetRegex(s);
    } else if (RE2::FullMatch(line.c_str(),patternReadID.c_str(),&match) )  {
      string s = trim(match);
      setReadIDRegex(s);
    } else if (RE2::FullMatch(line.c_str(),patternReadNumber.c_str(),&match) )  {
      string s = trim(match);
      setReadNumberRegex(s);
    } else if (RE2::FullMatch(line.c_str(),patternReadEnds.c_str(),&match) )  {
      string s = trim(match);
      setReadEnds(s);      
    } else if (RE2::FullMatch(line.c_str(),patternReferenceFastaFile.c_str(),&match) ) {
      string s = trim(match);
      setReferenceFastaFile(s);
    } else if (RE2::FullMatch(line.c_str(),patternOutputDir.c_str(),&match) ) {
      string s = trim(match);
      setOutputDir(s);
    } else if (RE2::FullMatch(line.c_str(),patternPrefix.c_str(),&match) )  {
      string s = trim(match);
      setPrefix(s);
    } else if (RE2::FullMatch(line.c_str(),patternAllowContigsRegex.c_str(),&match) ) {
      string s = trim(match);
      setAllowContigsRegex(s);      
    } else if (RE2::FullMatch(line.c_str(),patternMaxFragments.c_str(),&match) ) {
      setMaxFragments(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternReadPairSenseConfig.c_str(),&match) )  {
      setReadPairSenseConfig(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternDoInsertion.c_str(),&match) ) {
      string s = trim(match);
      bool ins1=toupper(s.at(0))=='T';
      setDoInsertions(ins1);         
    } else if (RE2::FullMatch(line.c_str(),patternPrintTextOut.c_str(),&match) ) {
      string s = trim(match);
      bool ins1=toupper(s.at(0))=='T';
      setPrintTextOut(ins1);
    } else if (RE2::FullMatch(line.c_str(),patternDoScanOnly.c_str(),&match) ) {
      string s = trim(match);
      bool ins1=toupper(s.at(0))=='T';
      setDoScanOnly(ins1);
    } else if (RE2::FullMatch(line.c_str(),patternAnchorfile.c_str(),&match) ) {
      string s = trim(match);
      setAnchorfile(s);
    } else if (RE2::FullMatch(line.c_str(),patternMinLength.c_str(),&match) ) {
      setMinLength(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternMaxMisMatches.c_str(),&match) ) {
      setMaxMisMatches(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternFractionMaxMisMatches.c_str(),&match) )  {
      setFractionMaxMisMatches((float)string2Double(match));
    } else if (RE2::FullMatch(line.c_str(),patternMinReadLength.c_str(),&match) ) {
      setMinReadLength(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternRDreferenceFile.c_str(),&match) )  {
      setRDreferenceFile(match);      
    } else if (RE2::FullMatch(line.c_str(),patternDepthBinSize.c_str(),&match) ) {
      setDepthBinSize(string2Int(match));      
    } else if (RE2::FullMatch(line.c_str(),patternWriteDepth.c_str(),&match) ) {
      string s = trim(match);
      bool ins1=toupper(s.at(0))=='T';
      setWriteDepth(ins1);  
    } else if (RE2::FullMatch(line.c_str(),patternRDminMedian.c_str(),&match) ) {
      setRDminMedian(string2Int(match));      
    } else if (RE2::FullMatch(line.c_str(),patternRDminExp.c_str(),&match) ) {
      setRDminExp(string2Int(match));      
    } else if (RE2::FullMatch(line.c_str(),patternDoMasking.c_str(),&match) )  {
      string s = trim(match);
      bool m1=toupper(s.at(0))=='T';
      setDoMasking(m1);         
    } else if (RE2::FullMatch(line.c_str(),patternRDminbin.c_str(),&match) ) {
      setRDminbin(string2Int(match));      
    } else if (RE2::FullMatch(line.c_str(),patternRDmaxGap.c_str(),&match) )  {
      setRDmaxGap(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternQmin.c_str(),&match) ) {
      setQmin(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternCNslosh.c_str(),&match) )  {
      setCNslosh(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternDupRemove.c_str(),&match) ) {
      setDupRemove(string2Int(match));
    } else if (RE2::FullMatch(line.c_str(),patternRDallow.c_str(),&match) ) {
      setRDallow(string2Int(match));            
    } /* else if (RE2::FullMatch(line.c_str(),patternDoRepeatCheck.c_str(),&match) ) {
      string s = trim(match);
      bool m1=toupper(s.at(0))=='T';
      setDoRepeatCheck(m1);         
    } */
	    else if (RE2::FullMatch(line.c_str(),patternMobileElements.c_str(),&match) ) {
      string s = trim(match);
      vector<string> me;
      //int Nme = 
      split(me,s,",");
      setMobileElements(me);
    } else if (RE2::FullMatch(line.c_str(),patternMobiMaskFile.c_str(),&match) ) {
      string s = trim(match);
      setMobiMaskFile(s);
 		} else if (RE2::FullMatch(line.c_str(),patternBamZA.c_str(),&match) ) {
      setBamZA(string2Int(match));
    }
  }
} 

RunControlParameters& RunControlParameters::operator=(const RunControlParameters &rhs)
{
   ParameterFile=rhs.ParameterFile;                   // file name
   FragmentLengthWindow=rhs.FragmentLengthWindow;     // fragment length window in % 
   FragmentLengthLo=rhs.FragmentLengthLo;             // low fragment length  
   FragmentLengthHi=rhs.FragmentLengthHi;             // high fragment length
   FragmentLength=rhs.FragmentLength;                 // nominal fragment length
   MinClustered=rhs.MinClustered;                     // minumum reads in cluster
   ClusteringLength=rhs.ClusteringLength;             // return cluster length parameter
   DepthLo=rhs.DepthLo;                               // threshold low depth of coverage to consider for CNV
   DepthHi=rhs.DepthHi;                               // threshold high depth of coverage to consider for CNV
   dbg=rhs.dbg;                                       // debug level 
   MaxFragments=rhs.MaxFragments;                     // stop loading alignments after this many reads from each end
   ReadSetRegex=rhs.ReadSetRegex;                     // Pair Set Name Regex string
   ReadIDRegex=rhs.ReadIDRegex;                       // Pair ID Regex string
   ReadNumberRegex=rhs.ReadNumberRegex;               // Pair Number Regex string 
   ReadEnd1=rhs.ReadEnd1;                             // Pair End 1 string 
   ReadEnd2=rhs.ReadEnd2;                             // Pair End 2 string 
   AllowContigsRegex=rhs.AllowContigsRegex;           // only process contigs with names consistent with this regex
   ReadPairSenseConfig=rhs.ReadPairSenseConfig;       // True if DNA sense if flipped within a paired-end
   ReferenceFastaFile=rhs.ReferenceFastaFile;         // Fasta reference file 
   doInsertions=rhs.doInsertions;                     // switch for Insertion detector
   ClusteringLengthRetro=rhs.ClusteringLengthRetro;   // return cluster length parameter for insertions
   MinClusteredRetro=rhs.MinClusteredRetro;           // minumum reads in cluster for insertions
   doPrintTextOut=rhs.doPrintTextOut;                 // switch for text output
   doScanOnly=rhs.doScanOnly;                         // switch for scanning input files for fragment length only
   OutputDir=rhs.OutputDir;                           // name of output directory
   Prefix=rhs.Prefix;                                 // prefix in name of output files  
   Anchorfile=rhs.Anchorfile;                         // Anchor info file
   writeDepth=rhs.writeDepth;   
   SpannerMode=rhs.SpannerMode;                 
   HistoGroups=rhs.HistoGroups;     
   MinLength=rhs.MinLength;                           // minumum event length 
   MaxMisMatches=rhs.MaxMisMatches;                   // maximum number of mismatches 
   FractionMaxMisMatches=rhs.FractionMaxMisMatches;   // fractional maximum number of mismatches 
   MinReadLength=rhs.MinReadLength;                   // minimum read length 
   RDreferenceFile=rhs.RDreferenceFile;
   DepthBinSize=rhs.DepthBinSize;
   RDminMedian=rhs.RDminMedian;
   RDminExp=rhs.RDminExp;
   doMasking=rhs.doMasking;                 
   RDmaxGap=rhs.RDmaxGap;
   RDminbin=rhs.RDminbin;
   Qmin=rhs.Qmin;
   CNslosh=rhs.CNslosh;
   DupRemove=rhs.DupRemove;
   RDallow=rhs.RDallow;
   doRepeatCheck=rhs.doRepeatCheck;
   MobileElements=rhs.MobileElements;
   MobiMaskFile=rhs.MobiMaskFile;
	 BamZA=rhs.BamZA;
   return *this;
}

// return parameter file name
string RunControlParameters::getParameterFile() const
{
	return ParameterFile;
}
// return low fragment length
double RunControlParameters::getFragmentLengthWindow() const
{
	return FragmentLengthWindow;
}
// return low fragment length
double RunControlParameters::getFragmentLengthLo() const
{
	return FragmentLengthLo;
}
// return high fragment length
double RunControlParameters::getFragmentLengthHi() const
{
	return FragmentLengthHi;
}
// nominal fragment length
double RunControlParameters::getFragmentLength() const
{
	return FragmentLength;
}
// return minumum reads in cluster
int RunControlParameters::getMinClustered() const
{
	return MinClustered;
}
// return minumum reads in cluster for insertions
int RunControlParameters::getMinClusteredRetro() const
{
	return MinClusteredRetro;
}
// return cluster length parameter   
int RunControlParameters::getClusteringLength() const
{
	return ClusteringLength;
}
// return cluster length parameter   
int RunControlParameters::getClusteringLengthRetro() const
{
	return ClusteringLengthRetro;
}
// return low depth ratio threshold
double RunControlParameters::getDepthLo() const
{
	return DepthLo;
}
// return high depth threshold
double RunControlParameters::getDepthHi() const
{
	return DepthHi;
}
// debug level 
int RunControlParameters::getDbg() const
{
	return dbg;
}
// regex expression for Read Set string
string RunControlParameters::getReadSetRegex() const
{
	return ReadSetRegex;
}
// regex expression for Read ID string
string RunControlParameters::getReadIDRegex() const
{
	return ReadIDRegex;
}
// get regex expression for Read Number string
string RunControlParameters::getReadNumberRegex() const
{
	return ReadNumberRegex;
}

// set file name
void RunControlParameters::setParameterFile(const string & filename) 
{
  ParameterFile=filename;
} 
// set Fragment Length Window
void RunControlParameters::setFragmentLengthWindow(const double FLW)
{
  FragmentLengthWindow=FLW;
} 
// set Fragment Length Low
void RunControlParameters::setFragmentLengthLo(const double FLL)
{
  FragmentLengthLo=FLL;
} 
// set Fragment Length High
void RunControlParameters::setFragmentLengthHi(const double FLH)
{
  FragmentLengthHi=FLH;
} 
// set nominal Fragment Length 
void RunControlParameters::setFragmentLength(const double FL)
{
  FragmentLength=FL;
} 
// set minumum reads in cluster
void RunControlParameters::setMinClustered(const int MRC)
{
  MinClustered=MRC;
} 
// set minumum reads in cluster
void RunControlParameters::setMinClusteredRetro(const int MRC)
{
  MinClusteredRetro=MRC;
} 
// set cluster length parameter 
void RunControlParameters::setClusteringLength(const int CL)
{
  ClusteringLength=CL;
} 
// set cluster length parameter for insertions
void RunControlParameters::setClusteringLengthRetro(const int CL)
{
  ClusteringLengthRetro=CL;
} 
// set Depth Low
void RunControlParameters::setDepthLo(const double DPL)
{
  DepthLo=DPL;
} 
// set Depth High
void RunControlParameters::setDepthHi(const double DPH)
{
  DepthHi=DPH;
} 

// set debug level 
void RunControlParameters::setDbg(const int i)
{
  dbg=i;
} 
// set regex expression for Read Set string
void RunControlParameters::setReadSetRegex(const string & s) 
{
  ReadSetRegex=s;
} 
// set regex expression for Read ID string
void RunControlParameters::setReadIDRegex(const string & s) 
{
  ReadIDRegex=s;
} 
// set regex expression for Read Number string
void RunControlParameters::setReadNumberRegex(const string & s) 
{
  ReadNumberRegex=s;
} 
// get expression for Read Number 1 swap string
string RunControlParameters::getReadEnd1() const
{
  return ReadEnd1;
}
// get expression for Read Number 2 swap string
string RunControlParameters::getReadEnd2() const
{
  return ReadEnd2;
}
// set expressions for Read swap strings
void RunControlParameters::setReadEnds(const string & s)
{
	string match2;
  if (RE2::FullMatch(s.c_str(),patternSplit.c_str(),&match,&match2) ) {
    ReadEnd1=match;
    ReadEnd2=match2;
  }
}	
// ReferenceFastaFile Fastafile name
string RunControlParameters::getReferenceFastaFile() const
{
  return ReferenceFastaFile;
}
// set ReferenceFastaFile Fastafile name
void  RunControlParameters::setReferenceFastaFile(const string & s)
{
  ReferenceFastaFile=s;
}
 // True if DNA sense if flipped within a paired-end
int RunControlParameters::getReadPairSenseConfig() const
{
   return ReadPairSenseConfig;
}
void RunControlParameters::setReadPairSenseConfig(const int c1)
{
  ReadPairSenseConfig=c1;
}  
// run Insertion detector
bool RunControlParameters::getDoInsertions() const 
{
  return doInsertions;
}  
void RunControlParameters::setDoInsertions(const bool ins1)
{
   doInsertions = ins1;
}
// set regexp for allowed contigs 
void RunControlParameters::setAllowContigsRegex(const string & s) 
{
  AllowContigsRegex=s;
} 
// set regexp for allowed contigs 
string  RunControlParameters::getAllowContigsRegex()  const
{
  return AllowContigsRegex;
} 
// set Max number of fragments
void RunControlParameters::setMaxFragments(const int i)
{
  MaxFragments=i;
} 
// get Max number of fragments
int RunControlParameters::getMaxFragments() const
{
  return MaxFragments;
} 
// scan only
bool RunControlParameters::getDoScanOnly() const 
{
  return doScanOnly;
}  
void RunControlParameters::setDoScanOnly(const bool ins1)
{
   doScanOnly = ins1;
}

// scan only
bool RunControlParameters::getPrintTextOut() const 
{
  return doPrintTextOut;
}  
void RunControlParameters::setPrintTextOut(const bool pto1)
{
   doPrintTextOut = pto1;
}

// output directory
string RunControlParameters::getOutputDir() const
{
  return OutputDir;
}
void  RunControlParameters::setOutputDir(const string & s)
{
  OutputDir=s;
}
// output prefix
string RunControlParameters::getPrefix() const
{
  return Prefix;
}
void  RunControlParameters::setPrefix(const string & s)
{
  Prefix=s;
}

C_HistoGroups RunControlParameters::getHistoGroups() const
{
  return HistoGroups;
}

void RunControlParameters::setHistoGroups(const C_HistoGroups & h1) 
{
    HistoGroups = h1;
}

 
// Anchor info file
string RunControlParameters::getAnchorfile() const
{
  return Anchorfile;
}
void  RunControlParameters::setAnchorfile(const string & s)
{
  Anchorfile=s;
}

// write depth files
bool RunControlParameters::getWriteDepth() const 
{
  return writeDepth;
}  
void RunControlParameters::setWriteDepth(const bool ins1)
{
   writeDepth = ins1;
}

// set minumum  event length
void RunControlParameters::setMinLength(const int L)
{
  MinLength=L;
} 
int RunControlParameters::getMinLength() const
{
	return MinLength;
}

// return depth bin size
int RunControlParameters::getDepthBinSize() const
{
	return DepthBinSize;
}

// set minumum read length
void RunControlParameters::setDepthBinSize(const int b)
{  
 DepthBinSize=b; 
} 

// return minumum read length 
int RunControlParameters::getMinReadLength() const
{
	return MinReadLength;
}
void RunControlParameters::setMinReadLength(const int M)
{  
 MinReadLength=M; 
} 

// set maximum number of mismatches 
void RunControlParameters::setMaxMisMatches(const int M)
{  
 MaxMisMatches=M; 
} 
// return maximum number of mismatches
int RunControlParameters::getMaxMisMatches() const
{
	return MaxMisMatches;
}
// set fractional maximum number of mismatches
void RunControlParameters::setFractionMaxMisMatches(const float F)
{  
  FractionMaxMisMatches=F; 
} 
// return fractional maximum number of mismatches
float RunControlParameters::getFractionMaxMisMatches() const
{
	return FractionMaxMisMatches;
}
// set RDreference Filename
void RunControlParameters::setRDreferenceFile(const string & S)
{  
  RDreferenceFile=S; 
} 
// return RDreference Filename
string RunControlParameters::getRDreferenceFile() const
{
	return RDreferenceFile;
}

// return minumum median read counts per bin (rebin)
int RunControlParameters::getRDminMedian() const
{
	return RDminMedian;
}
void RunControlParameters::setRDminMedian(const int M)
{  
 RDminMedian=M; 
} 
// return minumum expected read counts in each bin (clean)
float RunControlParameters::getRDminExp() const
{
	return RDminExp;
}
void RunControlParameters::setRDminExp(const float M)
{  
 RDminExp=M; 
} 
// run Insertion detector
bool RunControlParameters::getDoMasking() const 
{
  return doMasking;
}  
void RunControlParameters::setDoMasking(const bool m1)
{
   doMasking = m1;
}

int RunControlParameters::getRDminbin() const
{
	return RDminbin;
}
void RunControlParameters::setRDminbin(const int M)
{  
 RDminbin=M; 
} 
int RunControlParameters::getRDmaxGap() const
{
	return RDmaxGap;
}
void RunControlParameters::setRDmaxGap(const int M)
{  
 RDmaxGap=M; 
} 
int RunControlParameters::getRDallow() const
{
	return RDallow;
}
void RunControlParameters::setRDallow(const int M)
{  
 RDallow=M; 
} 


int RunControlParameters::getQmin() const
{
	return Qmin;
}
void RunControlParameters::setQmin(const int M)
{
 Qmin=M;
}

int RunControlParameters::getCNslosh() const
{
	return CNslosh;
}
void RunControlParameters::setCNslosh(const int M)
{
 CNslosh=M;
}

int RunControlParameters::getDupRemove() const
{
	return DupRemove;
}
void RunControlParameters::setDupRemove(const int M)
{
 DupRemove=M;
}

/*
 bool RunControlParameters::getDoRepeatCheck() const
{
	return doRepeatCheck;
}
void RunControlParameters::setDoRepeatCheck(const bool M)
{
 doRepeatCheck=M;
}
*/

// set regexp for Mobile Elements list
void RunControlParameters::setMobileElements(const vector<string> & s) 
{
  MobileElements=s;
} 
vector<string>  RunControlParameters::getMobileElements()  const
{
  return MobileElements;
} 

// set regexp for Mobile Elements list
void RunControlParameters::setMobiMaskFile(const string & s) 
{
  MobiMaskFile=s;
} 
string  RunControlParameters::getMobiMaskFile()  const
{
  return MobiMaskFile;
} 

// set bam ZA tag required
void RunControlParameters::setBamZA(const int i)
{
  BamZA=i;
} 
// get Bam ZA required flag
int RunControlParameters::getBamZA() const
{
  return BamZA;
} 


// SPanner mode
int RunControlParameters::getSpannerMode() const 
{
  return SpannerMode;
}  
void RunControlParameters::setSpannerMode(const int b1)
{
	SpannerMode = b1;
}


/*
void RunControlParameters::setFragmentLengthLimits() 
{
   if (this->fragHist.Ntot==0) {
     return;
   }
   // nominal mapping length stats
  this->FragmentLength = fragHist.mode;  
  double lowf = (1.0-FragmentLengthWindow/100.0)/2.0;  // half tail on low side
  this->FragmentLengthLo = fragHist.p2xTrim(lowf);
  double highf = (1.0-lowf); // half tail on high side     
  this->FragmentLengthHi = fragHist.p2xTrim(highf);
  this->MinLength = 3*int(fragHist.std);
}
*/

ostream &operator<<(ostream &output,  const RunControlParameters & p1)
{
    string line = "//================================================";    
    char s[100];
    int rc;
    time_t temp;
    struct tm *timeptr;
    temp = time(NULL);
    timeptr = localtime(&temp);
    rc = strftime(s,sizeof(s),"%A, %b %d:  %r", timeptr);
    output << line << endl;
    output << "// original file:\t " << p1.getParameterFile() << endl;
    output << s << endl;
    output << line << endl;
    output << "// \tNominal fragment length: " << endl;
    output << p1.spatternFL << "=" << p1.getFragmentLength()  << endl;
    output << "// \tLow fragment length: " << endl;
    output << p1.spatternFLLO << "=" << p1.getFragmentLengthLo()  << endl;
    output << "// \tHigh fragment length: " << endl;
    output << p1.spatternFLHI << "=" << p1.getFragmentLengthHi()  << endl;
    output << "// \tNominal fragment length window fraction: " << endl;
    output << p1.spatternFLWIN << "=" << p1.getFragmentLengthWindow()  << endl;
    output << "// \tMin clustered: " << endl;
    output << p1.spatternMinClus << "=" << p1.getMinClustered()  << endl;
    output << "// \tClustering length: " << endl;
    output << p1.spatternClusLen << "=" << p1.getClusteringLength() << endl;
    output << "//\tLow depth threshold: " << endl;
    output << p1.spatternDPLO << "=" << p1.getDepthLo() << endl;
    output << "//\tHigh depth threshold: " << endl;
    output << p1.spatternDPHI << "=" << p1.getDepthHi() << endl;
    output << "// \tMin length: " << endl;
    output << p1.spatternMinLength << "=" << p1.getMinLength()  << endl;
    output << "// \tMaximum number of mismatches: " << endl;
    output << p1.spatternMaxMisMatches << "=" << p1.getMaxMisMatches()  << endl;
    output << "// \tFractional Maximum number of mismatches: " << endl;
    output << p1.spatternFractionMaxMisMatches << "=" << p1.getFractionMaxMisMatches()  << endl;
    output << "// \tMinimum Read length: " << endl;
    output << p1.spatternMinReadLength << "=" << p1.getMinReadLength()  << endl;
    output << "//\tMaximum number of fragments: " << endl;
    output << p1.spatternMaxFragments << "=" << p1.getMaxFragments()  << endl;
    output << "//\tPrintTextOut: " << endl;
    output << p1.spatternPrintTextOut << "=" << (p1.getPrintTextOut()? "T": "F")  << endl;
    output << "//\tScan-only: " << endl;
    output << p1.spatternDoScanOnly << "=" << (p1.getDoScanOnly()? "T": "F")  << endl;
    output << "//\tDetect Insertion SV events: " << endl;
    output << p1.spatternDoInsertion << "=" << (p1.getDoInsertions()? "T": "F")  << endl; 
    output << "//\tMask file: " << endl;
    output << p1.spatternReferenceFastaFile << "=""" << p1.getReferenceFastaFile() << """" << endl;
    output << "//\tRead set selection regex: " << endl;
    output << p1.spatternReadSet << "=""" << p1.getReadSetRegex() << """" << endl;
    output << "//\tRead ID  regex: " << endl;
    output << p1.spatternReadID << "=""" << p1.getReadIDRegex() << """" << endl;
    output << "//\tRead number  regex: " << endl;
    output << p1.spatternReadNumber << "=""" << p1.getReadNumberRegex() << """" << endl;
    output << "//\tRead-ends swap strings  regex: " << endl;
    output << p1.spatternReadEnds << "=""" << p1.getReadEnd1() <<"/"<< p1.getReadEnd2() << """" << endl;
    output << "//\tAllowed contig names regex: " << endl;
    output << p1.spatternAllowContigsRegex << "=""" << p1.getAllowContigsRegex() << """" << endl;
    output << "//\tOutput area: " << endl;
    output << p1.spatternOutputDir << "=""" << p1.getOutputDir() << """" << endl;
    output << "//\tOutput file prefix: " << endl;
    output << p1.spatternPrefix << "=""" << p1.getPrefix() << """" << endl;
    output << "//\tAnchor file : " << endl;
    output << p1.spatternAnchorfile << "=""" << p1.getAnchorfile() << """" << endl;
    output << "//\tWrite Depth files: " << endl;
    output << p1.spatternWriteDepth << "=" << (p1.getWriteDepth()? "T": "F")  << endl;
    output << "//\tRDreference file: " << endl;
    output << p1.spatternRDreferenceFile << "=""" << p1.getRDreferenceFile()  << """" << endl;
    output << "//\tDepth bin size: " << endl;
    output << p1.spatternDepthBinSize << "=""" << p1.getDepthBinSize()  << """" << endl;
    output << "//\tRD min median read counts per bin: " << endl;
    output << p1.spatternRDminMedian << "=""" << p1.getRDminMedian()  << """" << endl;
    output << "//\tRD min expected read counts per bin: " << endl;
    output << p1.spatternRDminExp << "=""" << p1.getRDminExp()  << """" << endl;
    output << "//\tRepeat Masking: " << endl;
    output << p1.spatternDoMasking << "=" <<  (p1.getDoMasking()? "T": "F")  << endl;
    output << "//\tRD min bin CNV : " << endl;
    output << p1.spatternRDminbin << "=""" << p1.getRDminbin()  << """" << endl;
    output << "//\tRD max gap CNV : " << endl;
    output << p1.spatternRDmaxGap << "=""" << p1.getRDmaxGap()  << """" << endl;
    output << "//\tRD max allow CNV w/o test : " << endl;
    output << p1.spatternRDallow << "=""" << p1.getRDallow()  << """" << endl;
    output << "//\tQ min : " << endl;
    output << p1.spatternQmin << "=""" << p1.getQmin()  << """" << endl;
    output << "//\tCN slosh : " << endl;
    output << p1.spatternCNslosh << "=""" << p1.getCNslosh()  << """" << endl;
    output << "//\tDuplicate fragment removal : " << endl;
    output << p1.spatternDupRemove << "=""" << p1.getDupRemove()  << """" << endl;
    //output << "//\tRepeat mapping location checking: " << endl;
    //output << p1.spatternDoRepeatCheck << "=" <<  (p1.getDoRepeatCheck()? "T": "F")  << endl;  
    output << "//\tMobile Element list: " << endl;
    string me=p1.MobileElements[0];
    for (int i=1; i<int(p1.MobileElements.size());i++) {
       me=","+p1.MobileElements[i];
    }
    output << p1.spatternMobileElements << "=" <<  me  << endl;  
		output << p1.spatternMobiMaskFile << "=" <<  p1.getMobiMaskFile()  << endl;  
	  output << "//\tBam ZA : " << endl;
	  output << p1.spatternBamZA << "=""" << p1.getBamZA ()  << """" << endl;
	
    return output;
}


// print RunControlParameters object
void RunControlParameters::print() const
{
  cout << "\nRun control parameter file: " << getParameterFile() << endl;
  // non-default values only  
  RunControlParameters p1;
  if (p1.getFragmentLength()!=getFragmentLength()) {
    //cout << "\t//Nominal fragment length: " << endl;
    cout << "\t" << spatternFL << "=" << getFragmentLength()  << endl;
  }
  if (p1.getFragmentLengthLo()!=getFragmentLengthLo()) {
    //cout << "\t" << "// \tLow fragment length: " << endl;
    cout << "\t" << spatternFLLO << "=" << getFragmentLengthLo()  << endl;
  }
  if (p1.getFragmentLengthHi()!=getFragmentLengthHi()) {
    //cout << "\t"<< "// \tHigh fragment length: " << endl;
    cout << "\t"<< spatternFLHI << "=" << getFragmentLengthHi()  << endl;
  }
  if (fabs(p1.getFragmentLengthWindow()-getFragmentLengthWindow())>1e-5) {
    //cout << "\t//normal fragment length window fraction: " << endl;
    cout << "\t" << spatternFLWIN << "=" << getFragmentLengthWindow()  << endl;
  }
  if (p1.getMinClustered()!=getMinClustered()) {
    //cout << "\t"<< "// \tMin clustered: " << endl;
    cout << "\t"<< spatternMinClus << "=" << getMinClustered()  << endl;
  }
  if (p1.getClusteringLength() !=getClusteringLength() ) {
    //cout << "\t"<< "// \tClustering length: " << endl;
    cout << "\t"<< spatternClusLen << "=" << getClusteringLength() << endl;
  }
  if (p1.getDepthLo() !=getDepthLo()) {
    //cout << "\t" <<"//\tLow depth threshold: " << endl;
    cout << "\t" << spatternDPLO << "=" << getDepthLo() << endl;
  }
  if (p1.getDepthHi() !=getDepthHi()) {
    //cout << "\t" <<"//\tHigh depth threshold: " << endl;
    cout << "\t" << spatternDPHI << "=" << getDepthHi() << endl;
  }
  if (p1.getMinLength()  !=getMinLength() ) {
    //cout << "\t" <<"// \tMin length: " << endl;
    cout << "\t" << spatternMinLength << "=" << getMinLength()  << endl;
  }
  if (p1.getMaxMisMatches() !=getMaxMisMatches() ) {
    //cout << "\t" <<"// \tMaximum number of mismatches: " << endl;
    cout << "\t" << spatternMaxMisMatches << "=" << getMaxMisMatches()  << endl;
  }
  if (p1.getFractionMaxMisMatches() !=getFractionMaxMisMatches() ) {
    //cout << "\t" <<"// \tFractional Maximum number of mismatches: " << endl;
    cout << "\t" << spatternFractionMaxMisMatches << "=" << getFractionMaxMisMatches()  << endl;
  }
  if (p1.getMinReadLength() !=getMinReadLength()  ) {
    //cout << "\t" <<"// \tMinimum Read length: " << endl;
    cout << "\t" << spatternMinReadLength << "=" << getMinReadLength()  << endl;
  }
  if (p1.getMaxFragments() !=getMaxFragments()  ) {
    //cout << "\t" <<"//\tMaximum number of fragments: " << endl;
    cout << "\t" << spatternMaxFragments << "=" << getMaxFragments()  << endl;
  }
  if (p1.getPrintTextOut() !=getPrintTextOut() ) {
    //cout << "\t" <<"//\tPrintTextOut: " << endl;
    cout << "\t" << spatternPrintTextOut << "=" << (getPrintTextOut()? "T": "F")  << endl;
  }
  if (p1.getDoScanOnly() !=getDoScanOnly() ) {
    //cout << "\t" <<"//\tScan-only: " << endl;
    cout << "\t" << spatternDoScanOnly << "=" << (getDoScanOnly()? "T": "F")  << endl;
  }
  if (p1.getDoInsertions() !=getDoInsertions() ) {
    //cout << "\t" <<"//\tDetect Insertion SV events: " << endl;
    cout << "\t" << spatternDoInsertion << "=" << (getDoInsertions()? "T": "F")  << endl; 
  }
  if (p1.getReferenceFastaFile().compare(getReferenceFastaFile()) ) {
    //cout << "\t" <<"//\tMask file: " << endl;
    cout << "\t" << spatternReferenceFastaFile << "=""" << getReferenceFastaFile() << """" << endl;
  }
  if (p1.getReadSetRegex().compare(getReadSetRegex()) ) {
    //cout << "\t" <<"//\tRead set selection regex: " << endl;
    cout << "\t" <<p1.spatternReadSet << "=""" << getReadSetRegex() << """" << endl;
  }
  if (p1.getReadIDRegex().compare(getReadIDRegex()) ) {
    //cout << "\t" <<"//\tRead ID  regex: " << endl;
    cout << "\t" <<spatternReadID << "=""" << getReadIDRegex() << """" << endl;
  }
  if (p1.getReadNumberRegex().compare(getReadNumberRegex() ) ) {
    //cout << "\t" <<"//\tRead number  regex: " << endl;
    cout << "\t" <<spatternReadNumber << "=""" << getReadNumberRegex() << """" << endl;
  }
  if (p1.getReadEnd1().compare(getReadEnd1() ) ) {
    //cout << "\t" <<"//\tRead-ends swap strings  regex: " << endl;
    cout << "\t" <<p1.spatternReadEnds << "=""" << getReadEnd1() <<"/"<< getReadEnd2() << """" << endl;
  }
  if (p1.getAllowContigsRegex().compare(getAllowContigsRegex() ) ) {
    //cout << "\t" <<"//\tAllowed contig names regex: " << endl;
    cout << "\t" <<spatternAllowContigsRegex << "=""" << getAllowContigsRegex() << """" << endl;
  }
  if (p1.getOutputDir().compare(getOutputDir() ) ) {
    //cout << "\t" <<"//\tOutput area: " << endl;
    cout << "\t" <<spatternOutputDir << "=""" << getOutputDir() << """" << endl;
  }
  if (p1.getPrefix().compare(getPrefix()  ) ) { 
    //cout << "\t" <<"//\tOutput file prefix: " << endl;
    cout << "\t" <<spatternPrefix << "=""" << getPrefix() << """" << endl;
  }
  if (p1.getAnchorfile().compare(getAnchorfile() ) ) { 
    //cout << "\t" <<"//\tAnchor file : " << endl;
    cout << "\t" <<spatternAnchorfile << "=""" << getAnchorfile() << """" << endl;
  }
  if (p1.getWriteDepth()!=getWriteDepth() ) { 
    //cout << "\t" <<"//\tWrite Depth files: " << endl;
    cout << "\t" <<spatternWriteDepth << "=" << (getWriteDepth()? "T": "F")  << endl;
  }
  if (p1.getRDreferenceFile().compare(getRDreferenceFile() ) ) { 
    //cout << "\t" <<"//\tRDreference file: " << endl;
    cout << "\t" <<spatternRDreferenceFile << "=""" << getRDreferenceFile() << """"  << endl;
  }
  if (p1.getDepthBinSize()!=getDepthBinSize()  ) { 
    //cout << "\t" <<"//\tRDreference file: " << endl;
    cout << "\t" <<spatternDepthBinSize << "=" << getDepthBinSize()   << endl;
  }
  if (p1.getRDminMedian()!=getRDminMedian()  ) { 
    cout << "\t" <<spatternRDminMedian << "=" << getRDminMedian()   << endl;
  }
  if (p1.getRDminExp()!=getRDminExp()  ) { 
    cout << "\t" <<spatternRDminExp << "=" << getRDminExp()   << endl;
  }
  if (p1.getDoMasking()!=getDoMasking() ) { 
    cout << "\t" <<spatternDoMasking << "=" << (getDoMasking()? "T": "F")  << endl;
  }
  if (p1.getRDminbin()!=getRDminbin()  ) { 
    cout << "\t" <<spatternRDminbin << "=" << getRDminbin()   << endl;
  }
  if (p1.getRDmaxGap()!=getRDmaxGap()  ) { 
    cout << "\t" <<spatternRDmaxGap << "=" << getRDmaxGap()   << endl;
  }
  if (p1.getRDallow()!=getRDallow()  ) { 
    cout << "\t" <<spatternRDallow << "=" << getRDallow()   << endl;
  }
  if (p1.getQmin()!=getQmin()  ) {
    cout << "\t" <<spatternQmin << "=" << getQmin()   << endl;
  }
  if (p1.getCNslosh()!=getCNslosh()  ) {
    cout << "\t" <<spatternCNslosh << "=" << getCNslosh()   << endl;
  }
  if (p1.getDupRemove()!=getDupRemove()  ) {
    cout << "\t" <<spatternDupRemove << "=" << getDupRemove()   << endl;
  }
  /*
	if (p1.getDoRepeatCheck()!=getDoRepeatCheck()  ) {
    cout << "\t" <<spatternDoRepeatCheck << "=" << getDoRepeatCheck()   << endl;
  }
	*/
  if (p1.getMobileElements()!=getMobileElements() ){ 
    string me=MobileElements[0];
    for (int i=1; i<int(MobileElements.size());i++) {
       me=","+MobileElements[i];
    }
    cout << "\t" << p1.spatternMobileElements << "=" <<  me  << endl;  
  }
  if (p1.getMobiMaskFile().compare(getMobiMaskFile())!=0 ) {
    cout << "\t" <<spatternMobiMaskFile << "=" << getMobiMaskFile()   << endl;
  }
  
	if (p1.getBamZA()!=getBamZA()  ) {
    cout << "\t" <<spatternBamZA << "=" << getBamZA()   << endl;
  }
	
  cout << "\n" << flush;
}
