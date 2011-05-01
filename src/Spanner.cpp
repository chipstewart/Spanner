//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Spanner 
// Parse input files to create structures for detection of Structural Variations 
// (SV) from mapped paired-end reads 
// Copyright 2009 Chip Stewart, Boston College
// All rights reserved
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// standard includes
#include <iostream>
#include <fstream>
#include <string>

// "tclap" commandline parsing library
#include <tclap/CmdLine.h>

// private libraries
#include "LargeFileSupport.h"
#include "RunControlParameterFile.h"
#include "PairedData.h"
#include "SpanDet.h"
#include "DepthCnvDet.h"
#include "SpannerVersion.h"

// uses
using namespace std; 
using namespace TCLAP; 
using namespace __gnu_cxx;

int main (int argc, char *  argv[]) {    //int main (int argc, char * const argv[]) {
  
  //string  V = "3.44";
  string V = SpannerSvnVersion;
  //----------------------------------------------------------------------------
  // hi there
  //----------------------------------------------------------------------------  
  std::cout << "begin Spanner! \t Version " << V << "\n"; //int dbg =1;
  std::cout << "command line:\n"; 
  //----------------------------------------------------------------------------
  // dump command line to std out
  //----------------------------------------------------------------------------
  for(int i = 0; i < argc; i++) {
      cout <<  argv[i] << " ";
  }
  cout << endl;
  //----------------------------------------------------------------------------
  // Create new CmdLine object
  //----------------------------------------------------------------------------
  CmdLine cmd("Command line options for Spanner", ' ', V );    
  //----------------------------------------------------------------------------
  // command line argument parsing
  //----------------------------------------------------------------------------
  // run control parameter file
  ValueArg<string> cmd_rcp("R", "rcp", "name of parameter file", false, "none", "string", cmd);

  // output file area
  ValueArg<string> cmd_out("o", "outdir", "name of SV output area", false, ".", "string", cmd);

  // output file prefix
  ValueArg<string> cmd_prefix("p", "prefix", "output prefix", false, "", "string", cmd);
 
  // fasta file for GC content correction
  // ValueArg<string> cmd_fasta("F", "fasta", "name of reference fasta file", false, "", "string", cmd);
 
  // fragHist file from scan
  ValueArg<string> cmd_fragHist("f", "fragHist", "name of fragment length histogram file", false, "", "string", cmd);
 
  // Max Fragments
  ValueArg<int> cmd_maxfragments("x", "maxfrag", "maximum number of fragments", false, 0, "int", cmd);

  // aligner input file 1
  // ValueArg<string> cmd_in1("1", "in1", "name of first input file", false, "", "string", cmd);

  // aligner input file 2
  // ValueArg<string> cmd_in2("2", "in2", "name of second  input file", false, "", "string", cmd);

  // build input specifier
  ValueArg<string> cmd_input("i", "infile", "specify pre-built input files", false, "", "string", cmd);

  // Allowed Contigs 
  ValueArg<string> cmd_allowcontigs("c", "contigs", "regexp for allowed contigs", false, ".", "string", cmd);

  // Detection Algorithms  (always RP)
  // ValueArg<string> cmd_algorithms("A", "algorithms", "algorthms:RP,RD", false, "RP", "string", cmd);

  // CN slosh on selection criteria for PE detected CNV's
  // ValueArg<int> cmd_cnslosh("S", "cnslosh", "CN slosh for PE CNVs", false, 0, "int", cmd);

 // Duplicate Removal 
  ValueArg<int> cmd_uniquify("u", "uniquify", "duplicate fragment removal", false, 1, "int", cmd);

  // scan input files to extract fragment length, read length, and repeat stats
  SwitchArg cmd_scan("s", "scan", "scan fragment length only", false);
  cmd.add( cmd_scan);
  
  // parse input files to build span data structure/files
  SwitchArg cmd_build("b", "build", "build span files", false);
  cmd.add( cmd_build);
  
  // paired read sense convention (from bam file)
  // SwitchArg cmd_flip("4", "454PE", "454 paired read sense flip relative to Illumina PE", false);
  // cmd.add( cmd_flip);
  
  // write depth of coverage files
  // ValueArg<int> cmd_depth("D", "depth", "make depth of coverage  bin size", false,0,"int",cmd);
  
  // anchor info  file
  ValueArg<string> cmd_Anchorfile("a", "anchorfile", "name of anchor info file", false, "", "string", cmd);

  // anchor info  file
  // ValueArg<string> cmd_RDreference("r", "RDreference", "name of read depth reference file", false, "", "string", cmd);
 
  // ascii text outout
  SwitchArg  cmd_textout("t", "text", "ascii text output", false);
  cmd.add( cmd_textout);

  // debug bits
	int DBGDefault=0;
  ValueArg<int> cmd_dbg("d", "debug", "debug: interval>0, RD<0", false, DBGDefault, "int", cmd);

  // mapping Q min 
	int QminDefault=20;
  ValueArg<int> cmd_qmin("Q", "qmin", "minimum mapping q value", false,QminDefault,"int",cmd);

  // alternate Fragment tail p-value min 
	double FragTailDefault=99.9;
  ValueArg<double> cmd_fragtail("w", "fragtailwindow", "p-value for frag window cut", false,FragTailDefault,"double",cmd);

  // Masking Mobile Element file
  ValueArg<string> cmd_MEIMask("M", "MASKFILE", "Event masking file", false, "", "string", cmd);

  // Mobile element classes
  //ValueArg<string> cmd_elements("E", "element", "names of inserted elements", false, "", "string", cmd);
  ValueArg<string> cmd_elements("S", "specialtag", "prefix of element names in special reference", false, "", "string", cmd);
	
	// mapping Q min 
	int BamZADefault=1;
  ValueArg<int> cmd_BamZA("Z", "Zmode", "bam access mode (1=ZAtag, 2=SortedByReadName) ", false,BamZADefault,"int",cmd);


  //----------------------------------------------------------------------------
  // parse command line and catch possible errors
  //----------------------------------------------------------------------------
  try {
    cmd.parse(argc,argv);
  } 
  catch ( ArgException& e ) { 
    clog << "ERROR: " << e.error() << " " << e.argId() << endl; 
  }
  
  //----------------------------------------------------------------------------
  // assign command line parameters
  //----------------------------------------------------------------------------  
  // name of run control parameter file
  string rcp = cmd_rcp.getValue();
  // name of reference genome fasta file
  // string fasta = cmd_fasta.getValue();
  // name of frag length histogram  file
  string fragHistFile = cmd_fragHist.getValue();
  // name of frag length histogram  file
  // string RDreferenceFile = cmd_RDreference.getValue();
  // name of Mob Masking file
  string MobiMaskFile = cmd_MEIMask.getValue();
  // element list for insertions
  string elements = cmd_elements.getValue();
  // output area
  string outarea = cmd_out.getValue();
  // prefix 
  string prefix = cmd_prefix.getValue();
  // text output flag 
  bool printText = cmd_textout.getValue();
  // report interval
  int dbg = cmd_dbg.getValue();
  // max fragments to extract from input alignment files
  int maxFrag = cmd_maxfragments.getValue();
  // CN selection slosh for PE detected CNVs
  //  int CNslosh = cmd_cnslosh.getValue();
  // Duplicate removal switch 
  int Uniquify = cmd_uniquify.getValue();
  // read sense flip switch 
  // bool readPairSenseFlip = cmd_flip.getValue();  
  // scan input files - no processing or *.span output
  bool  scan = cmd_scan.getValue();
  // name of anchor info file
  string anchorfile = cmd_Anchorfile.getValue();
  // Quality min 
  int Qmin  = cmd_qmin.getValue();
  // fragtailcut
  double fragtailcut = cmd_fragtail.getValue();
   
  //----------------------------------------------------------------------------
	// require ZA tag in bam
  //----------------------------------------------------------------------------
  int BamZA = cmd_BamZA.getValue();

  //----------------------------------------------------------------------------
  // build options
  //----------------------------------------------------------------------------
  bool build = cmd_build.getValue();

  //----------------------------------------------------------------------------
  // specify of pre-built input files
  //----------------------------------------------------------------------------
  string  in = cmd_input.getValue();
  
  if (in.size()<1) { 
	   cout << " hey ... choose inut file -i " << endl;
	   exit(99);
  }    

	//----------------------------------------------------------------------------
	// regexp for allowed contig names
	//----------------------------------------------------------------------------
	string allowContigs = cmd_allowcontigs.getValue();

  //----------------------------------------------------------------------------
  // read parameters from rcp file
  //----------------------------------------------------------------------------
  RunControlParameters pars(rcp);
  pars.setDbg(dbg);

	//----------------------------------------------------------------------------
	// overide Mob Insertion Masking file
	//----------------------------------------------------------------------------
	if (MobiMaskFile.length()>0) {
    pars.setMobiMaskFile(MobiMaskFile);
  }
  // overide element list 
  if (elements.length()>0) {
    vector<string> svect;
    split(svect,elements,",");
    pars.setMobileElements(svect);
  }
  
	//----------------------------------------------------------------------------
	//overide maxFrag if present on command line
	//----------------------------------------------------------------------------
	if (maxFrag>0) {
    pars.setMaxFragments(maxFrag);
  }
 
	//----------------------------------------------------------------------------
	//don't DupRemove unless build or scan
	//----------------------------------------------------------------------------
	if (!(build||scan)) {
    pars.setDupRemove(false);		
	}	

	//----------------------------------------------------------------------------
	//overide Uniquify if present on command line
	//----------------------------------------------------------------------------
	if (Uniquify!=1) {
    pars.setDupRemove(Uniquify);
  }

	//----------------------------------------------------------------------------
	//overide BamZA if present on command line
	//----------------------------------------------------------------------------
	if (BamZA!=BamZADefault) {
    pars.setBamZA(BamZA);
  }
	
	//set Qmin to zero for build 
  /*
	 if (build) {
    pars.setQmin(1);
  }
	*/
	//----------------------------------------------------------------------------
	//overide Qmin if present on command line
  //----------------------------------------------------------------------------
  if (Qmin!=QminDefault) {
    pars.setQmin(Qmin);
  }
	
  //----------------------------------------------------------------------------
  //overide fragtailcut if present on command line
  //----------------------------------------------------------------------------
  if (fabs(fragtailcut-pars.getFragmentLengthWindow())>1e-5) {
    pars.setFragmentLengthWindow(fragtailcut);
  }
	
  //----------------------------------------------------------------------------
  //overide allowContigs if present on command line
  //----------------------------------------------------------------------------
  if (allowContigs.length()>0) {
    pars.setAllowContigsRegex(allowContigs);
  }
	
  //----------------------------------------------------------------------------
  //overide prefix if present on command line
  //----------------------------------------------------------------------------
  if (prefix.length()>0) {
    pars.setPrefix(prefix);
  }
	
  //----------------------------------------------------------------------------
  //overide outarea if present on command line
  //----------------------------------------------------------------------------
  if (outarea.length()>0) {
    pars.setOutputDir(outarea);
  }
	
  //----------------------------------------------------------------------------
  //overide textout if present on command line
  //----------------------------------------------------------------------------
  if (printText) {
    pars.setPrintTextOut(true);
  }
	
  //----------------------------------------------------------------------------
  //overide scan if present on command line
  //----------------------------------------------------------------------------
  if (scan) {
    pars.setDoScanOnly(true);
  } else {
    // overide fragment length limits
    if (fragHistFile.length()>0) {
      C_HistoGroups hg(fragHistFile);
	    if (hg.Groups.size()>0) { 
				if (hg.Groups[0].h.count("LF")>0) {
				 pars.setHistoGroups(hg);
				}
			}
		}
  }
	
  //overide anchor info file if present on command line
  if (anchorfile.length()>0) {
    pars.setAnchorfile(anchorfile);
  }


  //----------------------------------------------------------------------------
  // input paired-end read data
  //----------------------------------------------------------------------------
  C_pairedfiles data(in, pars);
  //----------------------------------------------------------------------------
  //
  int Nset = data.set.size();
  // loop over contigs
  if (Nset>1) {
      cout << "\n number of sets: " << Nset << endl;
  }
  //  
  if (!scan) {  // fragment length scanning already done
		
    for (int set=0; set<Nset; set++) {
      string ss = data.set[set].getSetName();
			if (dbg>0) {
        cout << "\t " << set+1 << ". \t" << ss << endl;
			}
      //------------------------------------------------------------------------
      // sort fragments
      //------------------------------------------------------------------------
      data.set[set].sort();
      //------------------------------------------------------------------------
      // remove redundant fragments
      //------------------------------------------------------------------------
      if ( pars.getDupRemove()>0)  data.set[set].uniquify();
      //------------------------------------------------------------------------
      // Calc raw stats
      //------------------------------------------------------------------------
      data.set[set].calcStats();
      //------------------------------------------------------------------------
      // output
      //------------------------------------------------------------------------
      if (build) {
        data.set[set].write();
        if (printText) {
          data.set[set].printOut();
        }
      } else { 
        //------------------------------------------------------------------------
        // if the input is the combined data from a directory
        // then write out all the input spanner structures to the 
        // output area 
        //------------------------------------------------------------------------      
        if (data.inputType=='D') {
          data.set[set].write();
        }
        //------------------------------------------------------------------------
        // detect SV's
        //------------------------------------------------------------------------
				C_SpannerSV SV(data,pars);
      }   
    }
  }
  //----------------------------------------------------------------------------
  // report
  //----------------------------------------------------------------------------
  if (dbg>0) {
    clog << "Spanner completed." << endl;
  }
  return 0;
}


