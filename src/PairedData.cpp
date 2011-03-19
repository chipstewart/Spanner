/*
 *  PairedData.cpp
 *  Span1
 *
 *  Created by Chip Stewart on 7/24/08.
 *  Copyright 2008 Boston College. All rights reserved.
 *
 */

#include "PairedData.h"
//------------------------------------------------------------------------------
// function to enforce p0 and p1 to be bound within pos0 and pos1 
//------------------------------------------------------------------------------
bool boundlimit(int p0, int p1, int pos0, int pos1) {
  bool within= ((p0>=pos0)&&(p1<=pos1));
  return within;
  /*
  if ((p1<pos0)|(p0>pos1)) {
     return false;
  }
  // shift to region pos0 to pos1
  p0=p0-pos0;
  p1=p1-pos0;
  if (p0<0) { p0=0;}
  if (p1>(pos1-pos0)) {p1=pos1-pos0;}
  return true;
  */
}
//------------------------------------------------------------------------------
// readmap (basic info for one read alignment)
//------------------------------------------------------------------------------
C_readmap::C_readmap()   // Constructor
{
   pos = 0;
   anchor = 0;
   len = 0;
   sense = 0;
   q=0;
   mm=0;
	 q2=0;
	 nmap=0;
	 mob=" ";
}

C_readmap::C_readmap(const C_readmap &copyin)   // Copy constructor to handle pass by value.
{                             
   pos = copyin.pos;
   anchor= copyin.anchor;
   len = copyin.len;
   sense = copyin.sense;
   q = copyin.q;
   mm = copyin.mm;
	 q2=copyin.q2;
	 nmap=copyin.nmap;
	 mob=copyin.mob;
}

C_readmap::C_readmap(unsigned int pos1, unsigned short anchor1,unsigned short len1, char sense1, char q1, char mm1) {
   pos = pos1;
   anchor = anchor1;
   len = len1;
   sense = sense1;
   q = q1;
   mm = mm1;
}

ostream &operator<<(ostream &output, const C_readmap & x)
{
   output << x.anchor << "\t" << x.pos << "\t" << x.len 
          << "\t" << x.sense <<"\t" << int(x.q) <<"\t" << int(x.mm) ;
   return output;
}

C_readmap& C_readmap::operator=(const C_readmap &rhs)
{
   this->pos = rhs.pos;
   this->anchor = rhs.anchor;
   this->len = rhs.len;
   this->sense = rhs.sense;
   this->q = rhs.q;
	 this->mm = rhs.mm;
	 this->nmap = rhs.nmap;
	 this->mob = rhs.mob;
   return *this;
}

int C_readmap::operator==(const C_readmap &rhs) const
{
   if( this->pos != rhs.pos) return 0;
   if( this->anchor != rhs.anchor) return 0;
   if( this->len != rhs.len) return 0;
   if( this->sense != rhs.sense) return 0;
   if( this->q != rhs.q) return 0;
   if( this->mm != rhs.mm) return 0;
 	 if( this->q2 != rhs.q2) return 0;
	 if( this->nmap != rhs.nmap) return 0;
  return 1;
}

// This function is required for built-in STL list functions like sort
int C_readmap::operator<(const C_readmap &rhs) const
{
   if( this->anchor < rhs.anchor ) return 1;
   if( this->anchor == rhs.anchor && this->pos < rhs.pos ) return 1;
   if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->len < rhs.len) return 1;
   if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->len == rhs.len && this->sense < rhs.sense) return 1;
   if( this->anchor == rhs.anchor && this->pos == rhs.pos && this->len == rhs.len 
      && this->sense == rhs.sense && this->q < rhs.q) return 1;
	if( this->anchor == rhs.anchor && this->pos == rhs.pos && this->len == rhs.len 
		 && this->sense == rhs.sense && this->q == rhs.q &&  this->mm < rhs.mm) return 1;
	if( this->anchor == rhs.anchor && this->pos == rhs.pos && this->len == rhs.len 
		 && this->sense == rhs.sense && this->q == rhs.q &&  this->mm == rhs.mm &&  this->nmap < rhs.nmap) return 1;
   return 0;
}

//------------------------------------------------------------------------------
// readmaps  (vector of readmap objects)
//------------------------------------------------------------------------------
C_readmaps::C_readmaps() {
    Nalign  = 0;
    element = 0;
}

ostream &operator<<(ostream &output, const C_readmaps & rhs)
{
    int N=rhs.align.size();
    output << N << "\t " << rhs.element << endl;
    for (int i=0; i<N; i++) {
        output << "\t" << rhs.align[i] << endl;
    }    
    return output;
}
    
C_readmaps::C_readmaps(const C_readmaps &copyin)   
{
  align = copyin.align;                             
  Nalign = copyin.Nalign;                             
  element = copyin.element;
}

C_readmaps& C_readmaps::operator=(const C_readmaps &rhs)
{
   this->align = rhs.align;
   this->Nalign = rhs.Nalign;
   this->element = rhs.element;
   return *this;
}

//------------------------------------------------------------------------------
// pairedread (array of two readmaps)
//------------------------------------------------------------------------------
C_pairedread::C_pairedread() {              // constructor
  Nend=0;
  ReadGroupCode=0;
	Name="";
}
 
ostream &operator<<(ostream &output, const C_pairedread & rhs)
{
    output << rhs.read[0] << endl;
    output << rhs.read[1] << endl;
    return output;
}


//------------------------------------------------------------------------------
// anchorinfo
//------------------------------------------------------------------------------
C_anchorinfo::C_anchorinfo() {               // constructor
 source = "";
 use.clear();
}

//------------------------------------------------------------------------------
// anchorinfo copy constructor
//------------------------------------------------------------------------------
C_anchorinfo::C_anchorinfo(const C_anchorinfo & rhs) {               // constructor
  source = rhs.source;
  L=rhs.L;
  use=rhs.use;
  element=rhs.element;
  names=rhs.names;
}


C_anchorinfo::C_anchorinfo(string & anchorfile) {              // constructor
  // anchor file patterns
  string patternHead("number.+\\s+(\\d+)$");
  //string patternLine("(\\d+)\\.\\s+(\\S+)\\s+(\\d+)\\s*$");
	string patternLine("^\\s*(\\d+)\\.\\s+(\\w+)\\s+(\\d+)\\s*$");
  source=anchorfile;
  fstream file1;
  file1.open(anchorfile.c_str(), ios::in);	
  if (!file1) {
    cerr << "Unable to open anchor file: " << anchorfile << endl;
    exit(1);
  }
  string line;
  string match,match2,match3;
  int N=0;
  while (getline(file1, line)) {		
		//================
    if (RE2::FullMatch(line.c_str(),patternHead.c_str(),&match) ) {
        N = string2Int(match);
    } else if ( RE2::FullMatch(line.c_str(),patternLine.c_str(),&match,&match2,&match3)  ) {
        string cn1 = match2;
        int L1 = string2Int(match3);
        L[cn1]=L1;
        names.push_back(cn1);
    }
  }
  use.resize(N,true);
}

C_anchorinfo::C_anchorinfo(char t) {               // constructor for build 36.1
  source = "";
  if (t==0) {
    source = "AB";
    for (int i = 0; i<25; i++) {
        names.push_back(ABCHROMOSOMENAME[i]); 
        L[ABCHROMOSOMENAME[i]]=ABCHROMOSOMELENGTH[i];
    }
    cout << "\t ... using default ncbi build 36.1 anchor info" << endl;
    use.resize(25,true);
  }
}

//------------------------------------------------------------------------------
// Mosaik anchor info filler
//------------------------------------------------------------------------------
C_anchorinfo::C_anchorinfo(Mosaik::CAlignmentReader & ar1) {             
  source = "MSK";
  
  // retrieve the reference sequence data
  vector<Mosaik::ReferenceSequence> ref1=ar1.GetReferenceData();
  int Nref = ref1.size();

  // reset 
  names.clear();
  L.clear();
  use.clear();
  element.clear();
  
  // allocate names and use
  names.resize(Nref,"");
  use.resize(Nref,true);
  element.resize(Nref,0);

  // loop over ref
  for (int id=0; id<Nref; id++) {
	names[id]=ref1[id].Name;
    L[names[id]]=ref1[id].NumBases;	
    // turn off SV detection from "NT_" (random) contigs
    size_t found=names[id].find("NT_");
    if (found==0) use[id]=false;  
	found=names[id].find("NC_");
	if (found==0) use[id]=false;  
	found=names[id].find("GL0",0,3);
	if (found==0) use[id]=false;  
  }
	
  if (!anchorElements()) {
    cerr << " anchorelements not found " << endl;
  }
}

/*
//------------------------------------------------------------------------------
// Mosaik anchor info filler
//------------------------------------------------------------------------------
C_anchorinfo::C_anchorinfo(BamReader & br1) {             
  source = "BAM";
  string ht=br1.GetHeaderText();
  // retrieve the reference sequence data
  vector<Mosaik::ReferenceSequence> ref1=GetAnchorInfoFromBamHeader(ht);
  int Nref = ref1.size();

  // reset 
  names.clear();
  L.clear();
  use.clear();
  element.clear();
  
  // allocate names and use
  names.resize(Nref,"");
  use.resize(Nref,true);
  element.resize(Nref,0);

  // loop over ref
  for (int id=0; id<Nref; id++) {
		names[id]=ref1[id].Name;
    L[names[id]]=ref1[id].NumBases;	
  }
  if (!anchorElements()) {
    cerr << " anchorelements not found " << endl;
  }
}

vector<Mosaik::ReferenceSequence> C_anchorinfo::GetAnchorInfoFromBamHeader(string & ht) 
{
//@SQ	SN:1	LN:247249719
  vector<Mosaik::ReferenceSequence> References;
  string line,thing;
  stringstream stream1(ht); //stringstream::out | stringstream::in);
  //stream1 << ht;
  while( getline(stream1, line) ) {
      if (line.find("@SQ")!=string::npos) {      
        cout << line << "\n";
        Mosaik::ReferenceSequence RS1;
        stringstream stream2(line); 
        while( getline(stream2,thing,'\t') ) {
          string token="SN:";
          size_t f1=thing.find(token);
          if (f1!=string::npos) {      
            RS1.Name=thing.substr(f1+token.length());
          }
          token="LN:";
          f1=thing.find(token);
          if (f1!=string::npos) {      
            string nb=thing.substr(f1+token.length());
            RS1.NumBases=string2Int(nb);          
          }
        }
        References.push_back(RS1);
      }
  }
  return References;
}

*/

int C_anchorinfo::operator==(const C_anchorinfo &rhs) const
{
   if (this->L.size()!=rhs.L.size() ) return 0;
    
   /* another day
   for(size_t i=0; i<names.size(); ++i) {
      unsigned int L1 =  L[names[i]];
      if (rhs.L.count(names[i])>0) return 0;
      if (rhs.L[names[i]] != L1) return 0;
    }
    */
    
    return 1;
}


//------------------------------------------------------------------------------
// Bam anchor info filler
//------------------------------------------------------------------------------
//C_anchorinfo::C_anchorinfo(BamReader & ar1) {             
C_anchorinfo::C_anchorinfo(BamMultiReader & ar1) {             
	source = "BAM";
	names.clear();
	L.clear();
	use.clear();
	element.clear();
	
	int maxId = 0;
	RefVector refSeq1=ar1.GetReferenceData();
	
	int Nref = ar1.GetReferenceCount();
	
	names.resize(Nref,"");
	use.resize(Nref,0);
	element.resize(Nref,0);
	
	
	for (RefVector::const_iterator ri = refSeq1.begin();
			 ri != refSeq1.end(); ri++) {
		
		// retrieve reference info
		RefData rd = *ri;
		
		// get member data
		string refName = rd.RefName;
		
		L[rd.RefName]=rd.RefLength;
		names[maxId]=refName;
		// turn off detection from contigs that begin with NT_, NC_ or GL0
		size_t found1=refName.find("NT_",0,3);
		size_t found2=refName.find("NC_",0,3);
		size_t found3=refName.find("GL0",0,3);
		use[maxId]=(found1!=0)&&(found2!=0)&&(found3!=0);
		
		maxId++;
		
	}
	
	/*
	 if (!anchorElements()) {
	 cerr << " anchorelements not found " << endl;
	 }
	 
	 for (int i=0; i<use.size(); i++) {
	 if (use[i]) printf("%d\t%d\t%s\t%d\n",i+1,use[i],names[i].c_str(),L[names[i]]);
	 }
	 */
}

void C_anchorinfo::anchorlimit(string & allow) {         
  string patternAllowContigsRegex(allow);  
  // regex match
  string name=" ";
  for (size_t i = 0; i<names.size(); i++) {
      name= names[i];
      // skip this ContigName if not present in AllowContigs
      if (!  RE2::PartialMatch(name.c_str(),patternAllowContigsRegex.c_str())) {
        use[i]=false;
      }
  }
}

bool C_anchorinfo::anchorElements() {         

    string e0[] = {"moblist_ALU","moblist_L1","moblist_SVA","moblist_ERV"};
    vector<string> elements1(e0, e0 + 4);
    
    return(anchorElements(elements1));
}


//------------------------------------------------------------------------------
// add list of reference elements to anchorinfo  
//------------------------------------------------------------------------------
bool C_anchorinfo::anchorElements(vector<string> & elements1) {         
  
  // init
  bool done=false;
  string name=" ";
  string elem=" ";
  element.clear();
  element.resize(names.size(),-1);
  
  // loop over anchors
  for (size_t i = 0; i<names.size(); i++) {
      name= names[i].substr(0,10);
      
      // loop over elements
      for (size_t e = 0; e<elements1.size(); e++) {
        elem= elements1[e].substr(0,10);

        // set element index 
        if (name==elem) {
           element[i]=e;

           // do not analyze coverage within these elements
           use[i]=false;
           done=true;
           break;
        }

      }

   }
   return done;
}  

//------------------------------------------------------------------------------
// get minimum element anchor index - return 999 if no elements
//------------------------------------------------------------------------------
unsigned int C_anchorinfo::anchorMinElement() {
  unsigned int a1;
  for(a1=0; a1<element.size(); a1++) {
    if (element[a1]>=0) {
      return a1;
    }
  }
  a1=999;
  return a1;
}


//------------------------------------------------------------------------------
// get anchor index - return -1 if no anchor name found
//------------------------------------------------------------------------------
char C_anchorinfo::anchorIndex(string & name1){
	unsigned int a1;
	for(a1=0; a1<names.size(); a1++) {
		if (names[a1].compare(name1)==0) {
			return a1;
		}
	}
	a1=-1;
	return a1;
}

char anchorIndex(string &);   // return anchor index given anchor name 



ostream &operator<<(ostream &output,  C_anchorinfo & a)
{
    output << "number of anchors:\t" << a.names.size() << endl ;      
    for(size_t i=0; i<a.names.size(); ++i) {
        string s = a.names[i];
        unsigned int L = a.L[s];
        output << "\t" << i <<".\t" << s << "\t" << L << endl; //" " << a.use[i] << endl ;  
    }
    return output;
}

void C_anchorinfo::push_anchor (string & name1, unsigned int L1) {
  names.push_back(name1);
  L[name1]=L1;
}

void C_anchorinfo::printAnchorInfo(string & name) {
// format: optimized to loadFragments.m matlab script  
  // open output binary file. bomb if unable to open
  fstream output(name.c_str(), ios::out);
  if (!output) {
      cerr << "Unable to open anchor info output file: " << name << endl;
      return;
  }
  int V = 207;
  output << "Anchor Information, version " << V << endl;
  output << "\t #.\t Anchor\t Length\t use" << endl;
  output << *this;
  output << endl ;
  output.close();
}

//------------------------------------------------------------------------------
// Library info 
//------------------------------------------------------------------------------
C_libraryinfo::C_libraryinfo()                          // constructor default
{
    LM=0;
    tailcut=0;
    LMlow=0;
    LMhigh=0;    
    LR=0;
    LRmin=0;
    LRmax=0;
    NPair=0;
    NSingle=0;
    NPairRedundant=0;
    NSingleRedundant=0;
}

C_libraryinfo::C_libraryinfo(const C_libraryinfo & rhs)     // copy constructor
{
    LM=rhs.LM;
    tailcut=rhs.tailcut;
    LMlow=rhs.LMlow;
    LMhigh=rhs.LMhigh;    
    LR=rhs.LR;
    LRmin=rhs.LRmin;
    LRmax=rhs.LR;
    NPair=rhs.NPair;
    NSingle=rhs.NSingle;
    NPairRedundant=rhs.NPairRedundant;
    NSingleRedundant=rhs.NSingleRedundant;
    fragHist=rhs.fragHist; 
    readLengthHist=rhs.readLengthHist;       
    Info.MedianFragmentLength=rhs.Info.MedianFragmentLength;
    Info.ReadGroupCode=rhs.Info.ReadGroupCode;
    Info.SequencingTechnology=rhs.Info.SequencingTechnology;
    Info.CenterName=rhs.Info.CenterName;
    Info.Description=rhs.Info.Description;
    Info.LibraryName=rhs.Info.LibraryName;
    Info.PlatformUnit=rhs.Info.PlatformUnit;
    Info.ReadGroupID=rhs.Info.ReadGroupID;
    Info.SampleName=rhs.Info.SampleName;    
}

ostream &operator<<(ostream &output,  C_libraryinfo & lib)
{
    output << "Read Group:\t" << lib.Info.ReadGroupCode;
    output << "\t ID:\t" << lib.Info.ReadGroupID  << endl ;      
    output << "\t "<< lib.Info.SequencingTechnology;
    output << "\t "<< lib.Info.CenterName;
    output << "\t "<< lib.Info.Description << endl;
    output << "\t "<< lib.Info.LibraryName;
    output << "\t "<< lib.Info.PlatformUnit;
    output << "\t "<< lib.Info.SampleName;
    output << "\t Insert:\t" << lib.Info.MedianFragmentLength;
    output << "\t "<< lib.LM;
    output << "\t "<< lib.LMlow;
    output << "\t "<< lib.LMhigh << endl;
    output << "\t Read:\t"<< lib.LR;
    output << "\t "<< lib.LRmin;
    output << "\t "<< lib.LRmax << endl;
    output << "\t NPE:\t"<< lib.NPair;
    output << "\t Redundnat:\t"<< lib.NPairRedundant << endl;
    output << "\t NSE:\t"<< lib.NSingle;
    output << "\t Redundnat:\t"<< lib.NSingleRedundant << endl;
    return output;
}

//------------------------------------------------------------------------------
// Libraries  
//------------------------------------------------------------------------------
C_libraries::C_libraries()              // constructor 
{
}    

//------------------------------------------------------------------------------
// binary file input
//------------------------------------------------------------------------------
C_libraries::C_libraries(string & infilename)      // load lib info 
{
  fstream input(infilename.c_str(), ios::in  | ios::binary);
  if (!input) {
      cerr << "Unable to open library.span file: " << infilename << endl;
  }
  C_headerSpan h(input);

  C_librarymap::iterator it;
  C_libraryinfo lib1;    

  int N= h.N;
  //----------------------------------------------------------------------------
  // fixed length part of lib record
  //----------------------------------------------------------------------------  
  for(int i=0; i<N; i++)  {
    input.read(reinterpret_cast< char *>(&lib1.Info.ReadGroupCode), sizeof(unsigned int));
    input.read(reinterpret_cast< char *>(&lib1.Info.MedianFragmentLength), sizeof(unsigned int));
    input.read(reinterpret_cast< char *>(&lib1.Info.SequencingTechnology), sizeof(unsigned short));
    input.read(reinterpret_cast< char *>(&lib1.LM), sizeof(int));
    input.read(reinterpret_cast< char *>(&lib1.LMlow), sizeof(int));
    input.read(reinterpret_cast< char *>(&lib1.LMhigh), sizeof(int));
    input.read(reinterpret_cast< char *>(&lib1.tailcut), sizeof(double));
    input.read(reinterpret_cast< char *>(&lib1.LR), sizeof(double));
    input.read(reinterpret_cast< char *>(&lib1.LRmin), sizeof(int));
    input.read(reinterpret_cast< char *>(&lib1.LRmax), sizeof(int));
    input.read(reinterpret_cast< char *>(&lib1.NPair), sizeof(int));
    input.read(reinterpret_cast< char *>(&lib1.NSingle), sizeof(int));
    input.read(reinterpret_cast< char *>(&lib1.NPairRedundant), sizeof(int));
    input.read(reinterpret_cast< char *>(&lib1.NSingleRedundant), sizeof(int));
    this->libmap[lib1.Info.ReadGroupCode]=lib1;
  }
  //----------------------------------------------------------------------------
  // variable length part of lib record
  //----------------------------------------------------------------------------  
  int nchar = 0;
  // buffer for loading strings
  char buff[512];
  for ( it=this->libmap.begin() ; it != this->libmap.end(); it++ )
  {  
    lib1 = (*it).second; 
    input.read(reinterpret_cast< char *>(&lib1.Info.ReadGroupCode), sizeof(unsigned int));
    input.read(reinterpret_cast< char *>(&lib1.Info.MedianFragmentLength), sizeof(unsigned int));
    // ReadGroupID
    input.read(reinterpret_cast < char * > (&nchar), sizeof(int));
    input.read(buff, nchar);
    // convert c-style string buff to string contigName 
    buff[nchar]=0;
    lib1.Info.ReadGroupID=buff;
    // CenterName
    input.read(reinterpret_cast < char * > (&nchar), sizeof(int));
    input.read(buff, nchar);
    buff[nchar]=0;
    lib1.Info.CenterName=buff;
    // Description
    input.read(reinterpret_cast < char * > (&nchar), sizeof(int));
    input.read(buff, nchar);
    buff[nchar]=0;
    lib1.Info.Description=buff;
    // LibraryName
    input.read(reinterpret_cast < char * > (&nchar), sizeof(int));
    input.read(buff, nchar);
    buff[nchar]=0;
    lib1.Info.LibraryName=buff;
    // PlatformUnit
    input.read(reinterpret_cast < char * > (&nchar), sizeof(int));
    input.read(buff, nchar);
    buff[nchar]=0;
    lib1.Info.PlatformUnit=buff;
    // SampleName
    input.read(reinterpret_cast < char * > (&nchar), sizeof(int));
    input.read(buff, nchar);
    buff[nchar]=0;
    lib1.Info.SampleName=buff;
    // LF distribution
    input.read(reinterpret_cast < char * > (&nchar), sizeof(int));
    input.read(buff, nchar);
    buff[nchar]=0;
    lib1.fragHist.title=buff;
    input.read(reinterpret_cast< char *>(&lib1.fragHist.Ntot), sizeof(double));
    input.read(reinterpret_cast< char *>(&lib1.fragHist.mean), sizeof(double));
    input.read(reinterpret_cast< char *>(&lib1.fragHist.std), sizeof(double));
    input.read(reinterpret_cast< char *>(&lib1.fragHist.median), sizeof(double));
    input.read(reinterpret_cast< char *>(&lib1.fragHist.Nin), sizeof(double));
    input.read(reinterpret_cast< char *>(&lib1.fragHist.Nunder), sizeof(double));
    input.read(reinterpret_cast< char *>(&lib1.fragHist.Nover), sizeof(double));
    input.read(reinterpret_cast< char *>(&lib1.fragHist.Nbin), sizeof(int));
    lib1.fragHist.xc.resize(lib1.fragHist.Nbin);
    lib1.fragHist.n.resize(lib1.fragHist.Nbin);
    lib1.fragHist.c.resize(lib1.fragHist.Nbin);
    for (int i=0; i<lib1.fragHist.Nbin; i++) {
      input.read(reinterpret_cast< char *>(&lib1.fragHist.xc[i]), sizeof(double));
      input.read(reinterpret_cast< char *>(&lib1.fragHist.n[i]), sizeof(double));
      input.read(reinterpret_cast< char *>(&lib1.fragHist.c[i]), sizeof(double));
    }
    lib1.fragHist.dx=lib1.fragHist.xc[1]-lib1.fragHist.xc[0];      
    lib1.fragHist.xlow=lib1.fragHist.xc[0]-lib1.fragHist.dx/2.0;      
    lib1.fragHist.xhigh=lib1.fragHist.xc[lib1.fragHist.Nbin-1]+lib1.fragHist.dx/2.0;      
    lib1.fragHist.Finalize();         
 
    // LR distribution
    input.read(reinterpret_cast < char * > (&nchar), sizeof(int));
    input.read(buff, nchar);
    buff[nchar]=0;
    lib1.readLengthHist.title=buff;
    input.read(reinterpret_cast<  char *>(&lib1.readLengthHist.Ntot), sizeof(double));
    input.read(reinterpret_cast<  char *>(&lib1.readLengthHist.mean), sizeof(double));
    input.read(reinterpret_cast<  char *>(&lib1.readLengthHist.std), sizeof(double));
    input.read(reinterpret_cast<  char *>(&lib1.readLengthHist.median), sizeof(double));
    input.read(reinterpret_cast<  char *>(&lib1.readLengthHist.Nin), sizeof(double));
    input.read(reinterpret_cast<  char *>(&lib1.readLengthHist.Nunder), sizeof(double));
    input.read(reinterpret_cast<  char *>(&lib1.readLengthHist.Nover), sizeof(double));
    input.read(reinterpret_cast<  char *>(&lib1.readLengthHist.Nbin), sizeof(int));
    lib1.readLengthHist.xc.resize(lib1.readLengthHist.Nbin);
    lib1.readLengthHist.n.resize(lib1.readLengthHist.Nbin);
    lib1.readLengthHist.c.resize(lib1.readLengthHist.Nbin);
    for (int i=0; i<lib1.readLengthHist.Nbin; i++) {
        input.read(reinterpret_cast<  char *>(&lib1.readLengthHist.xc[i]), sizeof(double));
        input.read(reinterpret_cast<  char *>(&lib1.readLengthHist.n[i]), sizeof(double));
        input.read(reinterpret_cast<  char *>(&lib1.readLengthHist.c[i]), sizeof(double));
    }
    lib1.readLengthHist.dx=lib1.readLengthHist.xc[1]-lib1.readLengthHist.xc[0];      
    lib1.readLengthHist.xlow=lib1.readLengthHist.xc[0]-lib1.readLengthHist.dx/2.0;      
    lib1.readLengthHist.xhigh=lib1.readLengthHist.xc[lib1.readLengthHist.Nbin-1]+lib1.readLengthHist.dx/2.0;      
    lib1.readLengthHist.Finalize();         
         
    this->libmap[lib1.Info.ReadGroupCode]=lib1;
		ReadGroupID2Code[lib1.Info.ReadGroupID]=lib1.Info.ReadGroupCode;

  }
  //----------------------------------------------------------------------------
  // anchor info
  //----------------------------------------------------------------------------
  int Na;
  input.read(reinterpret_cast<  char *>(&Na), sizeof(int));
  anchorinfo.names.resize(Na);
  anchorinfo.use.resize(Na);
  for(size_t i=0; int(i)<Na; ++i) {
    input.read(reinterpret_cast < char * > (&nchar), sizeof(int));
    input.read(buff, nchar);
    buff[nchar]=0;
    anchorinfo.names[i]=buff;
    unsigned int L;
    input.read(reinterpret_cast<char *>(&L), sizeof(unsigned int));
    anchorinfo.L[anchorinfo.names[i]]=L;
    bool use;
    input.read(reinterpret_cast< char *>(&use), sizeof(bool));
    anchorinfo.use[i]=use;
  }
 
   input.close();

}

C_libraries::C_libraries( const C_libraries & rhs)              // constructor 
{
  libmap=rhs.libmap;
  anchorinfo=rhs.anchorinfo;
	ReadGroupID2Code=rhs.ReadGroupID2Code;

}    

C_libraries::C_libraries(Mosaik::CAlignmentReader & ar1) // constructor - load from Mosaik file    
{
  vector<Mosaik::ReadGroup> readGroups;
	readGroups=ar1.GetReadGroups();
  int N = readGroups.size();
  C_libraryinfo lib1;
  for (int i=0; i<N; i++) {
    unsigned int ReadGroupCode=readGroups[i].ReadGroupCode;
		libmap[ReadGroupCode].Info=readGroups[i];
		ReadGroupID2Code[readGroups[i].ReadGroupID]=ReadGroupCode;
  }
}    

//C_libraries::C_libraries(BamReader  & br1)              // constructor - load from BAM file    
C_libraries::C_libraries(BamMultiReader  & br1)              // constructor - load from BAM file    
{
  string ht=br1.GetHeaderText();
  vector<Mosaik::ReadGroup> readGroups=GetReadGroups(ht);
  int N = readGroups.size();
  if (N<1) {
     cerr << "Need read group info in BAM file header" << endl;
     exit(-1);
  }
  
  // C_libraryinfo lib1;
  for (int i=0; i<N; i++) {
    unsigned int ReadGroupCode=readGroups[i].ReadGroupCode;
		libmap[ReadGroupCode].Info=readGroups[i];
		ReadGroupID2Code[readGroups[i].ReadGroupID]=ReadGroupCode;
  }
}    

vector<Mosaik::ReadGroup> C_libraries::GetReadGroups(string & ht) 
{
// @RG	ID:ERR001607	PL:ILLUMINA	PU:1000G-mpimg-081010-2_3	LB:NA 19210.11	PI:200	SM:NA19210	CN:MPIMG
  vector<Mosaik::ReadGroup> readGroups;
  string line,thing;
  stringstream stream1(ht); //stringstream::out | stringstream::in);
  //stream1 << ht;
  while( getline(stream1, line) ) {
      if (line.find("@RG")!=string::npos) {      
        cout << line << "\n";
        Mosaik::ReadGroup RG1;
        stringstream stream2(line); //stringstream::out | stringstream::in);
        //string tab1="\t";
        while( getline(stream2,thing,'\t') ) {
          string token="ID:";
          size_t f1=thing.find(token);
          if (f1!=string::npos) {      
            RG1.ReadGroupID=thing.substr(f1+token.length());
          }
          token="PL:";
          f1=thing.find(token);
          if (f1!=string::npos) {   
            string technology=upperCase(thing.substr(f1+token.length()));
            if (technology.find("ILLUMINA")!=string::npos) { RG1.SequencingTechnology=ST_ILLUMINA; }
            if (technology.find("454")!=string::npos) { RG1.SequencingTechnology=ST_454; }
            if (technology.find("HELICOS")!=string::npos) { RG1.SequencingTechnology=ST_HELICOS; } 
            if (technology.find("SOLID")!=string::npos) { RG1.SequencingTechnology=ST_SOLID; }
          }
          token="LB:";
          f1=thing.find(token);
          if (f1!=string::npos) {      
            RG1.LibraryName=thing.substr(f1+token.length());
          }
          token="PU:";
          f1=thing.find(token);
          if (f1!=string::npos) {      
            RG1.PlatformUnit=thing.substr(f1+token.length());
          }
          token="CN:";
          f1=thing.find(token);
          if (f1!=string::npos) {      
            RG1.CenterName=thing.substr(f1+token.length());
          }
          token="SM:";
          f1=thing.find(token);
          if (f1!=string::npos) {      
            RG1.SampleName=thing.substr(f1+token.length());
          }
          token="PI:";
          f1=thing.find(token);
          if (f1!=string::npos) {      
            string fraglen=thing.substr(f1+token.length());
            RG1.MedianFragmentLength=string2Int(fraglen);          
          }
        }
        RG1.ReadGroupCode=CSHA1::GenerateReadGroupCode(RG1.ReadGroupID, RG1.SampleName);
        readGroups.push_back(RG1);
      }
  }
  return readGroups;
}

vector<char> C_libraries::getSequencingTechnology() 
{
	C_librarymap::iterator it;
	int N = libmap.size();
	
	vector<char> p;
	p.resize(N);
	int i=0;
	for ( it=libmap.begin() ; it != libmap.end(); it++ )
	{
		//unsigned int ReadGroupCode=(*it).first;
		C_libraryinfo lib1 = (*it).second; 
		p[i]=lib1.Info.SequencingTechnology;
		i++;
	}
    return p;
}

ostream &operator<<(ostream &output,  C_libraries & libs)
{
  C_librarymap::iterator it;
  int N = libs.libmap.size();

  output << " Number of Libraries:\t " <<  N<< endl;
  for ( it=libs.libmap.begin() ; it != libs.libmap.end(); it++ )
  {
    //unsigned int ReadGroupCode=(*it).first;
    C_libraryinfo lib1 = (*it).second; 
    output << lib1 << endl;
   }
    return output;
}

void C_libraries::printLibraryInfo(string & file1)      // print lib info 
{
  fstream output(file1.c_str(), ios::out);
  if (!output) {
      cerr << "Unable to open library info output file: " << file1 << endl;
      return;
  }
  output << *this << endl;
  output.close();
}

//------------------------------------------------------------------------------
// binary file output
//------------------------------------------------------------------------------
void C_libraries::writeLibraryInfo(string & outfilename, string & setName)      // write lib info 
{
  C_librarymap::iterator it;
  //int N = libmap.size();
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out | ios::binary);
  if (!output) {
		cerr << "Unable to open file: " << outfilename << endl;
		return;
  }
  C_headerSpan h;
  h.setName = setName;
  h.contigName = " ";
  string typeName="library.span";
  if  (outfilename.find(typeName)==string::npos)  {
    cerr << "bad file name:\t" << outfilename << endl;
    exit(-1);
  }
  h.typeName = typeName;
  h.reclen =   2*sizeof(unsigned int)+9*sizeof(int)+2*sizeof(double)+sizeof(short);
  h.N = this->libmap.size();
  h.write(output);
	
  //----------------------------------------------------------------------------
  // fixed length part of lib record
  //----------------------------------------------------------------------------  
  for ( it=libmap.begin() ; it != libmap.end(); it++ )
  {
    unsigned int ReadGroupCode=(*it).first;
    C_libraryinfo lib1 = (*it).second; 
    output.write(reinterpret_cast<const char *>(&ReadGroupCode), sizeof(unsigned int));
    output.write(reinterpret_cast<const char *>(&lib1.Info.MedianFragmentLength), sizeof(unsigned int));
    output.write(reinterpret_cast<const char *>(&lib1.Info.SequencingTechnology), sizeof(unsigned short));
    output.write(reinterpret_cast<const char *>(&lib1.LM), sizeof(int));
    output.write(reinterpret_cast<const char *>(&lib1.LMlow), sizeof(int));
    output.write(reinterpret_cast<const char *>(&lib1.LMhigh), sizeof(int));
    output.write(reinterpret_cast<const char *>(&lib1.tailcut), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.LR), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.LRmin), sizeof(int));
    output.write(reinterpret_cast<const char *>(&lib1.LRmax), sizeof(int));
    output.write(reinterpret_cast<const char *>(&lib1.NPair), sizeof(int));
    output.write(reinterpret_cast<const char *>(&lib1.NSingle), sizeof(int));
    output.write(reinterpret_cast<const char *>(&lib1.NPairRedundant), sizeof(int));
    output.write(reinterpret_cast<const char *>(&lib1.NSingleRedundant), sizeof(int));
  }
  //----------------------------------------------------------------------------
  // variable length part of lib record
  //----------------------------------------------------------------------------  
  for ( it=libmap.begin() ; it != libmap.end(); it++ )
  {  
    //unsigned int ReadGroupCode=(*it).first;
    C_libraryinfo lib1 = (*it).second; 
    output.write(reinterpret_cast<const char *>(&lib1.Info.ReadGroupCode), sizeof(unsigned int));
    output.write(reinterpret_cast<const char *>(&lib1.Info.MedianFragmentLength), sizeof(unsigned int));
    // ReadGroupID
    const char * cC = lib1.Info.ReadGroupID.c_str();    
    int cL = strlen(cC)+1;
    output.write(reinterpret_cast<const char *>(&cL), sizeof(int));
    output.write((const char *)(cC), cL * sizeof(char));
    // CenterName
    cC = lib1.Info.CenterName.c_str();    
    cL = strlen(cC)+1;
    output.write(reinterpret_cast<const char *>(&cL), sizeof(int));
    output.write((const char *)(cC), cL * sizeof(char));
    // Description
    cC = lib1.Info.Description.c_str();    
    cL = strlen(cC)+1;
    output.write(reinterpret_cast<const char *>(&cL), sizeof(int));
    output.write((const char *)(cC), cL * sizeof(char));
    // LibraryName
    cC = lib1.Info.LibraryName.c_str();    
    cL = strlen(cC)+1;
    output.write(reinterpret_cast<const char *>(&cL), sizeof(int));
    output.write((const char *)(cC), cL * sizeof(char));
    // PlatformUnit
    cC = lib1.Info.PlatformUnit.c_str();    
    cL = strlen(cC)+1;
    output.write(reinterpret_cast<const char *>(&cL), sizeof(int));
    output.write((const char *)(cC), cL * sizeof(char));
    // SampleName
    cC = lib1.Info.SampleName.c_str();    
    cL = strlen(cC)+1;
    output.write(reinterpret_cast<const char *>(&cL), sizeof(int));
    output.write((const char *)(cC), cL * sizeof(char));
    // LF distribution
    cC = lib1.fragHist.title.c_str();    
    cL = strlen(cC)+1;
    output.write(reinterpret_cast<const char *>(&cL), sizeof(int));
    output.write((const char *)(cC), cL * sizeof(char));
    output.write(reinterpret_cast<const char *>(&lib1.fragHist.Ntot), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.fragHist.mean), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.fragHist.std), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.fragHist.median), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.fragHist.Nin), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.fragHist.Nunder), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.fragHist.Nover), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.fragHist.Nbin), sizeof(int));
    for (int i=0; i<lib1.fragHist.Nbin; i++) {
			//if (lib1.fragHist.n[i]>0) {
			output.write(reinterpret_cast<const char *>(&lib1.fragHist.xc[i]), sizeof(double));
			output.write(reinterpret_cast<const char *>(&lib1.fragHist.n[i]), sizeof(double));
			output.write(reinterpret_cast<const char *>(&lib1.fragHist.c[i]), sizeof(double));
			//}
    }
    // LR distribution
    cC = lib1.readLengthHist.title.c_str();    
    cL = strlen(cC)+1;
    output.write(reinterpret_cast<const char *>(&cL), sizeof(int));
    output.write((const char *)(cC), cL * sizeof(char));
    output.write(reinterpret_cast<const char *>(&lib1.readLengthHist.Ntot), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.readLengthHist.mean), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.readLengthHist.std), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.readLengthHist.median), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.readLengthHist.Nin), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.readLengthHist.Nunder), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.readLengthHist.Nover), sizeof(double));
    output.write(reinterpret_cast<const char *>(&lib1.readLengthHist.Nbin), sizeof(int));
    for (int i=0; i<lib1.readLengthHist.Nbin; i++) {
			//if (lib1.readLengthHist.n[i]>0) {
			output.write(reinterpret_cast<const char *>(&lib1.readLengthHist.xc[i]), sizeof(double));
			output.write(reinterpret_cast<const char *>(&lib1.readLengthHist.n[i]), sizeof(double));
			output.write(reinterpret_cast<const char *>(&lib1.readLengthHist.c[i]), sizeof(double));
			//}
    }
  }
  //----------------------------------------------------------------------------
  // anchor info
  //----------------------------------------------------------------------------
  int Na=anchorinfo.names.size();
  output.write(reinterpret_cast<const char *>(&Na), sizeof(int));
  for(size_t i=0; int(i)<Na; ++i) {
    string s = anchorinfo.names[i];
    unsigned int L = anchorinfo.L[s];
    bool use = anchorinfo.use[i];
    const char * cC = s.c_str();    
    int cL = strlen(cC)+1;
    output.write(reinterpret_cast<const char *>(&cL), sizeof(int));
    output.write((const char *)(cC), cL * sizeof(char));
    output.write(reinterpret_cast<const char *>(&L), sizeof(unsigned int));
    output.write(reinterpret_cast<const char *>(&use), sizeof(bool));
  }
  output.close();
}

// max fragment length in this set of libraries
int C_libraries::maxLF()              
{
  C_librarymap::iterator it;
  int LMmax=0;
  for ( it=libmap.begin() ; it != libmap.end(); it++ )
  {
    //unsigned int ReadGroupCode=(*it).first;
    C_libraryinfo lib1 = (*it).second; 
    int LMmax1=lib1.LMhigh;
    if (LMmax1 > LMmax) LMmax=LMmax1;
  }
  return LMmax;
}    
    
// return tailcut from first library    
double C_libraries::getTailcut()
{
  C_librarymap::iterator it;
  double t=0;
  for ( it=libmap.begin() ; it != libmap.end(); it++ )
  {
    C_libraryinfo lib1 = (*it).second; 
    t =lib1.tailcut;
    if (t>0) break;
  }
  return t;
}

// set frag limits LMlow, LMhigh by new tailcut 
void C_libraries::resetFragLimits(double tailcut)
{
  C_librarymap::iterator it;
  double lowf = (1.0-tailcut/100.0)/2.0;  // half tail on low side
  double highf = (1.0-lowf); // half tail on high side     
  for ( it=libmap.begin() ; it != libmap.end(); it++ )
  {
    unsigned int ReadGroupCode=(*it).first;
    int LMlow  = libmap[ReadGroupCode].fragHist.p2xTrim(lowf);
    int LMhigh = libmap[ReadGroupCode].fragHist.p2xTrim(highf);
		if (tailcut<0) { 
			LMlow=100000;
			LMhigh=-100000;
		}
		
		// fill LM 
		if (libmap[ReadGroupCode].LM==0) {
			libmap[ReadGroupCode].LM = libmap[ReadGroupCode].fragHist.mode1;
		}
		
    libmap[ReadGroupCode].LMlow=LMlow;
    libmap[ReadGroupCode].LMhigh=LMhigh;
    libmap[ReadGroupCode].tailcut=tailcut;
		
		
  }
}

// read fractions for each sample
map<string, double, less<string> > C_libraries::readFractionSamples() 
{
  map<string, double, less<string> >  rX;
  C_librarymap::iterator it;
  double N=0;
  for ( it=libmap.begin() ; it != libmap.end(); it++ )
  {
    unsigned int ReadGroupCode=(*it).first;
    double NR1  = double(libmap[ReadGroupCode].NPair);
    string SAM  = libmap[ReadGroupCode].Info.SampleName;
    rX[SAM]+=NR1;
    N+=NR1;
  }
  map<string, double, less<string> >::iterator ir;
  for ( ir=rX.begin() ; ir != rX.end(); ir++ )
  {
    string SAM=(*ir).first;
    rX[SAM]=rX[SAM]/N;
  }
  return rX;
}
    

//------------------------------------------------------------------------------
// localpair  (paired-read on 1 anchor)
//------------------------------------------------------------------------------
C_localpair::C_localpair() {
   pos=0; 
   lm=0;
   anchor=0;
   len1=0;
   len2=0;
   orient='-';
   q1=0;
   q2=0;
   mm1=0;
   mm2=0;
   constrain = 0;
   ReadGroupCode=0;
}

C_localpair::C_localpair(const C_localpair &copyin)   // Copy constructor to handle pass by value.
{                             
   pos = copyin.pos;
   anchor= copyin.anchor;
   lm = copyin.lm;
   len1 = copyin.len1;
   len2 = copyin.len2;
   orient = copyin.orient;
   q1 = copyin.q1;
   q2 = copyin.q2;
   mm1 = copyin.mm1;
   mm2 = copyin.mm2;
   constrain = copyin.constrain;
   ReadGroupCode=copyin.ReadGroupCode;
}

C_localpair::C_localpair(unsigned int pos1, int lm1,unsigned short anchor1, 
   unsigned short len11,unsigned short len21, char orient1, char q11,char q21, 
   char mm11,char mm21,char constrain1,
   unsigned int ReadGroupCode1) {
   pos = pos1;
   anchor = anchor1;
   lm = lm1;
   len1 = len11;
   len2 = len21;
   orient = orient1;
   q1 = q11;
   q2 = q21;
   mm1 = mm11;
   mm2 = mm21;
   constrain = constrain1;
   ReadGroupCode=ReadGroupCode1;
}

C_localpair::C_localpair(C_pairedread & pair1, char constrain1) {
  unsigned char a0 = pair1.read[0].align[0].anchor;
  unsigned char a1 = pair1.read[1].align[0].anchor;
  lm=0;
  anchor=0;
  len1=0;
  len2=0;
  orient='-';
  q1=0;
  q2=0;
  mm1=0;
  mm2=0;
  constrain = 0;
  if (a0!=a1) { return; }
  unsigned int p0 = pair1.read[0].align[0].pos;
  unsigned int p1 = pair1.read[1].align[0].pos;
  unsigned short l0 = pair1.read[0].align[0].len;
  unsigned short l1 = pair1.read[1].align[0].len;
  unsigned char o0 = pair1.read[0].align[0].sense;
  unsigned char o1 = pair1.read[1].align[0].sense;
  unsigned char Q0 = pair1.read[0].align[0].q;
  unsigned char Q1 = pair1.read[1].align[0].q;
  unsigned char MM0 = pair1.read[0].align[0].mm;
  unsigned char MM1 = pair1.read[1].align[0].mm;
  ReadGroupCode    = pair1.ReadGroupCode;
  anchor=a0;
  // harmonic mean q = (q0*q1)/(q0+q1);  ???
  // geometric mean seems to behave ok ???
  constrain = constrain1;
  if ((o0=='F')&&(o1=='R')) {
      lm = p1+l1-p0;
      pos=p0;
      len1=l0;
      len2=l1;
      q1=Q0;
      q2=Q1;
      mm1=MM0;
      mm2=MM1;
    } else if ((o0=='R')&&(o1=='F')) {
      lm = p0+l0-p1;
      pos=p1;
      len1=l1;
      len2=l0;
      q1=Q1;
      q2=Q0;
      mm1=MM1;
      mm2=MM0;
      if (constrain1==1) { constrain=2;}
      if (constrain1==2) { constrain=1;}
    } else {
      orient=(o0=='F'? '>' : '<');
      if (p0<p1) {
        lm = p1+l1-p0;
        pos=p0;
        len1=l0;
        len2=l1;
        q1=Q0;
        q2=Q1;
        mm1=MM0;
        mm2=MM1;
      } else {
        lm = p0+l0-p1;
        pos=p1;
        len1=l1;
        len2=l0;
        q1=Q1;
        q2=Q0;
        mm1=MM1;
        mm2=MM0;        
        if (constrain1==1) { constrain=2;}
        if (constrain1==2) { constrain=1;}
      };
  }
}


ostream &operator<<(ostream &output, const C_localpair & x)
{
   output << x.anchor << "\t" << x.pos << "\t"  << x.lm 
          << "\t"  << x.orient << "\t"  << int(x.q1) << "\t " << int(x.q2) 
          << "\t "<< int(x.constrain) << "\t " << x.ReadGroupCode;
   return output;
}

C_localpair& C_localpair::operator=(const C_localpair &rhs)
{
   this->pos = rhs.pos;
   this->anchor = rhs.anchor;
   this->lm = rhs.lm;
   this->len1 = rhs.len1;
   this->len2 = rhs.len2;
   this->orient = rhs.orient;
   this->q1 = rhs.q1;
   this->q2 = rhs.q2;
   this->mm1 = rhs.mm1;
   this->mm2 = rhs.mm2;
   this->constrain=rhs.constrain;
   this->ReadGroupCode=rhs.ReadGroupCode;
   return *this;
}

int C_localpair::operator==(const C_localpair &rhs) const
{
   if( this->pos != rhs.pos) return 0;
   if( this->anchor != rhs.anchor) return 0;
   if( this->lm != rhs.lm) return 0;
   if( this->orient != rhs.orient) return 0;
   if( this->len1 != rhs.len1) return 0;
   if( this->len2 != rhs.len2) return 0;
   if( this->q1 != rhs.q1) return 0;
   if( this->q2 != rhs.q2) return 0;
   if( this->mm1 != rhs.mm1) return 0;
   if( this->mm2 != rhs.mm2) return 0;
   if( this->constrain != rhs.constrain) return 0;
   if( this->ReadGroupCode != rhs.ReadGroupCode) return 0;
   return 1;
}

int C_localpair::operator<(const C_localpair &rhs) const
{
   if( this->anchor < rhs.anchor ) return 1;
   //if( this->ReadGroupCode < rhs.ReadGroupCode ) return 1;
   if( this->anchor == rhs.anchor && this->pos < rhs.pos ) return 1;
   if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->orient < rhs.orient ) return 1;
   if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->orient == rhs.orient  && this->lm < rhs.lm) return 1;
   if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->orient == rhs.orient  && this->lm == rhs.lm
       && this->len1 < rhs.len1 ) return 1;
   if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->orient == rhs.orient  && this->lm == rhs.lm
       && this->len1 == rhs.len1 && this->len2 < rhs.len2 ) return 1;
   if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->orient == rhs.orient  && this->lm == rhs.lm
       && this->len1 == rhs.len1 && this->len2 == rhs.len2 && this->q1 < rhs.q1) return 1;
   return 0;
}


//-----------------------------------------------------------------------------
// crosspair  (paired-read across two anchors)
//-----------------------------------------------------------------------------
C_crosspair::C_crosspair() {
    for (int i=0; i<2; i++) {
      read[i].pos=0; 
      read[i].len=0;
      read[i].anchor=0;
      read[i].sense=0;
      read[i].q=0;
      read[i].mm=0;
    }
    ReadGroupCode=0;
}

C_crosspair::C_crosspair(const C_crosspair &copyin)   // Copy constructor to handle pass by value.
{                             
   read[0] = copyin.read[0];
   read[1] = copyin.read[1];
   ReadGroupCode = copyin.ReadGroupCode;
}

C_crosspair::C_crosspair(const C_readmap  & read0, const C_readmap  & read1, const unsigned int ReadGroupCode1) {
   read[0] = read0;
   read[1] = read1;
   ReadGroupCode = ReadGroupCode1;
}

ostream &operator<<(ostream &output, const C_crosspair & x)
{
   output << x.read[0] << "\t" << x.read[1] <<  "\t "<< x.ReadGroupCode << endl;
   return output;
}

C_crosspair& C_crosspair::operator=(const C_crosspair &rhs)
{
   read[0] = rhs.read[0];
   read[1] = rhs.read[1];
   ReadGroupCode=rhs.ReadGroupCode;
   return *this;
}

int C_crosspair::operator==(const C_crosspair &rhs) //const
{
   if( !(this->read[0] == rhs.read[0])) return 0;
   if( !(this->read[1] == rhs.read[1])) return 0;
   if( !(this->ReadGroupCode== rhs.ReadGroupCode)) return 0;
   return 1;
}

int C_crosspair::operator<(const C_crosspair &rhs) const
{
   if( this->read[0] < rhs.read[0] ) return 1;
   if( this->read[0] == rhs.read[0] && this->read[1] < rhs.read[1] ) return 1;
   return 0;
}


//-----------------------------------------------------------------------------
// umpair  (paired-read with 0 unique and 1 multiply mapped)
//-----------------------------------------------------------------------------
C_umpair::C_umpair() {
    for (int i=0; i<2; i++) {
      read[i].pos=0; 
      read[i].len=0;
      read[i].anchor=0;
      read[i].sense=0;
      read[i].q=0;
      read[i].mm=0;
    }
    nmap=0;
    nmapA=0;
    elements=0;
    ReadGroupCode = 0;
}

C_umpair::C_umpair(const C_umpair &copyin)   // Copy constructor to handle pass by value.
{                             
   read[0] = copyin.read[0];
   read[1] = copyin.read[1];
   nmap=copyin.nmap;
   nmapA=copyin.nmapA;
   elements=copyin.elements;
   ReadGroupCode=copyin.ReadGroupCode;
}

C_umpair::C_umpair(const C_readmap  & read0, const C_readmap  & read1, 
   int nmap1, int elements1, unsigned int ReadGroupCode1) {
   read[0] = read0;
   read[1] = read1;
   nmap = nmap1;
   elements=elements1;
   ReadGroupCode=ReadGroupCode1;
}

//------------------------------------------------------------------------------
// create UM pair record
//------------------------------------------------------------------------------
C_umpair::C_umpair(const C_pairedread  & pair1, const C_anchorinfo & anchor1) {

  //----------------------------------------------------------------------------
  ReadGroupCode=pair1.ReadGroupCode;
  int N0 = pair1.read[0].align[0].nmap;     
  int N1 = pair1.read[1].align[0].nmap;     
  // check that both ends have maps and that at least one end is unique hit
  bool ok=  ( (N0*N1)>0 ) && ( (N0==1)||(N1==1) );
  if (!ok) {
     cerr << "bad UM fragment " << N0 << " " << N1 << endl;
     exit(-1);
  }
  
  // unique end, multiple-mapped end
  int eu=(N0==1? 0:1);
  // number of multiple maps
  int NM=(N0<N1? N1 : N0);  
  // markers for elements
  bool M0=(pair1.read[0].element)>0;
  bool M1=(pair1.read[1].element)>0;
  // override eu/em to force em to have element
  if (M0||M1) {
    eu=(M0? 1: 0);
  }
  int em=(eu==0? 1:0);
  //----------------------------------------------------------------------------
  // UU UM's are something of an oddity - but can happen from pair constraint 
  //----------------------------------------------------------------------------
  if (NM==1) {
    if (M0 && M1) {
      // UU UM with both ends in an element should not happen 
      cerr << "double element UU UM pair\n " << pair1.read[0].align[0] << endl;
      cerr <<  pair1.read[1].align[0] << endl;
    } else if ((!M0) && (!M1) ) {
      // UU UM with neither end in an element should not happen 
      cerr << "zero element UU UM pair\n " << pair1.read[0].align[0] << endl;
      cerr <<  pair1.read[1].align[0] << endl;
    }
  }
  // unique end
  read[0] = pair1.read[eu].align[0];
  // multiple map end
  read[1] = pair1.read[em].align[0];  
  // Total mumber of mappings at M end
  nmapA = NM;
  
  // loop over alignments to count *only* element hits for nmap
  nmap=nmapA;
  /*
  for (int i=0; i<N0 ; i++) { 
       if (anchor1.element[pair1.read[0].align[i].anchor]>0) { nmap++; } 
  }
  for (int i=0; i<N1 ; i++) { 
       if (anchor1.element[pair1.read[1].align[i].anchor]>0) { nmap++; } 
  }
  */
	
  // element bits at M end  - upper bits should be zero....
  elements = pair1.read[em].element;
  if (pair1.read[eu].element>0) {
    int elem1=pair1.read[eu].element;
    elements=elements|(elem1<<8);  // if U end has elements put them above 256
  }
}

ostream &operator<<(ostream &output, const C_umpair & x)
{
   output << x.read[0] << "\t" << x.read[1] << "\t "<< x.ReadGroupCode;
   output << "\t" << x.nmap << "\t" << int2binary(x.elements) << endl;
   return output;
}

C_umpair& C_umpair::operator=(const C_umpair &rhs)
{
   read[0] = rhs.read[0];
   read[1] = rhs.read[1];
   nmap = rhs.nmap;
   elements=rhs.elements;
   ReadGroupCode=rhs.ReadGroupCode;
   return *this;
}

int C_umpair::operator==(const C_umpair &rhs) //const
{
   if( !(this->read[0] == rhs.read[0])) return 0;
   if( !(this->read[1] == rhs.read[1])) return 0;
   // if( !(this->nmap == rhs.nmap)) return 0;
   //if( !(this->elements == rhs.elements)) return 0;
   return 1;
}

int C_umpair::operator<(const C_umpair &rhs) const
{
   if( this->read[0] < rhs.read[0] ) return 1;
   if( this->read[0] == rhs.read[0] && this->read[1] < rhs.read[1] ) return 1;
   if( this->read[0] == rhs.read[0] && this->read[1] == rhs.read[1] && this->nmap < rhs.nmap) return 1;
   return 0;
}

//------------------------------------------------------------------------------
// C_singleEnd class for single end reads - adds ReadGroupCode to readmap class
//------------------------------------------------------------------------------
C_singleEnd::C_singleEnd()   // Constructor
{
   /*
   Read.pos = 0;
   Read.anchor = 0;
   Read.len = 0;
   Read.sense = 0;
   Read.q=0;
   Read.mm=0;
   */
   ReadGroupCode=0; 
}

C_singleEnd::C_singleEnd(const C_singleEnd &copyin)   // Copy constructor to handle pass by value.
{                             
   pos = copyin.pos;
   anchor= copyin.anchor;
   len = copyin.len;
   sense = copyin.sense;
   q = copyin.q;
   mm=copyin.mm;
   ReadGroupCode = copyin.ReadGroupCode;
}

C_singleEnd::C_singleEnd(unsigned int pos1, unsigned short anchor1,unsigned short len1, 
   char sense1, char q1, char mm1, unsigned int ReadGroupCode1) {
   pos = pos1;
   anchor = anchor1;
   len = len1;
   sense = sense1;
   q = q1;
   mm = mm1;
   ReadGroupCode = ReadGroupCode1;
}

C_singleEnd&  C_singleEnd::operator=(const C_singleEnd &rhs) 
{
   pos = rhs.pos;
   anchor= rhs.anchor;
   len = rhs.len;
   sense = rhs.sense;
   q = rhs.q;
   mm = rhs.mm;
   ReadGroupCode = rhs.ReadGroupCode;
   return *this;
}

ostream &operator<<(ostream &output, const C_singleEnd & x)
{
   output << x.anchor << "\t" << x.pos << "\t" << x.len 
          << "\t" << x.sense <<"\t" << int(x.q) 
          <<"\t" << int(x.mm) << "\t "<< x.ReadGroupCode ;
   return output;
}

//------------------------------------------------------------------------------
// UM pair constraint test - 
//------------------------------------------------------------------------------
bool  C_umpair::constrain(int LMlow, int LMhigh) {
    // bug out if reads on different chromosomes
    if( !(this->read[0].anchor == this->read[1].anchor)) return false;
    int lm = -100000;
    unsigned int p0 = read[0].pos;
    unsigned int p1 = read[1].pos;
    unsigned short l0 = read[0].len;
    unsigned short l1 = read[1].len;
    unsigned char o0 = read[0].sense;
    unsigned char o1 = read[1].sense;
    if ((o0=='F')&&(o1=='R')) {
      lm = p1+l1-p0;
    } else if ((o0=='R')&&(o1=='F')) {
      lm = p0+l0-p1;
    }
    if (lm>LMhigh) return false;
    if (lm<LMlow) return false;
    return true;
}
 
 
//-----------------------------------------------------------------------------
// depth of coverage class
//-----------------------------------------------------------------------------
C_depth::C_depth() {                            // constructor
  C_depth d1(0);
  *this = d1;
}

//-----------------------------------------------------------------------------
// depth of coverage class
//-----------------------------------------------------------------------------
C_depth::C_depth(int L1) {                            // constructor
  n.resize(L1,0); //=d1;
  pos0=0;
  pos1=L1-1;
  light=1;
  nzbins=0;
  nzMedian=0.0;   
}
void C_depth::calcStats() {
  calcStats(5001,-0.5,5000.5);
}

void C_depth::calcStats(int Nbin, float x1, float x2) {
  this->Stats.Initialize(Nbin,x1,x2);  
  this->Stats.h.setTitle("RD count of reads/base ");
  this->Stats.h.setXlabel("bases");
  int L = n.size();
  for (int p=0; p<L; p++) {
    this->Stats.Fill1(this->n[p]);  
  }
  this->Stats.Finalize();     
  L = this->Stats.h.Nbin;
  nzMedian=0;
  double cumu=0.0;
  nzbins=int(Stats.h.Ntot-Stats.h.Nunder-Stats.h.n[0]);
  for (int b=1; b<this->Stats.h.Nbin; b++) {
    cumu += int(Stats.h.n[b]);
    if (cumu> (nzbins/2) ) {
       nzMedian=Stats.h.xc[b];
       break;
    }
  }
}


ostream &operator<<(ostream &output, C_depth & d1)
{
   output << d1.name  << endl;
   output << d1.Stats  << endl;
   return output;
}

//-----------------------------------------------------------------------------
// anchor marking class
//-----------------------------------------------------------------------------
C_marker::C_marker(int L1) {                            // constructor
  x.resize(L1,false); //=x1;
}

//-----------------------------------------------------------------------------
// add fasta file repeat annotations to marking vector x
//-----------------------------------------------------------------------------
int C_depth::addMarks(const string & fastafile, const string & cname1, char X)
{
  if (n.size()<1) {
     return 0;
  }
  if (fastafile.length()==0) {
     return 0;
  }
  // boost regex matches
  string match,match1,match2;
  // Patterns to match
  string patternFastaHeader("^>(\\s*\\S+\\s*.*)$");
  string patternFastaName("^\\s*(\\S+)");
  string patternContigName(cname1);
  // target seq flag
  bool  thisSeq = false;
  // base  counter
  int baseCount = 0;
  int Nadded = 0;
  // sequence name
  string seqName;
  // sequence header
  string seqHeader;
  // input line
  string line;
  //----------------------------------------------------------------------------
  // open input FASTA DNA file
  //----------------------------------------------------------------------------
  ifstream dnaIn(fastafile.c_str(), ios::in);  
  if (!dnaIn) {
    cerr << "Unable to open fasta file: " << fastafile << endl;
    exit(1);
  }  
  while (getline(dnaIn, line)) {
    // header line (long format): register previous sequence and start new
    if ( RE2::FullMatch(line.c_str(),patternFastaHeader.c_str(),&match) ) {
      // retreive info for new sequence
      seqHeader = match;      
      // parse out sequence name
      seqName = seqHeader;
      if ( RE2::FullMatch(seqHeader.c_str(),patternFastaName.c_str(),&match1) ) {     
          seqName = match1;
      }
      thisSeq = RE2::FullMatch(seqName.c_str(),patternContigName.c_str(),&match2);
      baseCount = 0;
    }  else if (thisSeq) {
    
      size_t found=line.find_first_of(X);
      while (found!=string::npos)
      {
        int p = found + baseCount;
        if (p>=int(n.size()-1)) {
          cout << " fasta file " << fastafile << " contig " <<  cname1;
          cout << " exceeds length of contig " << n.size() << " oops " <<  endl;
          return Nadded;
        }
        n[p] = float(-1e10);
        Nadded++;
        found=line.find_first_of(X,found+1);
      }
      baseCount += line.size();
    }
  }  
  dnaIn.close();
  cout << " added " << Nadded << " marked postions from fasta file " << fastafile << " contig " <<  cname1;    
  return Nadded;
}
//-----------------------------------------------------------------------------
// calculate GC content fraction from reference fasta file 
//-----------------------------------------------------------------------------
int C_depth::GCcontent(const string & fastafile, const string & cname1)
{
  if (n.size()<1) {
     return 0;
  }
  if (fastafile.length()==0) {
     return 0;
  }
  // boost regex matches
  string match,match1,match2;
  // Patterns to match
  string patternFastaHeader("^>(\\s*\\S+\\s*.*)$");
  string patternFastaName("^\\s*(\\S+)");
  string patternContigName(cname1);
  // target seq flag
  bool  thisSeq = false;
  // base  counter
  int baseCount = 0;
  int Ngc = 0;
  // sequence name
  string seqName;
  // sequence header
  string seqHeader;
  // input line
  string line;
  light_t info;
  info.i64=light;
  int binsize=info.i16[0];
  //----------------------------------------------------------------------------
  // open input FASTA DNA file
  //----------------------------------------------------------------------------
  ifstream dnaIn(fastafile.c_str(), ios::in);  
  if (!dnaIn) {
    cerr << "Unable to open fasta file: " << fastafile << endl;
    exit(1);
  }  
  while (getline(dnaIn, line)) {
    // header line (long format): register previous sequence and start new
    //if (boost::regex_search(line, match, patternFastaHeader)) {
    char C1 = line[0];
    if (C1=='>') { 
      if (! RE2::FullMatch(line.c_str(),patternFastaHeader.c_str(),&match)) {
        cerr << "Bad Fasta Header for GC content: " << line << endl;
        exit(1);
      }
      // retreive info for new sequence
      seqHeader = match;      
      // parse out sequence name
      seqName = seqHeader;
      if (RE2::FullMatch(seqHeader.c_str(),patternFastaName.c_str(),&match1) ) {     
          seqName = match1;
      }
      thisSeq = RE2::FullMatch(seqName.c_str(),patternContigName.c_str(),&match2) ;
      baseCount = 0;
      
    }  else if (thisSeq) {
       size_t Nb=line.size();
       for (int q=0; q<int(Nb); q++) {
         int p = (q + baseCount)/binsize;
         if (p>=int(n.size()-1)) {
           cout << " fasta file " << fastafile << " contig " <<  cname1;
           cout << " bin exceeds length of contig " << n.size()  <<  endl;
           dnaIn.close();
           cout << " GC content: " << float(Ngc)/baseCount << endl;
           return Ngc;
         }
         float isgc =0;
         char base = toupper(line[q]);
         if ( (base=='G') or (base=='C') ) isgc=1.0;  
         n[p] += isgc/binsize;
         Ngc+=int(isgc);
       }
       baseCount += line.size();
    }
  }  
  dnaIn.close();
  cout << " GC content: " << float(Ngc)/baseCount << " fasta file " << fastafile << " contig " <<  cname1;    
  return Ngc;
}

//-----------------------------------------------------------------------------
// add fasta file repeat annotations to marking vector x
//-----------------------------------------------------------------------------
int C_marker::addMarks(const string & fastafile, const string & cname1, char X)
{
  if (x.size()<1) {
     return 0;
  }
  if (fastafile.length()==0) {
     return 0;
  }
  // boost regex matches
  string match,match1,match2;
  // Patterns to match
  string patternFastaHeader("^>(\\s*\\S+\\s*.*)$");
  string patternFastaName("^\\s*(\\S+)");
  string patternContigName(cname1);
  // target seq flag
  bool  thisSeq = false;
  // base  counter
  int baseCount = 0;
  int Nadded = 0;
  // sequence name
  string seqName;
  // sequence header
  string seqHeader;
  // input line
  string line;
  //----------------------------------------------------------------------------
  // open input FASTA DNA file
  //----------------------------------------------------------------------------
  ifstream dnaIn(fastafile.c_str(), ios::in);  
  if (!dnaIn) {
    cerr << "Unable to open fasta file: " << fastafile << endl;
    exit(1);
  }  
  while (getline(dnaIn, line)) {
    // header line (long format): register previous sequence and start new
    //if (   boost::regex_search(line, match, patternFastaHeader)) {
	if ( RE2::FullMatch(line.c_str(),patternFastaHeader.c_str(),&match) ) {
			// retreive info for new sequence
      seqHeader = match;      
      // parse out sequence name
      seqName = seqHeader;
      if ( RE2::FullMatch(seqHeader.c_str(),patternFastaName.c_str(),&match1) ) {     
          seqName = match1;
      }
      thisSeq = RE2::FullMatch(seqName.c_str(),patternContigName.c_str(),&match2)  ;
      baseCount = 0;
    }  else if (thisSeq) {
    
      size_t found=line.find_first_of(X);
      while (found!=string::npos)
      {
        int p = found + baseCount;
        if (p>=int(x.size()-1)) {
          cout << " fasta file " << fastafile << " contig " <<  cname1;
          cout << " exceeds length of contig " << x.size() << " oops " <<  endl;
          return Nadded;
        }
        x[p] = true;
        Nadded++;
        found=line.find_first_of(X,found+1);
      }
      baseCount += line.size();
    }
  }  
  dnaIn.close();
  cout << " added " << Nadded << " marked postions from fasta file " << fastafile << " contig " <<  cname1;    
  return Nadded;
}

//=============================================
// contig class (container for all PE mappings)
//=============================================
C_contig::C_contig(string & contigName1, int L1, bool doRepeatCheck) {
  contigName = contigName1;     // contig name
  Length = L1;                  // contig length;
  uniquified = 0;               // reset flag
  //depth.n.resize(L1,0);                         
  //starts.n.resize(L1,0);
  if ((L1>0)&&(doRepeatCheck)) {
    float zip=0;           
    repeat.n.resize(L1,zip);    // non-unique read mapping positions
  }
}

// fetch set name
string C_contig::getContigName() const
{
  return this->contigName;
}

// set set name vector
void C_contig::setContigName(string & s1)
{
  this->contigName = s1;
}

void C_contig::calcStats() {
     // init total count of unique reads 
     totalUniqueReads=0;   
     // init total count of repeat bases
     totalRepeatBases=0;  
     // init total count of telomeric non-accessable bases
     totalNoCovBases=0;     
     //-------------------------------------------------------------------------
     // Don't waste time making histogram of empty list...
     //-------------------------------------------------------------------------
     unsigned int L = repeat.n.size();  
     if (repeat.Stats.N>0) {
       this->repeat.Stats.Initialize(101,-0.5,1000.5);  
       repeat.Stats.h.setTitle("NA count of bases with repeats");
       repeat.Stats.h.setXlabel("rpt");
       for (unsigned int i=0; i<L; i++) {
          this->repeat.Stats.Fill1(this->repeat.n[i]);  
          // bump repeat base count
          totalRepeatBases+=(repeat.n[i]>0);
          // check for first base with a repeat
          if ((totalNoCovBases<1.0)&&(repeat.n[i]>0)) {
             totalNoCovBases=double(i)-1.0;
          }
       }
       this->repeat.Stats.Finalize();    
     } else {
       this->repeat.Stats.Initialize(2,-0.5,1.5);  
       repeat.Stats.h.setTitle("NA count of bases with repeats");
       repeat.Stats.h.setXlabel("rpt");
       this->repeat.Stats.Finalize();    
     }
     // check for last base with a repeat
     // --i should start at L-1, but it doesnt...
     
     if (L>0) { 
       for (unsigned int i=(L-1); i>0; --i) {
          totalNoCovBases++;
          if (repeat.n[i]>0) {
            break;
          }
       }
     } 
     
     this->pairStats.Initialize(10001,-0.5,10000.5);       
     pairStats.h.setTitle("LF Fragment paired-read mapping length");
     pairStats.h.setXlabel("LM");
     unsigned int N = localpairs.size();
     list<C_localpair>::iterator i;
     for(i=localpairs.begin(); i != localpairs.end(); ++i) {
        this->pairStats.Fill1((*i).lm);  
     }
     pairStats.Finalize();
     //
     totalUniqueReads+=2*localpairs.size();          

     int Na = anchors.L.size();
     //this->crossStats.Initialize(101,-0.5,100.5);       
     this->crossStats.Initialize(Na,-0.5,double(Na)-0.5);       
     crossStats.h.setBinLabels(anchors.names);
     crossStats.h.setXlabel("anchor");
     string title = "CA Cross pair "+contigName+" linked to anchors ";
     crossStats.h.setTitle(title);

     N = crosspairs.size();
     list<C_crosspair>::iterator j;
     for(j=crosspairs.begin(); j != crosspairs.end(); ++j) {
        //this->crossStats.Fill1((*j).read[0].anchor);  
        this->crossStats.Fill1((*j).read[1].anchor);  
     }
     crossStats.Finalize();
     //
     totalUniqueReads+=crosspairs.size();          
     totalUniqueReads+=dangle.size();          
     totalUniqueReads+=umpairs.size();          
     totalUniqueReads+=singleton.size();               
}

//===========================================================
// extract read length and contig length info from span files
//===========================================================
void C_contig::calcLengths() {
   //  Length = repeat.n.size(); 
     double LR0 = 0; 
     double LR1 = 0; 
     double LR2 = 0; 
     list<C_localpair>::iterator i;
     for(i=localpairs.begin(); i != localpairs.end(); ++i) {
        LR0 = LR0+1.0;
        LR1 = LR1+(*i).len1;  
        LR2 = LR2+(*i).len1*(*i).len1;  
     }
     aLR = LR1 / LR0;
     sLR = sqrt( ( LR2 / LR0 )-aLR*aLR) ;
}

int C_contig::setLengthFromAnchor() {      // set contig Length from anchor[contigName].L
     
     //  Length = repeat.n.size(); 
     
     size_t Np = localpairs.size();
     int La = anchors.L[contigName];
     if (Np>0) {
        Length=La;
     } else {
        Length=0;
     }
     
     // cout << contigName << ":\t " << Length << endl;
     
     return Length;
}


void C_contig::calcDepth(int Qmin) {
  // initialize depth arrays (this can take a while)
  read_depth.n.resize(Length,0);
  read_start.n.resize(Length,0);
  frag_depth.n.resize(Length,0);
  //
  // mapping quality minumum?  
  // int Qmin=pars.getQmin();
  //
  int p0=0;
  int p1=0;
  // fill depth from local pairs
  unsigned int N = localpairs.size();
  list<C_localpair>::iterator i;
  // loop over pairs
  for(i=localpairs.begin(); i != localpairs.end(); ++i) {
    //loop over ends for read depth
    for (int e=0; e<2; e++) {
      // skip ends that are constrained
      if (((e+1)&&(*i).constrain)>0) {continue;} 
      // 
      if ((e==0)&&((*i).q1<Qmin)) {continue;} 
      if ((e==1)&&((*i).q2<Qmin)) {continue;} 
      //
      switch (e) {
        case 0: 
          p0 = (*i).pos;
          p1 = p0+(*i).len1;
          break;
        case 1:
          p0 = (*i).pos+(*i).lm-(*i).len2;
          p1 = (*i).pos+(*i).lm;
      }
      read_start.n[p0]+=1;
      for (int p = p0; p<p1; p++) {
        if ((p<0)|(p>Length)) { 
          cerr << " bound problem in calcDepth pairs" << p << endl;
        } else {
          read_depth.n[p]+=1;
        }
      }
    }
    // fragment depth
    int lm = (*i).lm;
    if (lm>0) {
      p0 = (*i).pos;
      p1 = p0+(*i).lm;
    } else {
      p1 = (*i).pos;
      p0 = p1+(*i).lm;
    }
    for (int p = p0; p<p1; p++) {
      if ((p<0)|(p>Length)) { 
        cerr << " bound problem in calcDepth frag depth from pairs" << p << endl;
      } else {
        frag_depth.n[p]+=1;
      }
    }    
  }
  // fill depth from cross pairs (not so many...)
  N = crosspairs.size();
  list<C_crosspair>::iterator i1;
  // loop over pairs
  for(i1=crosspairs.begin(); i1 != crosspairs.end(); ++i1) {
    //first end for read depth
    if ((*i1).read[0].q<Qmin) {continue;} 
    p0 = (*i1).read[0].pos;
    p1 = p0+(*i1).read[0].len;
    read_start.n[p0]+=1;
    for (int p = p0; p<p1; p++) {
      if ((p<0)|(p>Length)) { 
        cerr << " bound problem in calcDepth cross" << p << endl;
      } else {
        read_depth.n[p]+=1;
      }
    }
  }
  // fill depth from dangle starts 
  N = dangle.size();
  list<C_singleEnd>::iterator i2; 
  for(i2=dangle.begin(); i2 != dangle.end(); ++i2) {
    //first end for read depth
    if ((*i2).q<Qmin) {continue;} 
    p0 = (*i2).pos;
    p1 = p0+(*i2).len;
    read_start.n[p0]+=1;
    for (int p = p0; p<p1; p++) {
      if ((p<0)|(p>Length)) { 
        cerr << " bound problem in calcDepth dangle" << p << endl;
      } else {
        read_depth.n[p]+=1;
      }
    }
  }
  // fill depth from umpairs 
  N = umpairs.size();
  list<C_umpair>::iterator i3;
  for(i3=umpairs.begin(); i3 != umpairs.end(); ++i3) {
    if ((*i3).read[0].q<Qmin) {continue;} 
    //first end for read depth
    p0 = (*i3).read[0].pos;
    p1 = p0+(*i3).read[0].len;
    read_start.n[p0]+=1;
    for (int p = p0; p<p1; p++) {
      if ((p<0)|(p>Length)) { 
        cerr << " bound problem in calcDepth umpairs " << p << endl;
      } else {
        read_depth.n[p]+=1;
      }
    }
  }
  // fill depth from unique ends 
  N = singleton.size();
  for(i2=singleton.begin(); i2 != singleton.end(); ++i2) {
    if ((*i2).q<Qmin) {continue;} 
    //first end for read depth
    p0 = (*i2).pos;
    p1 = p0+(*i2).len;
    read_start.n[p0]+=1;
    for (int p = p0; p<p1; p++) {
      if ((p<0)|(p>Length)) { 
        cerr << " bound problem in calcDepth singleton" << p << endl;
      } else {
        read_depth.n[p]+=1;
      }
    }
  }
  this->read_depth.Stats.Initialize(5001,-0.5,5000.5);  
  read_depth.Stats.h.setTitle("RD count of reads/base ");
  read_depth.Stats.h.setXlabel("bases");
  this->read_start.Stats.Initialize(5001,-0.5,5000.5);  
  read_start.Stats.h.setTitle("SD count of starts/base ");
  read_start.Stats.h.setXlabel("bases");
  this->frag_depth.Stats.Initialize(5001,-0.5,5000.5);  
  frag_depth.Stats.h.setTitle("FD count of fragments/base ");
  frag_depth.Stats.h.setXlabel("bases");
  for (int p=0; p<Length; p++) {
  //   if (this->repeat.n[i]==0) {
        this->read_depth.Stats.Fill1(this->read_depth.n[p]);  
        this->read_start.Stats.Fill1(this->read_start.n[p]);  
        this->frag_depth.Stats.Fill1(this->frag_depth.n[p]);  
 //    }
  }
  this->read_depth.Stats.Finalize();     
  this->read_start.Stats.Finalize();     
  this->frag_depth.Stats.Finalize();     
}

void C_contig::calcDepth(int P0, int P1, int Qmin) {
  // initialize depth arrays (this can take a while)
  if (P1>=Length) P1=Length-1;
  int L = P1-P0;
  read_depth.pos0=P0;
  read_depth.pos1=P1;
  read_start.pos0=P0;
  read_start.pos1=P1;
  frag_depth.pos0=P0;
  frag_depth.pos1=P1;
  read_depth.n.resize(L,0);
  read_start.n.resize(L,0);
  frag_depth.n.resize(L,0);
  int pos0=P0;
  int pos1=P1;
  int p0=0;
  int p1=0;
  // fill depth from local pairs
  unsigned int N = localpairs.size();
  list<C_localpair>::iterator i;
  // loop over pairs
  for(i=localpairs.begin(); i != localpairs.end(); ++i) {
    //loop over ends for read depth
    for (int e=0; e<2; e++) {
      switch (e) {
        case 0: 
          if ((*i).q1<Qmin) {continue;} 
          p0 = (*i).pos;
          p1 = p0+(*i).len1;
          break;
        case 1:
          if ((*i).q2<Qmin) {continue;} 
          p0 = (*i).pos+(*i).lm-(*i).len2;
          p1 = (*i).pos+(*i).lm;
      }
      if ((p0>=pos0)&&(p0<=pos1)) {read_start.n[p0-pos0]+=1;}
      if (!boundlimit(p0,p1,pos0,pos1)) {continue;}
      for (int p = p0; p<p1; p++) {
          read_depth.n[p]+=1;
      }
    }
    // fragment depth
    int lm = (*i).lm;
    if (lm>0) {
      p0 = (*i).pos;
      p1 = p0+(*i).lm;
    } else {
      p1 = (*i).pos;
      p0 = p1+(*i).lm;
    }
    if (!boundlimit(p0,p1,pos0,pos1)) {continue;}
    for (int p = p0; p<p1; p++) {
        frag_depth.n[p]+=1;
    }    
  }
  // fill depth from cross pairs (not so many...)
  N = crosspairs.size();
  list<C_crosspair>::iterator i1;
  // loop over pairs
  for(i1=crosspairs.begin(); i1 != crosspairs.end(); ++i1) {
    if ((*i1).read[0].q<Qmin) {continue;} 
    //first end for read depth
    p0 = (*i1).read[0].pos;
    p1 = p0+(*i1).read[0].len;
    if ((p0>=pos0)&&(p0<=pos1)) {read_start.n[p0-pos0]+=1;}
    if (!boundlimit(p0,p1,pos0,pos1)) {continue;}
    for (int p = p0; p<p1; p++) {
        read_depth.n[p]+=1;
    }
  }
  // fill depth from dangles 
  N = dangle.size();
  list<C_singleEnd>::iterator i2; 
  for(i2=dangle.begin(); i2 != dangle.end(); ++i2) {
    if ((*i2).q<Qmin) {continue;} 
    //first end for read depth
    p0 = (*i2).pos;
    p1 = p0+(*i2).len;
    if ((p0>=pos0)&&(p0<=pos1)) {read_start.n[p0-pos0]+=1;}
    if (!boundlimit(p0,p1,pos0,pos1)) {continue;}
    for (int p = p0; p<p1; p++) {
        read_depth.n[p]+=1;
    }
  }
  // fill depth from unique ends 
  N = umpairs.size();
  list<C_umpair>::iterator i3;
  for(i3=umpairs.begin(); i3 != umpairs.end(); ++i3) {
    if ((*i3).read[0].q<Qmin) {continue;} 
    //first end for read depth
    p0 = (*i3).read[0].pos;
    p1 = p0+(*i3).read[0].len;
    if ((p0>=pos0)&&(p0<=pos1)) {read_start.n[p0-pos0]+=1;}
    if (!boundlimit(p0,p1,pos0,pos1)) {continue;}
    for (int p = p0; p<p1; p++) {
        read_depth.n[p]+=1;
    }
  }
  // fill depth from unique ends 
  N = singleton.size();
  for(i2=singleton.begin(); i2 != singleton.end(); ++i2) {
    if ((*i2).q<Qmin) {continue;} 
    //first end for read depth
    p0 = (*i2).pos;
    p1 = p0+(*i2).len;
    if ((p0>=pos0)&&(p0<=pos1)) {read_start.n[p0-pos0]+=1;}
    if (!boundlimit(p0,p1,pos0,pos1)) {continue;}
    for (int p = p0; p<p1; p++) {
        read_depth.n[p]+=1;
    }
  }
  this->read_depth.Stats.Initialize(5001,-0.5,5000.5);  
  read_depth.Stats.h.setTitle("RD count of reads/base ");
  read_depth.Stats.h.setXlabel("bases");
  this->read_start.Stats.Initialize(5001,-0.5,5000.5);  
  read_start.Stats.h.setTitle("SD count of starts/base ");
  read_start.Stats.h.setXlabel("bases");
  this->frag_depth.Stats.Initialize(5001,-0.5,5000.5);  
  frag_depth.Stats.h.setTitle("FD count of fragments/base ");
  frag_depth.Stats.h.setXlabel("bases");
  for (int p=0; p<L; p++) {
 //    if (!this->repeat.x[i]) {
        this->read_depth.Stats.Fill1(this->read_depth.n[p]);  
        this->read_start.Stats.Fill1(this->read_start.n[p]);  
        this->frag_depth.Stats.Fill1(this->frag_depth.n[p]);  
 //    }
  }
  this->read_depth.Stats.Finalize();     
  this->read_start.Stats.Finalize();     
  this->frag_depth.Stats.Finalize();     
}

//------------------------------------------------------------------------------
// Calculate UU Fragment coverage for this contig
//------------------------------------------------------------------------------
void C_contig::calcFragDepth(int lmLow, int lmHigh) {
  // initialize depth arrays (this can take a while)
  calcFragDepth(0,Length,lmLow,lmHigh);
}

void C_contig::calcFragDepth(int P0, int P1, int lmLow, int lmHigh) {
  // initialize depth arrays (this can take a while)
  if (P1>=Length) P1=Length-1;
  int L = P1-P0;
  frag_depth.n.resize(L,0);
  frag_depth.pos0=P0;
  frag_depth.pos1=P1;
  //int lmLow = int(pars.getFragmentLengthLo());
  //int lmHigh = int(pars.getFragmentLengthHi());    
  int pos0=P0;
  int pos1=P1;
  int p0=0;
  int p1=0;
  // fill depth from local pairs
  //unsigned int N = localpairs.size();
  list<C_localpair>::iterator i;
  // loop over pairs
  for(i=localpairs.begin(); i != localpairs.end(); ++i) {
    // fragment depth
    int lm = (*i).lm;
    // demand usual fragment
    if ((lm<lmLow)|(lm>lmHigh)) {continue; }
    // add only the non-read part of the fragment
    p0 = (*i).pos+(*i).len1;
    p1 = p0+(*i).lm-(*i).len2;
    if (!boundlimit(p0,p1,pos0,pos1)) {continue;}
    for (int p = p0; p<p1; p++) {
        frag_depth.n[p]+=1;
    }    
  }
  this->frag_depth.Stats.Initialize(5001,-0.5,5000.5);  
  frag_depth.Stats.h.setTitle("FD count of fragments/base ");
  frag_depth.Stats.h.setXlabel("bases");
  for (int p=0; p<L; p++) {
        this->frag_depth.Stats.Fill1(this->frag_depth.n[p]);  
  }
  this->frag_depth.Stats.Finalize();     
}




void C_contig::calcStarts(int Qmin) {
  // check if read starts already done 
  if (read_start.Stats.N>0) return;
  // initialize depth arrays (this can take a while)
  read_start.n.resize(Length,0);
  int p0=0;
  // fill depth from local pairs
  unsigned int N = localpairs.size();
  list<C_localpair>::iterator i;
  for(i=localpairs.begin(); i != localpairs.end(); ++i) {
    //loop over ends for read depth
    for (int e=0; e<2; e++) {
      switch (e) {
        case 0: 
          if ((*i).q1<Qmin) {continue;} 
          p0 = (*i).pos;
          break;
        case 1:
          if ((*i).q1<Qmin) {continue;} 
          p0 = (*i).pos+(*i).lm-(*i).len2;
      }
      read_start.n[p0]+=1;
    }
  }
  // fill depth from cross pairs (not so many...)
  N = crosspairs.size();
  list<C_crosspair>::iterator i1;
  for(i1=crosspairs.begin(); i1 != crosspairs.end(); ++i1) {
    if ((*i1).read[0].q<Qmin) {continue;} 
    p0 = (*i1).read[0].pos;
    read_start.n[p0]+=1;
  }
  // fill depth from dangles
  N = dangle.size();
  list<C_singleEnd>::iterator i2; 
  for(i2=dangle.begin(); i2 != dangle.end(); ++i2) {
    if ((*i2).q<Qmin) {continue;} 
    p0 = (*i2).pos;
    read_start.n[p0]+=1;
  }
  // fill depth from umpairs
  N = umpairs.size();
  list<C_umpair>::iterator i3;
  for(i3=umpairs.begin(); i3 != umpairs.end(); ++i3) {
    if ((*i3).read[0].q<Qmin) {continue;} 
    p0 = (*i3).read[0].pos;
    read_start.n[p0]+=1;
  }
  // fill depth from unique ends 
  N = singleton.size();
  for(i2=singleton.begin(); i2 != singleton.end(); ++i2) {
    if ((*i2).q<Qmin) {continue;} 
    p0 = (*i2).pos;
    read_start.n[p0]+=1;
  }
  this->read_start.Stats.Initialize(5001,-0.5,5000.5);  
  read_start.Stats.h.setTitle("SD count of starts/base ");
  read_start.Stats.h.setXlabel("bases");
  for (int p=0; p<Length; p++) {
     //if (!this->repeat.x[i]) {
        this->read_start.Stats.Fill1(this->read_start.n[p]);  
     //} else {
     //   this->read_start.Stats.Fill1(-1);  
     //}
  }
  this->read_start.Stats.Finalize();     
}

C_depth  C_contig::countReads(int binsize) {
  // check if read starts already done 
  if (read_start.Stats.N>0) {
    C_depth c0;
    return c0;
  }
  int nbin = 1+(Length/binsize);
  // initialize depth arrays (this can take a while)
  C_depth c(nbin);
  // fill header info  
  light_t info;
  info.i64=0;  
  info.i16[0]=binsize;
  info.i16[1]=0;  
  c.light = info.i64;
  // position var
  int p0=0;
  // fill depth from local pairs
  unsigned int N = localpairs.size();
  list<C_localpair>::iterator i;
  for(i=localpairs.begin(); i != localpairs.end(); ++i) {
    //loop over ends for read depth
    for (int e=0; e<2; e++) {
      // skip an end if was found by constraint  
	  // if (((e+1)&((*i).constrain)>0) {continue;} ???? 
	  if (((e+1)&(*i).constrain)>0) {continue;} 
      switch (e) {
        case 0: 
          p0 = (*i).pos;
          break;
        case 1:
          p0 = (*i).pos+(*i).lm-(*i).len2;
      }
      p0=p0/binsize;
      c.n[p0]+=1;
    }
  }
  // fill depth from cross pairs (not so many...)
  N = crosspairs.size();
  list<C_crosspair>::iterator i1;
  for(i1=crosspairs.begin(); i1 != crosspairs.end(); ++i1) {
    p0 = (*i1).read[0].pos/binsize;
    c.n[p0]+=1;
  }
  // fill depth from dangles
  N = dangle.size();
  list<C_singleEnd>::iterator i2; 
  for(i2=dangle.begin(); i2 != dangle.end(); ++i2) {
    p0 = (*i2).pos/binsize;
    c.n[p0]+=1;
  }
  // fill depth from unique ends 
  N = umpairs.size();
  list<C_umpair>::iterator i3;
  for(i3=umpairs.begin(); i3 != umpairs.end(); ++i3) {
    p0 = (*i3).read[0].pos/binsize;
    c.n[p0]+=1;
  }
  // fill depth from singleton ends 
  N = singleton.size();
  for(i2=singleton.begin(); i2 != singleton.end(); ++i2) {
    p0 = (*i2).pos/binsize;
    c.n[p0]+=1;
  }
  // stats
  c.calcStats();
  /*
  c.Stats.Initialize(15001,-0.5,15000.5);  
  c.Stats.h.setTitle("UR count of starts/base ");
  c.Stats.h.setXlabel("bases");
  for (int p=0; p<nbin; p++) {
     c.Stats.Fill1(c.n[p]);  
  }
  c.Stats.Finalize();     
  */
  return c;
}


unsigned short C_contig::getAnchorIndex() {
// this seems to alaways return zero
/*
  list<C_localpair>::iterator i;
  i=localpairs.begin(); 
  unsigned short a = (*i).anchor;
*/
  unsigned short a=0;
  for (size_t i = 0; i<this->anchors.names.size(); i++) {
      string name= this->anchors.names[i];
      a++;
      if (contigName.compare(name)==0) break;
  }
  return a;
}

long long C_contig::countReads(int P0, int P1) {
  // initialize depth arrays (this can take a while)
  if (P1>=Length) P1=Length-1;
  if (P0<0) P0=0;
  int p0, p1;
  long long Nread = 0;
  // fill depth from local pairs
  unsigned int N = localpairs.size();
  list<C_localpair>::iterator i;
  // 
  cout << " countReads " << contigName << " " << P0 << " " << P1 << " " << N << endl; 
  // loop over pairs
  for(i=localpairs.begin(); i != localpairs.end(); ++i) {
    //loop over ends for read depth
    for (int e=0; e<2; e++) {
      switch (e) {
        case 0: 
          p0 = (*i).pos;
          p1 = p0+(*i).len1;
          break;
        case 1:
          p0 = (*i).pos+(*i).lm-(*i).len2;
          p1 = (*i).pos+(*i).lm;
      }
      if ((p0>=P0)&&(p0<=P1)) {
       Nread+=1;
      }
    }
  }
  //
  cout << " countReads pairs " << Nread << endl;
    
  // fill depth from cross pairs (not so many...)
  N = crosspairs.size();
  list<C_crosspair>::iterator i1;
  // loop over pairs
  for(i1=crosspairs.begin(); i1 != crosspairs.end(); ++i1) {
    //first end for read depth
    p0 = (*i1).read[0].pos;
    p1 = p0+(*i1).read[0].len;
    if ((p0>=P0)&&(p0<=P1)) {
      Nread+=1;
    } else if (p0>P1) {
      break;
    }
  }
  // fill depth from dangles
  N = dangle.size();
  list<C_singleEnd>::iterator i2; 
  for(i2=dangle.begin(); i2 != dangle.end(); ++i2) {
    //first end for read depth
    p0 = (*i2).pos;
    if ((p0>=P0)&&(p0<=P1)) {
      Nread+=1;
    } else if (p0>P1) {
      break;
    }
  }
  // fill depth from unique ends 
  N = umpairs.size();
  list<C_umpair>::iterator i3;
  for(i3=umpairs.begin(); i3 != umpairs.end(); ++i3) {
    //first end for read depth
    p0 = (*i3).read[0].pos;
    p1 = p0+(*i3).read[0].len;
    if ((p0>=P0)&&(p0<=P1)) {
      Nread+=1;
    } else if (p0>P1) {
      break;
    }
  }
  // fill depth from unique ends 
  N = singleton.size();
  for(i2=singleton.begin(); i2 != singleton.end(); ++i2) {
    //first end for read depth
    p0 = (*i2).pos;
    p1 = p0+(*i2).len;
    if ((p0>=P0)&&(p0<=P1)) {
      Nread+=1;
    } else if (p0>P1) {
      break;
    }
  }
  return Nread;
}


int C_contig::howFar(int p0, int Nnu) {
  //---------------------------------------------------------------------------
  // returns end position of region starting at P0 
  // with Nnu non-repeat bases
  //---------------------------------------------------------------------------
  if (p0<0) p0=0;
  int inu = 0;
  int n;
  if (repeat.n.size()<1)  {
     n=p0+Nnu;
     n = (n>=Length? Length-1: n);
     return n;
  }
  // loop over repeats until reaching Nnu non-repeats
  for(n=p0; n < (Length-1); n++) {
      inu+=1-int(repeat.n[n]>0);
      if (inu == Nnu) { break; }
  }
  return n;
}

void C_contig::sort() {

    // check if already done, bug out if done
    if (uniquified>0) return;

    localpairs.sort();               
    crosspairs.sort();               
    dangle.sort();         
    singleton.sort();
    umpairs.sort();    
    uniquified=1;
    
}

void C_contig::uniquify() {

    // check if already done, bug out if done
    if (uniquified>1) return;    

    int N0 = localpairs.size();               
    localpairs.unique(isRedundantPair);               
    int NU = localpairs.size();
    printf(" remove %d of %d (%5.2f%%) localpairs\n",N0-NU,N0,100.*double(N0-NU)/N0);               
    N0 = crosspairs.size();               
    crosspairs.unique(isRedundantCross);               
    NU = crosspairs.size();               
    printf(" remove %d of %d (%5.2f%%) crosspairs\n",N0-NU,N0,100.*double(N0-NU)/N0);               
    //-------------------------------------------------------------
    // single reads redundant? 
    // remove only for variable read length platforms (454 / Helicos)
    //-------------------------------------------------------------     
    if ( (sLR/aLR)>0.25 ) {  // large variable read length
      N0 = dangle.size();               
      dangle.unique(isRedundantRead);               
      NU = dangle.size();               
      printf(" remove %d of %d (%5.2f%%) dangles\n",N0-NU,N0,100.*double(N0-NU)/N0);               
      N0 = singleton.size();               
      singleton.unique(isRedundantRead);               
      NU = singleton.size();               
      printf(" remove %d of %d (%5.2f%%) singletons\n",N0-NU,N0,100.*double(N0-NU)/N0);               
      N0 = umpairs.size();               
      umpairs.unique(isRedundantMulti);               
      NU = umpairs.size();               
      printf(" remove %d of %d (%5.2f%%) U-M\n",N0-NU,N0,100.*double(N0-NU)/N0);               
    }
    uniquified=2;
    
}

/*
void C_contig::uniquifyBam() {

    // check if already done, bug out if done    
    if (uniquified>2) return;

    int N0 = localpairs.size();               
    localpairs.unique(isRedundantPairBam);               
    int NU = localpairs.size();
    printf(" remove %d of %d (%5.2f%%) localpairs\n",N0-NU,N0,100.*double(N0-NU)/N0);               
    N0 = crosspairs.size();               
    crosspairs.unique(isRedundantCross);               
    NU = crosspairs.size();               
    printf(" remove %d of %d (%5.2f%%) crosspairs\n",N0-NU,N0,100.*double(N0-NU)/N0);               
    //-------------------------------------------------------------
    // single reads redundant? 
    // remove only for variable read length platforms (454 / Helicos)
    //-------------------------------------------------------------     
    if ( (sLR/aLR)>0.25 ) {  // large variable read length
      N0 = dangle.size();               
      dangle.unique(isRedundantRead);               
      NU = dangle.size();               
      printf(" remove %d of %d (%5.2f%%) dangles\n",N0-NU,N0,100.*double(N0-NU)/N0);               
      N0 = singleton.size();               
      singleton.unique(isRedundantRead);               
      NU = singleton.size();               
      printf(" remove %d of %d (%5.2f%%) singletons\n",N0-NU,N0,100.*double(N0-NU)/N0);               
      N0 = umpairs.size();               
      umpairs.unique(isRedundantMulti);               
      NU = umpairs.size();               
      printf(" remove %d of %d (%5.2f%%) U-M\n",N0-NU,N0,100.*double(N0-NU)/N0);               
    }
    uniquified=3;
}
*/

// I/O function for contigset
void C_contig::writePairs(string & outfilename) //const
{
  // format: optimized to loadFragments.m matlab script  
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out | ios::binary);
  if (!output) {
      cerr << "Unable to open file: " << outfilename << endl;
      return;
  }
  C_headerSpan h;
  h.setName = setName;
  h.contigName = contigName;
  int t =0;
  while(outfilename.find(h.spanext[t])==string::npos) t++;
  h.typeName = h.spanext[t];
  h.reclen =   3*sizeof(int)+2*sizeof(short)+6*sizeof(char);
  h.N = this->localpairs.size();
  h.write(output);
  list<C_localpair>::iterator i;
  for(i=localpairs.begin(); i != localpairs.end(); ++i) {
    unsigned int pos = (*i).pos;  
    int lm = (*i).lm;  
    char o = (*i).orient;  
    char q1 = (*i).q1;  
    char q2 = (*i).q2;  
    char mm1 = (*i).mm1;  
    char mm2 = (*i).mm2;  
    short len1 = (*i).len1;  
    short len2 = (*i).len2;  
    char constrain = (*i).constrain;  
    unsigned int ReadGroupCode0 = (*i).ReadGroupCode;  
    
    output.write(reinterpret_cast<const char *>(&pos), sizeof(int));
    output.write(reinterpret_cast<const char *>(&lm), sizeof(int));
    output.write(reinterpret_cast<const char *>(&o), sizeof(char));
    output.write(reinterpret_cast<const char *>(&len1), sizeof(short));
    output.write(reinterpret_cast<const char *>(&len2), sizeof(short));
    output.write(reinterpret_cast<const char *>(&q1), sizeof(char));    
    output.write(reinterpret_cast<const char *>(&q2), sizeof(char));    
    output.write(reinterpret_cast<const char *>(&mm1), sizeof(char));    
    output.write(reinterpret_cast<const char *>(&mm2), sizeof(char));    
    output.write(reinterpret_cast<const char *>(&constrain), sizeof(char));    
    output.write(reinterpret_cast<const char *>(&ReadGroupCode0), sizeof(unsigned int));    
  }
  output.close();
}

// I/O function for contigset
void C_contig::printPairs(string & outfilename) //const
{
  // format: optimized to loadFragments.m matlab script  
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out);
  if (!output) {
      cerr << "Unable to open file: " << outfilename << endl;
      return;
  }
  int V = 207;
  output << "Local Pair list, version " << V << endl;
  output << "Set " << setName << "\tContig " << contigName << endl;
  unsigned int nf = this->localpairs.size();
  output << "Number of Pairs " << nf << endl;
  output << " anchor  position lengthFragment orientation   q1 q2  mm1 mm2 ReadGroup "<< endl;
  //
  list<C_localpair>::iterator i;
  for(i=localpairs.begin(); i != localpairs.end(); ++i) {
    output << " " << (*i) << endl;
  }
  output.close();
}


// I/O function for contigset
void C_contig::printEnd(string & outfilename,  list<C_singleEnd> &  reads) //const
{
  // format: optimized to loadFragments.m matlab script  
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out);
  if (!output) {
      cerr << "Unable to open End output file: " << outfilename << endl;
      return;
  }
  int V = 207;
  output << "Read-end list, version " << V << endl;
  output << "Set " << setName << "\tContig " << contigName << endl;
  unsigned int ne = reads.size();
  output << "Number of Ends " << ne << endl;
  output << " anchor  position length sense   quality mm ReadGroup"<< endl;
  //
  list<C_singleEnd>::iterator i;
  for(i=reads.begin(); i != reads.end(); ++i) {
    output << " " << (*i) << endl;
  }
  output.close();
}


// I/O function for crosspairs
void C_contig::printCross(string & outfilename) //const
{
  // format: optimized to loadFragments.m matlab script  
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out);
  if (!output) {
      cerr << "Unable to open Cross output file: " << outfilename << endl;
      return;
  }
  int V = 207;
  output << "Cross pair list, version " << V << endl;
  output << "Set " << setName << "\tContig " << contigName << endl;
  unsigned int nc = crosspairs.size();
  output << "Number of Cross pairs " << nc << endl;
  output << " anchor0  position0 length0 sense0   quality0 mm0";
  output << " anchor1  position1 length1 sense1   quality1 mm1";
  output << " readGroup " << endl;
  //
  list<C_crosspair>::iterator i;
  for(i=crosspairs.begin(); i != crosspairs.end(); ++i) {
    output << " " << (*i).read[0] << "\t " << (*i).read[1]  << endl;
  }
  output.close();
}
// I/O function for crosspairs
void C_contig::printMulti(string & outfilename, list<C_umpair> &  um) //const
{
  // format: optimized to loadFragments.m matlab script  
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out);
  if (!output) {
      cerr << "Unable to open Retro output file: " << outfilename << endl;
      return;
  }
  int V = 207;
  output << "UM pair list, version " << V << endl;
  output << "Set " << setName << "\tContig " << contigName << endl;
  unsigned int nm = um.size();
  output << "Number of U-M pairs " << nm << endl;
  output << " anchor0  position0 length0 sense0   quality0 mm0 ";
  output << " anchor1  position1 length1 sense1   quality1 mm1 ";
  output << " Nmap  ReadGroup"<< endl;
  //
  list<C_umpair>::iterator i;
  for(i=um.begin(); i != um.end(); ++i) {
    output << " " << (*i).read[0] << "\t " << (*i).read[1]; 
    output << "\t " << (*i).nmap <<  "\t " << int2binary((*i).elements) << endl;
  }
  output.close();
}

// I/O function for contigset
void C_contig::printStats(string & outfilename) //const
{
  // format: optimized to loadFragments.m matlab script  
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out);
  if (!output) {
      cerr << "Unable to open Stats output file: " << outfilename << endl;
      return;
  }
  int V = 207;
  output << "stats, version " << V << endl;
  output << "Set " << setName << "\tContig " << contigName << endl;
  //
  output << "repeats " << endl;
  output <<   repeat.Stats << endl;
  output <<   repeat.Stats.h << endl;
  output << "pairs " << endl;
  output <<   pairStats << endl;
  output <<   pairStats.h << endl;
  output << "cross " << endl;
  output <<   crossStats << endl;
  output <<   crossStats.h << endl;
  output.close();
}


// I/O function for contigset
void C_contig::writeCross(string & outfilename) //const
{
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out | ios::binary);
  if (!output) {
      cerr << "Unable to open Cross output file: " << outfilename << endl;
      return;
  }
  C_headerSpan h;
  h.typeName = "cross";
  h.setName = setName;
  h.contigName = contigName;
  int t =0;
  while(outfilename.find(h.spanext[t])==string::npos) t++;
  h.typeName = h.spanext[t];
  h.reclen =   2*(sizeof(int)+2*sizeof(short)+3*sizeof(char))+sizeof(int);
  h.N = this->crosspairs.size();
  h.write(output);  // write version 
  //
  list<C_crosspair>::iterator i;
  for(i=crosspairs.begin(); i != crosspairs.end(); ++i) {
    for (int e=0; e<2; e++) {
      unsigned int p0 = (*i).read[e].pos;
      unsigned short l0 = (*i).read[e].len;          // length of this read aligment 
      unsigned short a0 = (*i).read[e].anchor;    // anchor index  
      char  s0 = (*i).read[e].sense;    // anchor index  
      char  q0 = (*i).read[e].q;    // anchor index  
      char  mm0 = (*i).read[e].mm;    // anchor index  
      //position
      output.write(reinterpret_cast<const char *>(&p0), sizeof(int));
      //length
      output.write(reinterpret_cast<const char *>(&l0), sizeof(short));
      output.write(reinterpret_cast<const char *>(&a0), sizeof(short));
      output.write(reinterpret_cast<const char *>(&s0), sizeof(char));
      output.write(reinterpret_cast<const char *>(&q0), sizeof(char));
      output.write(reinterpret_cast<const char *>(&mm0), sizeof(char));
    }
    unsigned int ReadGroupCode0 = (*i).ReadGroupCode;
    output.write(reinterpret_cast<const char *>(&ReadGroupCode0), sizeof(int));
    
  }
  output.close();
}

void C_contig::writeEnd(string & outfilename, list<C_singleEnd> &  reads) 
  {
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out | ios::binary);
  if (!output) {
      cerr << "Unable to open End output file: " << outfilename << endl;
      return;
  }
  C_headerSpan h;
  h.typeName = "End";
  h.setName = setName;
  h.contigName = contigName;
  int t =0;
  while(outfilename.find(h.spanext[t])==string::npos) t++;
  h.typeName = h.spanext[t];
  //
  h.reclen = 2*sizeof(int)+2*sizeof(short)+3*sizeof(char);
  h.N = reads.size();
  h.write(output);

  list<C_singleEnd>::iterator i;
  for(i=reads.begin(); i != reads.end(); ++i) {
    unsigned int p0 = (*i).pos;
    unsigned short l0 = (*i).len;                 // length of this read aligment 
    unsigned short a0 = (*i).anchor;              // anchor index  
    char  s0 = (*i).sense;                        // sense  
    char  q0 = (*i).q;                            // quality index  
    char  mm0 = (*i).mm;                          // mismatch  
    unsigned int ReadGroupCode0 = (*i).ReadGroupCode;    // library  index  
    //position
    output.write(reinterpret_cast<const char *>(&p0), sizeof(int));
    //length
    output.write(reinterpret_cast<const char *>(&l0), sizeof(short));
    output.write(reinterpret_cast<const char *>(&a0), sizeof(short));
    output.write(reinterpret_cast<const char *>(&s0), sizeof(char));
    output.write(reinterpret_cast<const char *>(&q0), sizeof(char));
    output.write(reinterpret_cast<const char *>(&mm0), sizeof(char));
    output.write(reinterpret_cast<const char *>(&ReadGroupCode0), sizeof(int));
  }
  output.close();
}

void C_contig::writeMulti(string & outfilename, list<C_umpair> &  um) 
  {
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out | ios::binary);
  if (!output) {
      cerr << "Unable to open Retro output file: " << outfilename << endl;
      return;
  }
  C_headerSpan h;
  h.typeName = "Multi";
  h.setName = setName;
  h.contigName = contigName;
  int Nt=h.spanext.size(),t;
  for (t=0; t<Nt; t++) {
    if (outfilename.find(h.spanext[t])!=string::npos) { break;}
  }
  if (t<Nt) { 
    h.typeName = h.spanext[t];
  } else {
    h.typeName = "multi.span";
  }  
  // two reads + nmap
  h.reclen =3*sizeof(int)+2*(sizeof(int)+2*sizeof(short)+3*sizeof(char));
  h.N = um.size();
  h.write(output);
  // loop
  list<C_umpair>::iterator i;
  for(i=um.begin(); i != um.end(); ++i) {
    for (int e=0; e<2; e++) {
      unsigned int p0 = (*i).read[e].pos;
      unsigned short l0 = (*i).read[e].len;       // length of this read aligment 
      unsigned short a0 = (*i).read[e].anchor;    // anchor index  
      char  s0 = (*i).read[e].sense;              // sense 
      char  q0 = (*i).read[e].q;                  // quality
      char mm0 = (*i).read[e].mm;                 // mm  
      //position
      output.write(reinterpret_cast<const char *>(&p0), sizeof(int));
      //length
      output.write(reinterpret_cast<const char *>(&l0), sizeof(short));
      output.write(reinterpret_cast<const char *>(&a0), sizeof(short));
      output.write(reinterpret_cast<const char *>(&s0), sizeof(char));
      output.write(reinterpret_cast<const char *>(&q0), sizeof(char));
      output.write(reinterpret_cast<const char *>(&mm0), sizeof(char));
    }
    int nmap= (*i).nmap;                          //nmapped positions
    output.write(reinterpret_cast<const char *>(&nmap), sizeof(int));
    int elements= (*i).elements;                  //nmapped positions
    output.write(reinterpret_cast<const char *>(&elements), sizeof(int));
    unsigned int ReadGroupCode= (*i).ReadGroupCode;   //library index 
    output.write(reinterpret_cast<const char *>(&ReadGroupCode), sizeof(int));
  }  
  output.close();
}

// I/O function for depth of coverage
void C_contig::writeMarker(string & outfilename, C_marker & m1) 
{
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out | ios::binary);
  if (!output) {
      cerr << "Unable to open Marker output file: " << outfilename << endl;
      return;
  }
  C_headerSpan h;
  h.setName = setName;
  h.contigName = contigName;
  // find typeName in spanext (build file extensions)
  for (int t=0; t< int(h.spanext.size()); t++) {
    if (!outfilename.find(h.spanext[t])==string::npos) break;
    h.typeName = h.spanext[t];
  }
  // find typeName in spanextc (cluster file extensions)
  for (int t=0; t<int(h.spanextc.size()); t++) {
    if (!outfilename.find(h.spanextc[t])==string::npos) break;
    h.typeName = h.spanextc[t];
  }
  if (h.typeName.size()==0) { 
    cerr << " unknown marker file extension " << outfilename << endl; 
  }
  h.reclen = sizeof(bool);
  h.N = m1.x.size();
  h.write(output);
  // old 204
  //copy(m1.x.begin(), m1.x.end(), ostreambuf_iterator<char>(output));
  int dp = DPCOMPRESS;
  char ibuf[dp];
  char cbuf[dp];
  int L = h.N;
  for (int p0=0; p0<L; p0+=dp)	{
    int d1 = dp;
    if ((p0+d1)>L) d1 = L-p0;
    for (int p=0; p<d1; p++) {
      ibuf[p]=m1.x[p0+p];
    }
    int nc = fastlz_compress(ibuf,d1, cbuf);	
    output.write(reinterpret_cast<const char *>(&p0), sizeof(int));
    output.write(reinterpret_cast<const char *>(&d1), sizeof(int));
    output.write(reinterpret_cast<const char *>(&nc), sizeof(int));
    output.write(reinterpret_cast<const char *>(cbuf), nc);
  }
  output.close(); 
}

// I/O function for depth of coverage
void C_contig::writeDepth(string & outfilename, C_depth & d1) 
{
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out | ios::binary);
  if (!output) {
      cerr << "Unable to open Depth output file: " << outfilename << endl;
      return;
  }
  // cout << " depth " <<endl;
  // cout << d1 << endl << endl;
  C_headerSpan h;
  h.setName = setName;
  h.contigName = contigName;
  h.light = d1.light;
  // find typeName in spanext (build file extensions)
  for (int t=0; t<int(h.spanext.size()); t++) {
    if (outfilename.find(h.spanext[t])!=string::npos) {
      h.typeName = h.spanext[t];
    }
  }
  // find typeName in spanextc (cluster file extensions)
  for (int t=0; t<int(h.spanextc.size()); t++) {
    if (outfilename.find(h.spanextc[t])==string::npos) {
      h.typeName = h.spanextc[t];
    }
  }
  if (h.typeName.size()==0) { 
    cerr << " unknown write depth file extension " << outfilename << endl; 
    size_t k = outfilename.find(contigName);
    if (k!=string::npos) {
       h.typeName = outfilename.substr(k+1);
    } 
    cerr << " unknown write depth file extension " << outfilename << endl; 
    cerr << " use " << h.typeName << endl; 

  }
  //
  h.reclen = sizeof(float);
  h.N = d1.n.size();
  h.write(output);
  int L=d1.n.size();
  //old version < 206
  /*
  for (int i=0; i<L; i++)	{
    output.write(reinterpret_cast<const char *>(&d1.n[i]), sizeof(float));
  }
  */
  int dp = DPCOMPRESS;
  float ibuf[dp];
  char  cbuf[dp*sizeof(float)];
  for (int p0=0; p0<L; p0+=dp)	{
    int dp1 = dp;
    if ((p0+dp1)>L) dp1 = L-p0;
    for (int p=0; p<dp1; p++) {
      ibuf[p]=d1.n[p0+p];
    }
    int nc = fastlz_compress(ibuf,dp1*sizeof(float), cbuf);	
    output.write(reinterpret_cast<const char *>(&p0), sizeof(int));
    output.write(reinterpret_cast<const char *>(&dp1), sizeof(int));
    output.write(reinterpret_cast<const char *>(&nc), sizeof(int));
    output.write(reinterpret_cast<const char *>(cbuf), nc);
  }
  output.close(); 
}

// I/O function for depth of coverage
void C_contig::writeDepth(string & outfilename, C_depth & d1, int binsize, int totbin) 
{
  // open output binary file. bomb if unable to open
  fstream output(outfilename.c_str(), ios::out | ios::binary);
  if (!output) {
      cerr << "Unable to open Depth output file: " << outfilename << endl;
      return;
  }
  cout << " depth " <<endl;
  cout << d1 << endl << endl;
  C_headerSpan h;
  h.setName = setName;
  h.contigName = contigName;
  light_t info;
  info.i64=0;  
  info.i16[0]=binsize;
  info.i16[1]=totbin;  
  h.light = info.i64;

  // find typeName in spanext (build file extensions)
  for (int t=0; t<int(h.spanext.size()); t++) {
    if (!outfilename.find(h.spanext[t])==string::npos) break;
    h.typeName = h.spanext[t];
  }
  // find typeName in spanextc (cluster file extensions)
  for (int t=0; t<int(h.spanextc.size()); t++) {
    if (!outfilename.find(h.spanextc[t])==string::npos) break;
    h.typeName = h.spanextc[t];
  }
  if (h.typeName.size()==0) { 
    cerr << " unknown write depth file extension " << outfilename << endl; 
  }
  //
  h.reclen = sizeof(float);
  h.N = d1.n.size();
  h.write(output);
  int L=d1.n.size();
  //old version < 206
  /*
  for (int i=0; i<L; i++)	{
    output.write(reinterpret_cast<const char *>(&d1.n[i]), sizeof(float));
  }
  */
  int dp = DPCOMPRESS;
  float ibuf[dp];
  char  cbuf[dp*sizeof(float)];
  for (int p0=0; p0<L; p0+=dp)	{
    int dp1 = dp;
    if ((p0+dp1)>L) dp1 = L-p0;
    for (int p=0; p<dp1; p++) {
      ibuf[p]=d1.n[p0+p];
    }
    int nc = fastlz_compress(ibuf,dp1*sizeof(float), cbuf);	
    output.write(reinterpret_cast<const char *>(&p0), sizeof(int));
    output.write(reinterpret_cast<const char *>(&dp1), sizeof(int));
    output.write(reinterpret_cast<const char *>(&nc), sizeof(int));
    output.write(reinterpret_cast<const char *>(cbuf), nc);
  }
  output.close(); 
}


// I/O function for depth of coverage
C_depth  C_contig::loadDepth(string & infilename) 
{
   C_depth x1;
  // open output binary file. bomb if unable to open
  fstream input(infilename.c_str(), ios::in  | ios::binary);
  if (!input) {
      cerr << "Unable to open Depth input file: " << infilename << endl;
      return x1;
  }
  C_headerSpan  h(input); 
  int L=h.N;
  C_depth d1(L);
  d1.light=h.light;
  float  n1=0;
  d1.Stats.N=0;
  if (h.V<206) {  //non-compressed array 
    for (int i=0; i<L; i++)	{
      input.read(reinterpret_cast < char * > (&n1), sizeof(n1));
      d1.n[i]=n1;
      d1.Stats.N+=int(n1);
    }
  } else { // compressed 
    float ibuf[DPCOMPRESS];
    char  cbuf[DPCOMPRESS*sizeof(float)];
    int p1=0; 
    int p0, dp, nc;
    while (p1<L) {
      input.read(reinterpret_cast < char * > (&p0), sizeof(p0));
      input.read(reinterpret_cast < char * > (&dp), sizeof(dp));
      input.read(reinterpret_cast < char * > (&nc), sizeof(nc));
      input.read(reinterpret_cast < char * > (cbuf), nc);
      int nd = fastlz_decompress(cbuf,nc, ibuf,int(dp*1.1*sizeof(float)));	
      p1 = p0+dp;
      if (nd!=(dp*int(sizeof(float)))) {
        cerr << "Unable to decompress Depth input file: " << infilename << endl;
        return x1;
      }
      for (int p=0; p<dp; p++) {
        d1.n[p0+p]=ibuf[p];
        d1.Stats.N+=int(d1.n[p0+p]);
      }
    }
  }
  input.close(); 
  return d1;
}

C_marker   C_contig::loadMarker(string & infilename) 
{
  C_marker x1;
  // open output binary file. bomb if unable to open
  fstream input(infilename.c_str(), ios::in  | ios::binary);
  if (!input) {
      cerr << "Unable to open Marker input file: " << infilename << endl;
      return x1;
  }
  C_headerSpan h(input);
  int L = h.N;
  C_marker m1(L);
  bool b1=false;
  if (h.V<206) {  //non-compressed array 
    for (int i=0; i<L; i++)	{
      input.read(reinterpret_cast < char * > (&b1), sizeof(b1));
      m1.x[i]=b1;
    }
  } else { // compressed buffer

    char ibuf[DPCOMPRESS];
    char  cbuf[DPCOMPRESS];
    int p1=0; 
    int p0, dp, nc;
    while (p1<L) {
      input.read(reinterpret_cast < char * > (&p0), sizeof(p0));
      input.read(reinterpret_cast < char * > (&dp), sizeof(dp));
      input.read(reinterpret_cast < char * > (&nc), sizeof(nc));
      input.read(reinterpret_cast < char * > (cbuf), nc);
      int nd = fastlz_decompress(cbuf,nc, ibuf,dp);	
      p1 = p0+dp;
      if (nd!=dp) {
        cerr << "Unable to decompress Marker input file: " << infilename << endl;
        return x1;
      }
      for (int p=0; p<dp; p++) {
        m1.x[p0+p]=ibuf[p];
      }
    }
  }
  input.close(); 
  return m1;
}

void  C_contig::loadRepeat(string & infilename) 
{
  C_depth r1=loadDepth(infilename);
  if (repeat.n.size()==0) {
     repeat = r1;
  } else {
     if (repeat.n.size()==r1.n.size()) {
        for (int i=0; i<int(r1.n.size()); i++) {
          repeat.n[i]+=r1.n[i];
        }
        // stats??
     } else { 
       cerr << " mixed length repeat.span " << endl;
       cerr << infilename << "\t length " << r1.n.size() << endl;
       cerr << "existing repeat length " << repeat.n.size() << endl;
       exit(-1);
    }
  }
}

void  C_contig::loadPairs(string & infilename) 
{
  fstream input(infilename.c_str(), ios::in  | ios::binary);
  if (!input) {
      cerr << "Unable to open local pairs input file: " << infilename << endl;
  }
  C_headerSpan h(input);

  int N= h.N;
  contigName=h.contigName;
  setName=h.setName;
  // read data
  for(int i=0; i<N; i++)  {
    C_localpair p1;
    input.read(reinterpret_cast < char * > (&p1.pos), sizeof(int));
    input.read(reinterpret_cast < char * > (&p1.lm), sizeof(int));
    input.read(reinterpret_cast < char * > (&p1.orient), sizeof(char));
    input.read(reinterpret_cast < char * > (&p1.len1), sizeof(short));
    input.read(reinterpret_cast < char * > (&p1.len2), sizeof(short));
    input.read(reinterpret_cast < char * > (&p1.q1), sizeof(char));
    if (h.V>206) {
        input.read(reinterpret_cast < char * > (&p1.q2), sizeof(char));
        input.read(reinterpret_cast < char * > (&p1.mm1), sizeof(char));
        input.read(reinterpret_cast < char * > (&p1.mm2), sizeof(char));
    } else {
      p1.q2=0;
      p1.mm1=0;
      p1.mm2=0;
    }
    if (h.V>203) {
        input.read(reinterpret_cast < char * > (&p1.constrain), sizeof(char));
    } else {
        p1.constrain=0;
    }
    if (h.V>206) {
        input.read(reinterpret_cast < char * > (&p1.ReadGroupCode), sizeof(int));
    } else {
        p1.ReadGroupCode=0;
    }
    this->localpairs.push_back(p1);
  }
  input.close();
}

void  C_contig::loadCross(string & infilename) 
{
  fstream input(infilename.c_str(), ios::in  | ios::binary);
  if (!input) {
      cerr << "Unable to open Cross input file: " << infilename << endl;
  }
  C_headerSpan h(input);

  int N = h.N;
  //
  for(int i=0; i<N; i++)  {
    C_crosspair c1;
    //read 0
    input.read(reinterpret_cast < char * > (&c1.read[0].pos), sizeof(int));
    input.read(reinterpret_cast < char * > (&c1.read[0].len), sizeof(short));
    input.read(reinterpret_cast < char * > (&c1.read[0].anchor), sizeof(short));
    input.read(reinterpret_cast < char * > (&c1.read[0].sense), sizeof(char));
    input.read(reinterpret_cast < char * > (&c1.read[0].q), sizeof(char));
    if (h.V>206) {
        input.read(reinterpret_cast < char * > (&c1.read[0].mm), sizeof(char));
    } else {
       c1.read[0].mm=0;
    }
    //read 1
    input.read(reinterpret_cast < char * > (&c1.read[1].pos), sizeof(int));
    input.read(reinterpret_cast < char * > (&c1.read[1].len), sizeof(short));
    input.read(reinterpret_cast < char * > (&c1.read[1].anchor), sizeof(short));
    input.read(reinterpret_cast < char * > (&c1.read[1].sense), sizeof(char));
    input.read(reinterpret_cast < char * > (&c1.read[1].q), sizeof(char));
    if (h.V>206) {
        input.read(reinterpret_cast < char * > (&c1.read[1].mm), sizeof(char));
        input.read(reinterpret_cast < char * > (&c1.ReadGroupCode), sizeof(int));
    } else {
       c1.read[1].mm=0;
       c1.ReadGroupCode=0;
    }
    this->crosspairs.push_back(c1);
  }
  input.close();
}                  

void  C_contig::loadMultipairs(string & infilename) 
{
  list<C_umpair> x1 =loadMulti(infilename);
  umpairs.merge(x1);
}

void  C_contig::loadDangle(string & infilename) 
{
  list<C_singleEnd> x1 =loadEnd(infilename);
  dangle.merge(x1);
}


void  C_contig::loadSingleton(string & infilename) 
{
  list<C_singleEnd> x1 =loadEnd(infilename);
  singleton.merge(x1);
}

bool C_contig::isRedundantPair(const C_localpair &p1, const C_localpair &p2) 
{
   if( p1.pos != p2.pos) return false;
   if( p1.anchor != p2.anchor) return false;
   if( p1.lm != p2.lm) return false;
   if( p1.orient != p2.orient) return false;
   return true;
}


/*
//-----------------------------------------------------------------------------
// special test for Bam fragments which are expected to be redundant in pairs
// except for items like mm1 mm2 q1 q2
// this function overwrites reciprocal information from each end such that 
// p1 and p2 are really redundant 
//-----------------------------------------------------------------------------
bool C_contig::isRedundantPairBam(C_localpair &p1, C_localpair &p2) 
{
   if( p1.pos != p2.pos) return false;
   if( p1.anchor != p2.anchor) return false;
   if( p1.lm != p2.lm) return false;
   if( p1.orient != p2.orient) return false;
   // impose redundant information
   // determine which end is the "mate" with incomplete info
   bool fixp1m1=((p1.mm1==99)&&(p1.mm2!=99));
   bool fixp2m1=((p2.mm1==99)&&(p2.mm2!=99));
   bool fixp1m2=((p1.mm2==99)&&(p1.mm1!=99));
   bool fixp2m2=((p2.mm2==99)&&(p2.mm1!=99));
   // p1m2 and p2m1 are ok
   if (fixp1m1) {
      p1.q1=p2.q1;
      p1.len1=p2.len1;
      p1.mm1=p2.mm1;
      return true;      
   } 
   if (fixp2m2) {
      p2.q2=p1.q2;
      p2.len2=p1.len2;
      p2.mm2=p1.mm2;
      return true;      
   }
    // p1m1 and p2m2 are ok
   if (fixp1m2) {
      p1.q2=p2.q2;
      p1.len2=p2.len2;
      p1.mm2=p2.mm2;
      return true;      
   }
   if (fixp2m1) {
      p2.q1=p1.q1;
      p2.len1=p1.len1;
      p2.mm1=p1.mm1;
      return true;      
   } 
   // other combinations problematic
   cerr << "problem Bam redundancy" << endl;
   cerr << p1 << endl;
   cerr << p2 << endl;
   return true;
}

 */


bool C_contig::isRedundantRead(const C_singleEnd &p1, const C_singleEnd &p2)
{
   if( !(p1.pos == p2.pos)) return false;
   if( !(p1.anchor == p2.anchor)) return false;
   if( !(p1.len == p2.len)) return false;
   return true;
}

bool C_contig::isRedundantMulti(const C_umpair &p1, const C_umpair &p2)
{
   if( !(p1.read[0].pos == p2.read[0].pos)) return false;
   if( !(p1.read[0].anchor == p2.read[0].anchor)) return false;
   if( !(p1.read[0].len == p2.read[0].len)) return false;
   // check number of mappings on other end... if within 10% then redundant
   double xnm=2.0*fabs(double(p1.nmap-p2.nmap))/double(p1.nmap+p2.nmap);   
   if( xnm>0.1) return false;
   return true;
}

bool C_contig::isRedundantCross(const C_crosspair &p1, const C_crosspair &p2)
{
   if( !(p1.read[0].pos == p2.read[0].pos)) return false;
   if( !(p1.read[1].pos == p2.read[1].pos)) return false;
   if( !(p1.read[0].anchor == p2.read[0].anchor)) return false;
   if( !(p1.read[1].anchor == p2.read[1].anchor)) return false;
   if( !(p1.read[0].sense == p2.read[0].sense)) return false;
   if( !(p1.read[1].sense == p2.read[1].sense)) return false;
   return true;
}
    
 
list<C_singleEnd>   C_contig::loadEnd(string &   infilename) 
{
  list<C_singleEnd> x1;
  fstream input(infilename.c_str(), ios::in|ios::binary);
  if (!input) {
      cerr << "Unable to open read end input file: " << infilename << endl;
      return x1;
  }
  C_headerSpan h(input);
  int N = h.N;
  for(int i=0; i<N; i++)  {
    C_singleEnd r1;
    //position
    input.read(reinterpret_cast < char * > (&r1.pos), sizeof(int));
    //length
    input.read(reinterpret_cast < char * > (&r1.len), sizeof(short));
    input.read(reinterpret_cast < char * > (&r1.anchor), sizeof(short));
    input.read(reinterpret_cast < char * > (&r1.sense), sizeof(char));
    input.read(reinterpret_cast < char * > (&r1.q), sizeof(char));
    if (h.V>206) {
      input.read(reinterpret_cast < char * > (&r1.mm), sizeof(char));
      input.read(reinterpret_cast < char * > (&r1.ReadGroupCode), sizeof(int));
    } else {
       r1.mm=0;
       r1.ReadGroupCode=0;
    }
    x1.push_back(r1);
  }
  input.close();
  return x1;
}        

list<C_umpair>   C_contig::loadMulti(string &   infilename) 
{
  list<C_umpair> x1;
  fstream input(infilename.c_str(), ios::in|ios::binary);
  if (!input) {
      cerr << "Unable to open read retro input file: " << infilename << endl;
      return x1;
  }
  C_headerSpan h(input);
  int reclen0 =sizeof(int)+2*(sizeof(int)+2*sizeof(short)+2*sizeof(char));
  //int reclen1 =reclen0+sizeof(int);
  int N = h.N;
  for(int i=0; i<N; i++)  {
    C_umpair c1;
    //read 0 - unique end
    input.read(reinterpret_cast < char * > (&c1.read[0].pos), sizeof(int));
    input.read(reinterpret_cast < char * > (&c1.read[0].len), sizeof(short));
    input.read(reinterpret_cast < char * > (&c1.read[0].anchor), sizeof(short));
    input.read(reinterpret_cast < char * > (&c1.read[0].sense), sizeof(char));
    input.read(reinterpret_cast < char * > (&c1.read[0].q), sizeof(char));
    if (h.V>206) {
        input.read(reinterpret_cast < char * > (&c1.read[0].mm), sizeof(char));
    } else {
      c1.read[0].mm=0;
    }
    //read 1 - multiple map
    input.read(reinterpret_cast < char * > (&c1.read[1].pos), sizeof(int));
    input.read(reinterpret_cast < char * > (&c1.read[1].len), sizeof(short));
    input.read(reinterpret_cast < char * > (&c1.read[1].anchor), sizeof(short));
    input.read(reinterpret_cast < char * > (&c1.read[1].sense), sizeof(char));
    input.read(reinterpret_cast < char * > (&c1.read[1].q), sizeof(char));
    if (h.V>206) {
        input.read(reinterpret_cast < char * > (&c1.read[1].mm), sizeof(char));
    } else {
      c1.read[1].mm=0;
    }
    // nmap
    input.read(reinterpret_cast < char * > (&c1.nmap), sizeof(int));
    // elements
    if (int(h.reclen)>int(reclen0)) {
      input.read(reinterpret_cast < char * > (&c1.elements), sizeof(int));
    }
    if (h.V>206) {
        input.read(reinterpret_cast < char * > (&c1.ReadGroupCode), sizeof(int));
    } else {
      c1.ReadGroupCode=0;
    }
    x1.push_back(c1);
  }
  input.close();
  return x1;
}        
    
C_set::C_set(string & setName1, RunControlParameters & pars1) {
  setSetName(setName1);  
  this->fileName = "";
  this->pars=pars1;
  this->Nfrag = 0;  
  this->Npair = 0;  
  this->Npair_00 = 0;    
  this->Npair_01 = 0;  
  this->Npair_0N = 0;
  this->Npair_11 = 0;
  this->Npair_11o = 0;
  this->Npair_1N = 0;
  this->Npair_NN = 0;
	
  // Initialize stats 
  // labels: 0=no map; 1=one-map; N=multi-map;
  string pairCountLab1[] = {"0-0","0-1","0-N","X","1-1","1-N","X","X","N-N"};
  vector<string> pairCountLab(pairCountLab1, pairCountLab1 + 9);
  // labels: F1R2=(first in pos order)(next); F=forward sense, R=reverse;  1=mate 1, 2=mate 2;
	
  string pairModelLab1[] = {"F1F2","F1R2","R1F2","R1R2","F2F1","F2R1","R2F1","R2R1"};
  vector<string> pairModelLab(pairModelLab1, pairModelLab1 + 8);
  
  repeatStats.Initialize(1001,-0.5,1000.5);  
  repeatStats.h.setTitle("NA read mapping multiplicity");
  
  lengthStats.Initialize(1001,-0.5,1000.5);  
  lengthStats.h.setTitle("LR read length");

  fragStats.Initialize(10100,-99.5,10000.5);  
  fragStats.h.setTitle("LF fragment mapping length");

  spanStats.Initialize(10100,-99.5,10000.5);   
  spanStats.h.setTitle("SL span mapping length");
	
  qStats.Initialize(101,-0.5,100.5); 
  qStats.h.setTitle("RQ read map quality");
  
  pairCountStats.Initialize(9,-0.5,8.5); 
  pairCountStats.h.setTitle("RC pair multiplicity combinations");
  pairCountStats.h.setXlabel("combo");
  pairCountStats.h.setBinLabels(pairCountLab);

  pairModelStats.Initialize(8,0.5,8.5); 
  pairModelStats.h.setTitle("PM pair model combinations");
  pairModelStats.h.setXlabel("model");
  pairModelStats.h.setBinLabels(pairModelLab);
	
  refStats.Initialize(101,-0.5,100.5); 
  refStats.h.setTitle("RS reference ID");
	
	
  // Debug
  /*
  this->Npair_ALU = 0;
  this->Npair_ALUC = 0;
  */
  //
  elementMinAnchor=255; // init min element anchor index

}

void C_set::initContigs() {
  for (size_t i = 0; i<this->anchors.names.size(); i++) {
      string name= this->anchors.names[i];
      int L = 0;
      if (this->anchors.use[i]) { L = this->anchors.L[name]; }

      // turn off repeat checking by setting L<0
  		bool doRepeatCheck = false; //pars.getDoRepeatCheck();
      
      C_contig contig1(name,L, doRepeatCheck);
      contig1.setName=setName;
      this->contig[name]=contig1;
      // add anchor info to each contig object
      this->contig[name].anchors=this->anchors;
  }
}

void C_set::processRepeat(C_readmaps & read1, int MinMap=2) {
  //if (!pars.getDoRepeatCheck()) return;
  int N = read1.align.size();     
  if (N<MinMap) {return;}
  for (int i = 0 ; i<N; i++) {
    unsigned int p = read1.align[i].pos;
    unsigned short a = read1.align[i].anchor;
    string name = anchors.names[a];
    if (this->contig[name].Length>0) {
        contig[name].repeat.n[p]++;
    }
  }
}

void C_set::processSingleton(C_readmaps & read1,unsigned int ReadGroupCode) {
  int N = read1.align.size();     
  if (N!=1) {return;}
  unsigned short a = read1.align[0].anchor;
  string name = anchors.names[a];
  if (this->contig[name].Length>0) {
      C_singleEnd r1;
      r1.pos=read1.align[0].pos;
      r1.len=read1.align[0].len;
      r1.anchor=read1.align[0].anchor;
      r1.sense=read1.align[0].sense;
      r1.q=read1.align[0].q;
      r1.mm=read1.align[0].mm;
      r1.ReadGroupCode=ReadGroupCode;
      contig[name].singleton.push_back(r1);
  }
}


//------------------------------------------------------------------------------
// sort out fragment types to fill spanner records 
//------------------------------------------------------------------------------
void C_set::processPair(C_pairedread & pair1, char constrain) {

  //int N0 = pair1.read[0].align[0].nmap;     
  //int N1 = pair1.read[1].align[0].nmap;   
  int N0 = pair1.read[0].Nalign;     
  int N1 = pair1.read[1].Nalign;   
  
  //----------------------------------------------------------------------------
  // bump N if the unique hit is ANY element
  // prevent mobile-element hits from going into UU pairs
  // want them to go to UM pairs
  //----------------------------------------------------------------------------
  if ( (pair1.read[0].element!=0)&&(N0==1) ) {
     N0++;
  }
  if ( (pair1.read[1].element!=0)&&(N1==1) ) {
     N1++;
  }
  //Npair++;

  //----------------------------------------------------------------------------
  //  U0 pairs
  //----------------------------------------------------------------------------
  if ((N0+N1)==1) {
    //=========================
    // 0 1 / 1 0  dangling reads
    //=========================
    Npair_01++;
    int e = (N0==1 ? 0: 1);
		if (pair1.read[e].align.size()>0) {
			unsigned short a = pair1.read[e].align[0].anchor;
			string name = anchors.names[a];
			if (this->contig[name].Length>0) {
        //char s = pair1.read[e].align[0].sense;
        C_singleEnd r1;
        r1.pos=pair1.read[e].align[0].pos;
        r1.len=pair1.read[e].align[0].len;
        r1.anchor=pair1.read[e].align[0].anchor;
        r1.sense=pair1.read[e].align[0].sense;
        r1.q=pair1.read[e].align[0].q;
        r1.mm=pair1.read[e].align[0].mm;
        r1.ReadGroupCode=pair1.ReadGroupCode;
        contig[name].dangle.push_back(r1);
			}
		}
		

  //----------------------------------------------------------------------------
  //  UM pairs
  //---------------------------------------------------------------------------    
  } else if ( (N0+N1)>2 && (N0==1 || N1==1) ) {
    //--------------------------------------------------------------------------
    // unique read on one end - many on other end 1 N
    //--------------------------------------------------------------------------
    Npair_1N++;
    int e = (N0==1 ? 0: 1);
    unsigned short a = pair1.read[e].align[0].anchor;
    string name = anchors.names[a];
    if (this->contig[name].Length>0) {
    
        //----------------------------------------------------------------------
        // retro pairs - Starts "F" - Ends "R"
        //----------------------------------------------------------------------
        //char s = pair1.read[e].align[0].sense;
        //int em = (e==1 ? 0: 1);
        //int nmap = pair1.read[em].align.size();
        //--------------------------------------------------------------------------
        
        // Debug
        /*
        if ( ((pair1.read[0].element&1)>0) || ((pair1.read[1].element&1)>0)  ) {
            Npair_ALU++;
            int lm = Fraglength(pair1);
            if (abs(lm-140)<70) { 
              Npair_ALUC++;
            } else {
              // cerr << " missed aluy " << lm  <<  " Nalu " << Npair_ALU << " Naluc " << Npair_ALUC << endl;
            }
        }
        */
          
        C_umpair r1(pair1,anchors);
        contig[name].umpairs.push_back(r1);
  
        /*
        if (s=='F') {
            contig[name].retroStart.push_back(r1);
        } else {
            contig[name].retroEnd.push_back(r1);
        }
        */
    }    
    
  //----------------------------------------------------------------------------
  // UU pairs
  //----------------------------------------------------------------------------
  } else if ((N0*N1)==1) {
    //=========================
    // unique pairs 
    //=========================
    Npair_11++;
    unsigned short a0 = pair1.read[0].align[0].anchor;
    unsigned short a1 = pair1.read[1].align[0].anchor;
    if (a0==a1) {
      //===========
      // local pair
      //===========
      C_localpair localpair1(pair1,constrain);
      string name = anchors.names[a0];
      if (this->contig[name].Length>0) {
        contig[name].localpairs.push_back(localpair1);
      }
    } else {
      //==================
      // cross-anchor pair
      //==================
      for (int e=0; e<2; e++) {
        unsigned short a = pair1.read[e].align[0].anchor;
        string name = anchors.names[a];
        int e1 = (e==0? 1: 0);
        if (this->contig[name].Length>0) {
            C_crosspair cross1(pair1.read[e].align[0],pair1.read[e1].align[0],pair1.ReadGroupCode);
            contig[name].crosspairs.push_back(cross1);
        }
      }
    } 
  } else if (N0==0 && N1==0) {
     //-------------------------------------------------------------------------
     // No mappings on either end 
     //-------------------------------------------------------------------------
     Npair_00++;
  } else if ( (N0*N1==0)&& ((N0+N1)>1) ) {
     //-------------------------------------------------------------------------
     // many mappings on one end - none on other end
     //-------------------------------------------------------------------------
     Npair_0N++;
  } else {
     //-------------------------------------------------------------------------
     // many mappings on both ends
     //-------------------------------------------------------------------------
     Npair_NN++;
  }
}

int C_set::Fraglength(C_pairedread & pair1) {
	int lm = -100000;
	unsigned int p0 = pair1.read[0].align[0].pos;
	unsigned int p1 = pair1.read[1].align[0].pos;
	unsigned short l0 = pair1.read[0].align[0].len;
	unsigned short l1 = pair1.read[1].align[0].len;
	unsigned char c0 = pair1.read[0].align[0].anchor;
	unsigned char c1 = pair1.read[1].align[0].anchor;
	unsigned char o0 = pair1.read[0].align[0].sense;
	unsigned char o1 = pair1.read[1].align[0].sense;
	if (c0==c1) {
		if ((o0=='F')&&(o1=='R')) {
			lm = p1+l1-p0-1;
		} else if ((o0=='R')&&(o1=='F')) {
			lm = p0+l0-p1-1;
		}
	}
	return lm;
}

int C_set::spanLength(C_pairedread & pair1) {
	int sl = -100000;
	unsigned int p0 = pair1.read[0].align[0].pos;
	unsigned int p1 = pair1.read[1].align[0].pos;
	unsigned short l0 = pair1.read[0].align[0].len;
	unsigned short l1 = pair1.read[1].align[0].len;
	unsigned char c0 = pair1.read[0].align[0].anchor;
	unsigned char c1 = pair1.read[1].align[0].anchor;
	unsigned char o0 = pair1.read[0].align[0].sense;
	unsigned char o1 = pair1.read[1].align[0].sense;
	if (c0==c1) {
		if ((o0=='F')&&(o1=='R')) {
			sl = p1-(p0+l0);
		} else if ((o0=='R')&&(o1=='F')) {
			sl = p0-(p1+l1);
		}
	}
	return sl;
}

//----------------------------------------------------------
// Michael's models... 20100803_Models in MosiakSort.pptx
//----------------------------------------------------------
int C_set::pairModel(C_pairedread & pair1) {
	int m = -1;
	
	unsigned int p0 = pair1.read[0].align[0].pos;
	unsigned int p1 = pair1.read[1].align[0].pos;
	unsigned char c0 = pair1.read[0].align[0].anchor;
	unsigned char c1 = pair1.read[1].align[0].anchor;
	unsigned char o0 = pair1.read[0].align[0].sense;
	unsigned char o1 = pair1.read[1].align[0].sense;
	//  3'/5' ENDS, F/R Polarity, 1/2 read order
	// {"F1F2","F1R2","R1F2","R1R2","F2F1","F2R1","R2F1","R2R1"};
	if (c0==c1) {
		if ((o0=='F')&&(o1=='R')&&(p0<=p1)) {
			m=2; // Illumina PE  F1R2
		} else if ((o0=='F')&&(o1=='F')&&(p0<=p1))  {
			m=1; // SOLiD PE  F1F2
		} else if ((o0=='R')&&(o1=='F')&&(p0<=p1))  {
			m=3; // R1F2
		} else if ((o0=='R')&&(o1=='R')&&(p0<=p1))  {
			m=4;  // 454 Mate Pairs R1R2
	  } else if ((o0=='F')&&(o1=='F')&&(p0>p1))  {
			m=5;  // 454 Mate Pairs  F2F1
		} else if ((o0=='R')&&(o1=='F')&&(p0>p1))  {
			m=6;  // Illumina PE  R2F1
		}	else if ((o0=='F')&&(o1=='R')&&(p0<p1)) {
			m=7; // F2R1				
		} else if ((o0=='R')&&(o1=='R')&&(p0>p1))  {
			m=8; // SOLiD PE  R2R1
		}  
	}
	return m;
}


void C_set::calcStats() {
  C_contigs::const_iterator iterContig;   
  for (iterContig = contig.begin();	iterContig != contig.end(); iterContig++) {
      string name = iterContig->first;
      if (contig[name].Length>0) { 
        cout << "calc stats for contig " << name << endl;
        // fill depth & repeat vector before uniquifying pairs
        contig[name].calcStats();
      }
  }
}

void C_set::uniquify() {
  if (pars.getDupRemove()<1) return;
  C_contigs::const_iterator iterContig;   
  for (iterContig = contig.begin();	iterContig != contig.end(); iterContig++) {
      string name = iterContig->first;
      if (contig[name].Length>0) { 
        cout << "uniquify contig " << name << endl;
        // fill depth & repeat vector before uniquifying pairs
        contig[name].uniquify();
      }
  }
}

void C_set::sort() {
  C_contigs::const_iterator iterContig;   
  for (iterContig = contig.begin();	iterContig != contig.end(); iterContig++) {
      string name = iterContig->first;
      if (contig[name].Length>0) { 
        cout << "sort contig " << name << endl;
        // fill depth & repeat vector before uniquifying pairs
        contig[name].sort();
      }
  }
}
void C_set::countReads(int binsize) {
  C_contigs::const_iterator iterContig;   
  for (iterContig = contig.begin();	iterContig != contig.end(); iterContig++) {
      string name = iterContig->first;
      if (contig[name].Length>0) { 
        cout << "calc Depth contig " << name << endl;
        contig[name].read_counts = contig[name].countReads(binsize);
      }
  }
}
 
void C_set::calcDepth() {
  C_contigs::const_iterator iterContig;   
  int Qmin = pars.getQmin();
  for (iterContig = contig.begin();	iterContig != contig.end(); iterContig++) {
      string name = iterContig->first;
      if (contig[name].Length>0) { 
        cout << "calc Depth contig " << name << endl;
        contig[name].calcDepth(Qmin);
      }
  }
}


void C_set::write() {
  C_contigs::const_iterator iterContig;   
  string area = pars.getOutputDir();
  int k = 1+area.find_last_of("/");
  string subdir = area.substr(k);
  if (area.size()>0) { area = area+"/";}
  string prefix = pars.getPrefix();
  if (prefix.size()>0) { prefix = prefix+".";}
  
        
  string sn1 = setName;
  // dont need redundant setName if subdirectory is already the setName
  // if (sn1==subdir) sn1="";
  string fn1 = fileName;
  size_t f1 = fn1.find(sn1); 
  if ((fn1.size()>0)&&(f1==string::npos)) sn1=fn1+"."+sn1;      
  if (sn1.size()>0) sn1 = sn1+".";

  string anchorfile = area+prefix+sn1+"anchors.txt";
  anchors.printAnchorInfo(anchorfile);          
  
  cout << "write span output  " <<  endl;
    
    
  string fname = area+prefix+sn1+"library.span";
  cout << "\t " << fname << "\t " << libraries.libmap.size() << endl;
  libraries.writeLibraryInfo(fname,setName);
  
  for (iterContig = contig.begin();	iterContig != contig.end(); iterContig++) {
      string cn1 = iterContig->first;      
      string basename = area+prefix+sn1+cn1;
      if (setName==cn1) basename = area+prefix+cn1;

      // replace evil character "|" with benign "_" 
      size_t found=basename.find("|");
      while (found!=string::npos) {
        basename.replace(found,1,"_");
        found=basename.find("|");
      }
  
      if (contig[cn1].Length>0) { 
        //cout << "write output for contig " << cn1 << endl;
        fname = basename+".pair.span";
        cout << "\t " << fname << "\t " << contig[cn1].localpairs.size() << endl;
        contig[cn1].writePairs(fname);
        fname = basename+".cross.span";
        cout << "\t " << fname << "\t " << contig[cn1].crosspairs.size() << endl;
        contig[cn1].writeCross(fname);
        fname = basename+".repeat.span";
        cout << "\t " << fname << "\t " << contig[cn1].repeat.n.size() << endl;
        contig[cn1].writeDepth(fname, contig[cn1].repeat);
        fname = basename+".dangle.span";
        cout << "\t " << fname << "\t " << contig[cn1].dangle.size() << endl;
        contig[cn1].writeEnd(fname, contig[cn1].dangle);
        fname = basename+".multi.span";
        cout << "\t " << fname << "\t " << contig[cn1].umpairs.size() << endl;
        contig[cn1].writeMulti(fname, contig[cn1].umpairs);
      }
   }
}

void C_set::printOut() {
  C_contigs::const_iterator iterContig;   
  string area = pars.getOutputDir();
  int k = 1+area.find_last_of("/");
  string subdir = area.substr(k);
  if (area.size()>0) { area = area+"/";}
  string prefix = pars.getPrefix();
  if (prefix.size()>0) { prefix = prefix+".";}
  
  string sn1 = setName;
  // dont need redundant setName if subdirectory is already the setName
  // if (sn1==subdir) sn1="";
  string fn1 = fileName;
  size_t f1 = fn1.find(sn1); 
  if ((fn1.size()>0)&&(f1==string::npos)) sn1=fn1+"."+sn1;    
  if (sn1.size()>0) { sn1 = sn1+".";}
  
  for (iterContig = contig.begin();	iterContig != contig.end(); iterContig++) {
      string cn1 = iterContig->first;      
      string basename = area+prefix+sn1+cn1;  
      if (setName==cn1) basename = area+prefix+cn1;

      
      
      // replace evil character "|" with benign "_" 
      size_t found=basename.find("|");
      while (found!=string::npos) {
        basename.replace(found,1,"_");
        found=basename.find("|");
      }
      
      if (contig[cn1].Length>0) { 
        cout << "print output for contig " << cn1 << endl;
        string fname = basename+".pair.span.txt";
        cout << "\t " << fname << "\t " << contig[cn1].localpairs.size() << endl;
        contig[cn1].printPairs(fname);
        fname = basename+".cross.span.txt";
        cout << "\t " << fname << "\t " << contig[cn1].crosspairs.size() << endl;
        contig[cn1].printCross(fname);
        /*
        fname = basename+".repeat.span";
        cout << "\t " << fname << endl;
        contig[cn1].writeMarker(fname, contig[cn1].repeat);
        */
        fname = basename+".dangle.span.txt";
        cout << "\t " << fname << "\t " << contig[cn1].dangle.size() << endl;
        contig[cn1].printEnd(fname, contig[cn1].dangle);
		fname = basename+".stat.span.txt";
        cout << "\t " << fname << endl;
        contig[cn1].printStats(fname);
        //
        fname = basename+".multi.span.txt";
        cout << "\t " << fname << "\t " << contig[cn1].umpairs.size() << endl;
        contig[cn1].printMulti(fname, contig[cn1].umpairs);
      }
   }
}

void C_set::writeDepth() {
  C_contigs::const_iterator iterContig;   
  string area = pars.getOutputDir();
  int k = 1+area.find_last_of("/");
  string subdir = area.substr(k);
  if (area.size()>0) { area = area+"/";}
  string prefix = pars.getPrefix();
  if (prefix.size()>0) { prefix = prefix+".";}
  
  string sn1 = setName;
  
  // dont need redundant setName if subdirectory is already the setName
  //if (sn1==subdir) sn1="";
  
  string fn1 = fileName;
  size_t f1 = fn1.find(sn1); 
  if ((fn1.size()>0)&&(f1==string::npos)) sn1=fn1+"."+sn1;    
  if (sn1.size()>0) { sn1 = sn1+".";}


  for (iterContig = contig.begin();	iterContig != contig.end(); iterContig++) {
      string cn1 = iterContig->first;      
      string basename = area+prefix+sn1+cn1;
      if (setName==cn1) basename = area+prefix+cn1;
    
      // replace evil character "|" with benign "_" 
      size_t found=basename.find("|");
      while (found!=string::npos) {
        basename.replace(found,1,"_");
        found=basename.find("|");
      }
      
      if (contig[cn1].Length>0) { 
        cout << "write output for read depth " << cn1 << endl;
        string fname = basename+".read.depth.span";
        cout << "\t " << fname << "\t " << contig[cn1].read_depth.n.size() << endl;
        contig[cn1].writeDepth(fname, contig[cn1].read_depth);
        cout << "write output for read start " << cn1 << endl;
        fname = basename+".read.start.span";
        cout << "\t " << fname << "\t " << contig[cn1].read_start.n.size() << endl;
        contig[cn1].writeDepth(fname, contig[cn1].read_start);
        cout << "write output for frag depth " << cn1 << endl;
        fname = basename+".frag.depth.span";
        cout << "\t " << fname << "\t " << contig[cn1].frag_depth.n.size() << endl;
        contig[cn1].writeDepth(fname, contig[cn1].frag_depth);
       }
  }
}

void C_set::writeCountReads() {
  C_contigs::const_iterator iterContig;   
  string area = pars.getOutputDir();
  int k = 1+area.find_last_of("/");
  string subdir = area.substr(k);
  if (area.size()>0) { area = area+"/";}
  string prefix = pars.getPrefix();
  if (prefix.size()>0) { prefix = prefix+".";}
  
  // input filename + setname 
  string sn1 = setName;
  // dont need redundant setName if subdirectory is already the setName
  // if (sn1==subdir) sn1="";
  string fn1 = fileName;
  size_t f1 = fn1.find(sn1); 
  if ((fn1.size()>0)&&(f1==string::npos)) sn1=fn1+"."+sn1;   
  if (sn1.size()>0) { sn1 = sn1+".";}


  for (iterContig = contig.begin();	iterContig != contig.end(); iterContig++) {
      string cn1 = iterContig->first;      
      string basename = area+prefix+sn1+cn1;
      if (setName==cn1) basename = area+prefix+cn1;

      // replace evil character "|" with benign "_" 
      size_t found=basename.find("|");
      while (found!=string::npos) {
        basename.replace(found,1,"_");
        found=basename.find("|");
      }
      
      if (contig[cn1].Length>0) { 
        cout << "write output for read depth " << cn1 << endl;
        string fname = basename+".read.count.span";
        cout << "\t " << fname << "\t " << contig[cn1].read_counts.n.size() << endl;
        contig[cn1].writeDepth(fname, contig[cn1].read_counts);
       }
  }
}

//------------------------------------------------------------------------------
// code for paired-end constraint  + put closest map pos first in M map list
// unless M is an element,then put the longest L element in M map list. 
//------------------------------------------------------------------------------
C_pairedread C_set::resolvePairConstraint(C_pairedread & pair0) {
    int Nalign0 = pair0.read[0].align.size();
    int Nalign1 = pair0.read[1].align.size();
    
    //--------------------------------------------------------------------------
    // one library in build phase.... get from pars  
	  // this is a problem for bam files with multiple libraries 
    //--------------------------------------------------------------------------
    int lmLow = int(pars.getFragmentLengthLo());
    int lmHigh = int(pars.getFragmentLengthHi());    
    double fraglen = pars.getFragmentLength();    
 
    //--------------------------------------------------------------------------
    // widen resolve window for first M map - to LF 'ballpark' 
    //--------------------------------------------------------------------------
    lmHigh=(2*lmHigh)-int(fraglen);
    lmLow=(2*lmLow)-int(fraglen);

    int i0,i1;
    
    // blank output uu constrained pair
    C_pairedread pair1;
    
    if ((Nalign0*Nalign1)==0) {
      return pair0;                     // missing end - singleton
    } else if ((Nalign0*Nalign1)==1) {
      return pair0;                     // already unique - no constraint needed
    } else {
      //------------------------------------------------------------------------
      // check combinations for a single combination meeting the constraint
      //------------------------------------------------------------------------
      int Nok = 0;
      // init lm closest to average lm
      double dlmin=1e20;
      // init indices to best lm
      int    m0 = -1;
      int    m1 = -1;
      int   m0p =  0;
      int   m1p =  0;
      // init indices to longest element
      int    mm0 = -1;
      int    mm1 = -1;
      int   mmL0 =  0;
      int   mmL1 =  0;
      // int target lm
      //double fraglen=0.5*double(lmLow+lmHigh);      
      // loop over end 0
      for (i0=0; i0<Nalign0; i0++) {
          unsigned char c0 = pair0.read[0].align[i0].anchor;
          if (c0>=elementMinAnchor) { 
            string aname= anchors.names[c0];
            int mmlen= anchors.L[aname]; 
            if (mmlen>mmL0) {
              mmL0=mmlen;
              mm0=i0;
            }
            //continue; 
          }
          unsigned int p0 = pair0.read[0].align[i0].pos;
          unsigned short l0 = pair0.read[0].align[i0].len;
          unsigned char o0 = pair0.read[0].align[i0].sense;
          
          /* missing constraint debug
          if ( (p0>29837)&&(p0<29897)&&(Nalign0==1)&&(Nalign1>1) ) {
             cerr << p0 << endl;
          }
          */

          
          // loop over end 1
          for (i1=0; i1<Nalign1; i1++) {
            double lm = 1e30; 
            unsigned char c1 = pair0.read[1].align[i1].anchor;
            if (c1>=elementMinAnchor) { 
              string aname= anchors.names[c1];
              int mmlen= anchors.L[aname]; 
              if (mmlen>mmL1) {
                mmL1=mmlen;
                mm1=i1;
              }
              //continue; 
            }
            unsigned int p1 = pair0.read[1].align[i1].pos;
            unsigned short l1 = pair0.read[1].align[i1].len;
            unsigned char o1 = pair0.read[1].align[i1].sense;
            if ((o0=='F')&&(o1=='R')&&(c0==c1)) {
                lm = p1+l1-p0;
            } else if ((o0=='R')&&(o1=='F')&&(c0==c1)) {
                lm = p0+l0-p1;
            }
            
            // fragment constrained
            bool ok=(lm>lmLow)&&(lm<lmHigh); 
            if (ok) {
               Nok++;
               // make a pair record with the first constrained maps
               if (Nok==1) {
                  pair1.read[0].align.push_back(pair0.read[0].align[i0]);
                  pair1.read[0].element=pair0.read[0].element;
                  pair1.read[1].align.push_back(pair0.read[1].align[i1]);
                  pair1.read[1].element=pair0.read[1].element;
               }
            }
            // check for closest lm to fraglen
            double dl = fabs(lm-fraglen);
            if (dl<dlmin) {
               m0=i0;
               m1=i1;
               m0p=p0;
               m1p=p1;
               dlmin=dl;
            }
          }
          // if (Nok>1) break;
          // get out of loop if more than one combination is ok and this is MM fragment
          // keep going if this is a UM fragment
          if ( (Nalign0>1)&&(Nalign1>1)&&(Nok>1) ) break;
      }
      
      // contrain collapes UM to UU or MM to UU
      // bool MM = (Nalign0>1)&&(Nalign1>1);
      // bool UM = ((Nalign0==1)||(Nalign1==1))&&( (Nalign0*Nalign1)>1) ;
      bool noelem=(pair0.read[0].element==0)&&(pair0.read[1].element==0);
      bool crunch = (Nok==1) && noelem;   // no elements allowed in constrained UU
      
      /*
      // debug alu constraint
      if ((!noelem)&&( (Nalign1==1)||(Nalign0==1) ) ) {
        if (((pair0.read[0].element&1)>0)&&(Nalign1==1)&&(Nalign0>1)&&(dlmin>300.0)) {
             cerr << "alu fail" << endl;
        }
        if (((pair0.read[1].element&1)>0)&&(Nalign0==1)&&(Nalign1>1)&&(dlmin>300.0)) {
             cerr << "alu fail" << endl;
        }
      }
      */
      //------------------------------------------------------------------------
      // never crunch - just change order...
      //------------------------------------------------------------------------
      crunch = false;
                  
      // crunch constrained fragments
      if ( crunch ){
        pair1.read[0].Nalign=1;
        pair1.read[1].Nalign=1;
        return pair1;        
      } else { 
         // complain if too many combinations are within contraint
         if (Nok>200) {
            cerr << " more than 200 combinations contrained in resolvePairContraint! " << endl;
         }

         //---------------------------------------------------------------------
         // mobile element hit - put largest reference hit first if not resolved
         //---------------------------------------------------------------------
         if (((mm1+mm0)>-2)&&(Nok==0))  {
            if (mm0>=0) {
                pair1=pair0;
                C_readmap a1=pair1.read[0].align[mm0];
                pair1.read[0].align.erase(pair1.read[0].align.begin()+mm0);
                pair1.read[0].align.insert(pair1.read[0].align.begin(),a1);
            }
            if (mm1>=0) {
                pair1=pair0;
                C_readmap a1=pair1.read[1].align[mm1];
                pair1.read[1].align.erase(pair1.read[1].align.begin()+mm1);
                pair1.read[1].align.insert(pair1.read[1].align.begin(),a1);
            }         
            return pair1;
         }

         //---------------------------------------------------------------------
         // no best lm to put at begining of alignment vector
         //---------------------------------------------------------------------         
         if (m0<0) {
            return pair0;
         } else {
            //------------------------------------------------------------------
            // swap m0 and m1 alignments to the front of the vector...
            //------------------------------------------------------------------
            pair1=pair0;
            C_readmap a1=pair1.read[0].align[m0];
            // debug
            if (int(a1.pos)!=m0p) {
                cerr << m0p << endl;
            }    
            pair1.read[0].align.erase(pair1.read[0].align.begin()+m0);
            pair1.read[0].align.insert(pair1.read[0].align.begin(),a1);
            a1=pair1.read[1].align[m1];
            // debug
            if (int(a1.pos)!=m1p) {
                cerr << m1p << endl;
            }    
            pair1.read[1].align.erase(pair1.read[1].align.begin()+m1);
            pair1.read[1].align.insert(pair1.read[1].align.begin(),a1);
            return pair1;
         }
      } 
    }
}

/*
C_pairedread C_set::resolvePairConstraint(C_pairedread & pair0) {
    int Nalign0 = pair0.read[0].align.size();
    int Nalign1 = pair0.read[1].align.size();
    double lm0 = pars.getFragmentLength();
    int lmLow = int(pars.getFragmentLengthLo());
    int lmHigh = int(pars.getFragmentLengthHi());    
    int i0,i1;
    // blank output object
    C_pairedread pair1;
    if ((Nalign0*Nalign1)==0) {
      return pair0;  // missing end - singleton
    } else if ((Nalign0*Nalign1)==1) {
      return pair0;  // already unqiue - no constraint needed
    } else {
      //------------------------------------------------------------------------
      // check combinations for a single combination meeting the constraint
      //------------------------------------------------------------------------
      bool reduce = false;
      struct tfit {
        double d ;
        int lm;
        int i0;
        int i1;
      };
      double big = 1.0e30;
      tfit best = { big, -1,-1}; 
      tfit next = { big, -1,-1};
      for (i0=0; i0<Nalign0; i0++) {
          unsigned char c0 = pair0.read[0].align[i0].anchor;
          unsigned int p0 = pair0.read[0].align[i0].pos;
          unsigned short l0 = pair0.read[0].align[i0].len;
          unsigned char o0 = pair0.read[0].align[i0].sense;
          for (i1=0; i1<Nalign1; i1++) {
            double lm = big*2.0; 
            unsigned char c1 = pair0.read[1].align[i1].anchor;
            unsigned int p1 = pair0.read[1].align[i1].pos;
            unsigned short l1 = pair0.read[1].align[i1].len;
            unsigned char o1 = pair0.read[1].align[i1].sense;
            if ((o0=='F')&(o1=='R')&(c0==c1)) {
                lm = p1+l1-p0;
            } else if ((o0=='R')&(o1=='F')&(c0==c1)) {
                lm = p0+l0-p1;
            }
            double d=sqrt(pow(lm-lm0,2));
            if (c0==c1) {
              if (d<best.d) {
                next.d=best.d;
                next.lm=best.lm;                
                next.i0=best.i0;
                next.i1=best.i1;
                best.d=d;
                best.lm=int(lm);                
                best.i0=i0;
                best.i1=i1;
              } else if (d<next.d) {
                next.d=d;
                next.lm=int(lm);
                next.i0=i0;
                next.i1=i1;
              }
            }
          }            
      }
      if (best.d<big) {        
        bool ok1=(best.lm>lmLow)&(best.lm<lmHigh); 
        bool ok2=(next.lm>lmLow)&(next.lm<lmHigh); 
        if (ok1 &(!ok2)) { // constraint 
          pair1.read[0].align.push_back(pair0.read[0].align[best.i0]);
          pair1.read[1].align.push_back(pair0.read[1].align[best.i1]);
          reduce = true;
        }
      }
      if (reduce) {
        return pair1;
      } else { 
        return pair0;
      } 
    }
}
*/
ostream &operator<<(ostream &output, const C_set & set1)
{
   output << set1.getSetName()  << ":Nfrag=" << set1.Nfrag  << ",Npair=" << set1.Npair << "\t " ;
   return output;
}

// fetch set name
string C_set::getSetName() const
{
  return this->setName;
}

// set set name vector
void C_set::setSetName(string & s1)
{
  this->setName = s1;
  if (s1[s1.size()-1]=='.') {
    this->setName.erase(s1.size()-1);
  }
}

// fetch file name
string C_set::getFileName() const
{
  return this->fileName;
}
// set file name 
void C_set::setFileName(string & s1)
{
  this->fileName = s1;
}

//------------------------------------------------------------------------------
 // fill in redundant read part of library info
//------------------------------------------------------------------------------
void C_set::calcLibraryRedundancy(C_libraryinfo & lib1) 
{
  C_contigs::const_iterator iterContig;   
  lib1.NPair=0;
  lib1.NSingle=0;
  lib1.NPairRedundant=0;
  lib1.NSingleRedundant=0;
  
  // cout << "calculate fragment redundancy " << endl;
  
  for (iterContig = contig.begin();	iterContig != contig.end(); iterContig++) {
      string name = iterContig->first;
      if (contig[name].Length>0) { 
        
        cout << "\t sort contig " << name << endl;
        // fill depth & repeat vector before uniquifying pairs
        contig[name].calcLengths();
        C_contig c1 = contig[name];
        lib1.NPair+=c1.localpairs.size();
        lib1.NPair+=c1.crosspairs.size();
        lib1.NSingle+=c1.dangle.size();
        lib1.NSingle+=c1.singleton.size();
        lib1.NSingle+=c1.umpairs.size();        
        c1.sort();
        c1.uniquify();
        lib1.NPairRedundant+=c1.localpairs.size();
        lib1.NPairRedundant+=c1.crosspairs.size();
        lib1.NSingleRedundant+=c1.dangle.size();
        lib1.NSingleRedundant+=c1.singleton.size();
        lib1.NSingleRedundant+=c1.umpairs.size();        
      }
    }
    // count the redundant fragments rather than the unique ones....
    lib1.NPairRedundant=lib1.NPair-lib1.NPairRedundant;
    lib1.NSingleRedundant=lib1.NSingle-lib1.NSingleRedundant;
}

/*
//------------------------------------------------------------------------------
// remove redundant fragments arising from Bam format conventions and 
// fill in redundant read part of library info
//------------------------------------------------------------------------------
void C_set::fixBamRedundancy(C_libraryinfo & lib1) 
{
  C_contigs::const_iterator iterContig;   
  lib1.NPair=0;
  lib1.NSingle=0;
  lib1.NPairRedundant=0;
  lib1.NSingleRedundant=0;
  
  cout << "fix Bam record fragment redundancy " << endl;
  
  for (iterContig = contig.begin();	iterContig != contig.end(); iterContig++) {
      string name = iterContig->first;
      if (contig[name].Length>0) { 
        
        cout << "\t sort contig " << name << endl;
        // fill depth & repeat vector before uniquifying pairs
        contig[name].calcLengths();
        lib1.NPair+=contig[name].localpairs.size();
        lib1.NPair+=contig[name].crosspairs.size();
        lib1.NSingle+=contig[name].dangle.size();
        lib1.NSingle+=contig[name].singleton.size();
        lib1.NSingle+=contig[name].umpairs.size();        
        contig[name].sort();
        contig[name].uniquifyBam();
        lib1.NPairRedundant+=contig[name].localpairs.size();
        lib1.NPairRedundant+=contig[name].crosspairs.size();
        lib1.NSingleRedundant+=contig[name].dangle.size();
        lib1.NSingleRedundant+=contig[name].singleton.size();
        lib1.NSingleRedundant+=contig[name].umpairs.size();        
      }
    }
    // count the redundant fragments rather than the unique ones....
    lib1.NPairRedundant=lib1.NPair-lib1.NPairRedundant;
    lib1.NSingleRedundant=lib1.NSingle-lib1.NSingleRedundant;
}


 */

//------------------------------------------------------------------------------
// PairedFilesClass is the container class for paired-end reads files
//------------------------------------------------------------------------------
C_pairedfiles::C_pairedfiles( string  & file1, RunControlParameters & pars)
{
	time(&tprev);         // init clock
	Nbad=0;
	NbadPos=0;
	MateMode=0;
	string file2="";
	this->inputcheck(file1,pars);
	if (this->inputType=='S') {
		this->loadSpanner(file1,pars);
	} else if ((this->inputType=='Z')||(this->inputType=='B'))  {
		this->loadBam(file1,pars);
		//this->testMultiMapBam(file1,pars);
		//this->loadBamSortedByName(file1,pars);
	} else if (this->inputType=='D') {
		this->loadSpannerDirectory(file1,pars);
	} else {
		cerr << "inputType \t" << this->inputType << "  not implemented " << endl;
	}
	// summaryLog();
	if (pars.getDoScanOnly()) {
		printScanResults();
	}  
}

float C_pairedfiles::elapsedtime() {
  time_t t;
  time(&t);
  float et = difftime (t,tprev);
  //float et=float(t-tprev)*1.0/float(CLOCKS_PER_SEC);
  return et;
}


void C_pairedfiles::summaryLog() {
  // loop over existing sets
  cout << endl;
  for(int s=0; s < int(this->setNames.size()); s++) {
    cout << set[s] << endl;
    //cout << "read length" << endl;
    //cout << set[s].lengthStats << endl;   
	HistObj h = set[s].lengthStats.h.collapse(20);
	cout << h << endl;
    //cout << set[s].lengthStats.h << endl;   
    //cout << "alignments per read" << endl;
    //cout << set[s].repeatStats << endl;   
    h = set[s].repeatStats.h.collapse(20);
    cout << h << endl;
    //cout << "fragment length" << endl;
    //cout << set[s].fragStats << endl;
    h = set[s].fragStats.h.collapse(20);
    cout << h << endl;
    //cout << "map quality" << endl;
    //cout << set[s].qStats << endl;   
    cout << set[s].qStats.h << endl;   
    //cout << "pair stats" << endl;
    //cout << set[s].pairStats << endl;   
    cout << set[s].pairCountStats.h << endl;   
	  
	cout << set[s].pairModelStats.h << endl;   

  }   
}

//------------------------------------------------------------------------------
// print scan result file
//------------------------------------------------------------------------------
void C_pairedfiles::printScanResults() {
	
  string filename = set[0].getFileName();
	
  string outfilename=set[0].pars.getOutputDir()+"/"+set[0].pars.getPrefix()
	+filename+".stats";
  ofstream outfile(outfilename.c_str(), ios::out);
  if (!outfile) {
		cerr << "unable to open output stats file: " << outfilename << endl;
		return;
  }
	
	
  for(int s=0; s < int(this->setNames.size()); s++) {
	  string setname = set[s].getSetName();
	  int tech = 0;
	  string ReadGroupID="";
	  string sample="";
	  C_librarymap::iterator it;
	  //int N = set[s].libraries.libmap.size();	 
	  for ( it=set[s].libraries.libmap.begin() ; it != set[s].libraries.libmap.end(); it++ )
	  {
		  //unsigned int ReadGroupCode=(*it).first;
		  C_libraryinfo lib1 = (*it).second; 		  
		  ReadGroupID = lib1.Info.ReadGroupID;
		  sample = lib1.Info.SampleName;
		  tech  = lib1.Info.SequencingTechnology;
		  // always get first read group
		  //if (setname.compare(ReadGroupID) == 0) break;
		  string PlatformTech="unknown";
		  switch (tech) {
			  case 1:
				  PlatformTech = "454";
				  break;
			  case 4:
				  PlatformTech = "Illumina";
				  break;
			  case 8:
				  PlatformTech = "PacificBiosciences";
				  break;
			  case 16:
				  PlatformTech = "SOLiD";
				  break;
		  }

		  outfile << "@RG\tID:" << setname << "\tSM:" << sample <<"\tPL:" << PlatformTech << endl;
      }
	  
	  outfile << "\n";

	  outfile << set[s].fragStats.h << "\n";
	  //outfile << "read length histogram\n";
	  outfile << set[s].lengthStats.h << "\n";
	  //outfile << "repeat multiplicity histogram\n";
	  outfile << set[s].repeatStats.h << "\n";
	  //outfile << "mapping quality histogram\n";
	  outfile << set[s].qStats.h << "\n";
	  //outfile << "pair stat histogram\n";
	  outfile << set[s].pairCountStats.h << "\n";

	  outfile << set[s].pairModelStats.h << "\n";

	  outfile << set[s].spanStats.h << "\n";
	  
	  outfile << set[s].refStats.h << "\n";

  }
  outfile.close(); 
}


//==================
// check input type
//==================
void C_pairedfiles::inputcheck( string & in1, RunControlParameters & pars)
{
	// Assembly info line: for ACE file
	//string patternAS("^AS\\s+(\\d+)\\s+(\\d+)");
	// fasta header line AB map file
	//boost::regex patternAB("^>\\d+\\_\\d+\\_\\d+\\_");
	//string patternAB("^>\\S+\\_[F|R]3");
	// fasta header line Mosiak RAR map file
	string patternMSK("^MSKAA");
	// fasta header line Helicos csv map file
	//string patternHS("^Reference_ID");
	// attempt parsing set of Spanner files if there is no in2 filename
	fstream file1;
  
  
	C_headerSpan h1;
	if (checkSpannerDirectory(in1))  {
          inputType='D';  // spanner files      
          return;
	}
	
	if (checkBamDirectory(in1))  {
		inputType='B';  // spanner files      
		return;
	}
	
	if (checkSpannerFiles(in1)) {
		//if (headers.size()!=6) return; 
		if (headers.size()!=4) return; 
         bool ok = true;
         C_headers::iterator i;
         for(i=headers.begin(); i != headers.end(); ++i) {
            string filename = i->first;
            h1 = i->second;
            //cout  << filename << " \t " << h1 << endl;
            ok = ok && ((h1.V>=201)|(h1.V<=207));
         }
         if (ok) { 
            // if all Spanner files are ok, then put the setName into vector as element 0
            setNames.clear();
            setNames.push_back(h1.setName);
            // input type "S" for Spanner....
            inputType='S';  // spanner files      
            return;
         }
      }
      //==================================
      // must be a file of some sort
      // check if file exists
      //==================================
      file1.open(in1.c_str(), ios::in);	
      if (!file1) {
        cerr << "Unable to open file: " << in1 << endl;
        exit(1);
      }
      file1.close();
  
  //==================================
  // check non-Spanner input files....
  //==================================
  file1.open(in1.c_str(), ios::in);	
  if (!file1) {
    cerr << "Unable to open file: " << in1 << endl;
    exit(1);
  }
  string line;
  string match;
  inputType='X';
  int n = 0;
  while ( (getline(file1, line).good()) && (inputType=='X') && (n<4) ) {		
	//================
    // line test
	//================
	if (RE2::FullMatch(line.c_str(),patternMSK.c_str(),&match) ) {
        inputType='M'; 
    };	
    n++;	
    // skip comment lines in header (Helicos...)
    if (line[0]=='#') {
      n--;
    }
	};
  file1.close();
  if (inputType=='X') {
    if (checkBamFile(in1))  {
        inputType='Z';  // Bam files      
    }
  }
}

//------------------------------------------------------------------------------
// extract read mapping info from native Spanner files
//------------------------------------------------------------------------------
void C_pairedfiles::loadSpannerDirectory( string  & dir1,  RunControlParameters & pars)
{
  C_headers::iterator hi;
 
  int Nfile = SpannerFileNames.size();
  for (int i=0; i<Nfile; i++) {
     string f = dir1+"/"+SpannerFileNames[i];
     headers.clear();
     if (checkSpannerFiles(f)) {
        hi=headers.begin();
        string filename = hi->first;
        int k = 1+dir1.find_last_of("/");
        string sets = dir1.substr(k);
        headers[filename].setName=sets;
        loadSpanner(f, pars);
     }
  }
  cout << "-----------------------------------------------------------------\n";
  //cout << "final stats: "<< dir1 << endl;
  //summaryLog();
  //cout << "-----------------------------------------------------------------\n";

}

//------------------------------------------------------------------------------
// extract read mapping info from native Spanner files - one contig at a time
//------------------------------------------------------------------------------
void C_pairedfiles::loadSpanner( string  & file1,  RunControlParameters & pars)
{
  C_headers::iterator hi;
  if (headers.size()==0) {
    cerr << " bad file " << file1 << endl;
    return ;
  }
  hi=headers.begin();
  string filename = hi->first;
  C_headerSpan h1 = hi->second;
  string setname = h1.setName;
  string contigname = h1.contigName;
  if (set.size()==0) {
    C_set set1(setname, pars); 
    set.push_back(set1);
  }
  // initialize C_contig c1;
  set[0].contig[contigname];
  // copy anchor info to set and contig objects
  set[0].anchors=anchors;
  set[0].contig[contigname].anchors=anchors;
  // setName in contig should be setName in set
  set[0].contig[contigname].setName=set[0].getSetName();
  //if (headers.size()!=4) return;    
	if (headers.size()!=4) return;  

  // libraries loaded in checkSpannerFiles
  set[0].libraries=libraries;
   
  //----------------------------------------------------------------------------
  // if tailcut not pars.getFragmentLengthWindow() then revise LMlow, LMhigh 
  //----------------------------------------------------------------------------
  double FragmentLengthWindow = pars.getFragmentLengthWindow();
  int lmmax = set[0].libraries.maxLF();
  if (lmmax>0) { 
    double tailcut = set[0].libraries.getTailcut();
    if (fabs(tailcut-FragmentLengthWindow)>1e-5)  {
       set[0].libraries.resetFragLimits(FragmentLengthWindow);
    }
  }

  
   
  for(hi=headers.begin(); hi != headers.end(); ++hi) {
      filename = hi->first;
      h1 = hi->second;
      cout  << filename << " \t " << endl;
      if (h1.V<203) break; 
      switch (h1.type)  {
        case 0: 
          set[0].contig[contigname].loadPairs(filename);
          break;
        case 1:
          set[0].contig[contigname].loadCross(filename);
          break;
        case 100:
          set[0].contig[contigname].loadSingleton(filename);
          break;
        case 2:
          set[0].contig[contigname].loadDangle(filename);
          break;
        case 3:
          set[0].contig[contigname].loadMultipairs(filename);
          break;
        case 101:
          set[0].contig[contigname].loadRepeat(filename);
          break;
      }
  }
  //---------------------------------------------------------------------------
   set[0].contig[contigname].setLengthFromAnchor(); 
  //---------------------------------------------------------------------------
  // Calculate average & stdev read lengths for this contig
  //---------------------------------------------------------------------------
  set[0].contig[contigname].calcLengths();
  //---------------------------------------------------------------------------
  // make sure contig set name is same as set name (loadSpannerDir)
  //---------------------------------------------------------------------------
  set[0].contig[contigname].setName=set[0].getSetName();
  //---------------------------------------------------------------------------
  // if mask file exists with contigname - add mask info to repeat Marker object
  //---------------------------------------------------------------------------
  string ReferenceFastaFilefile = pars.getReferenceFastaFile(); 

  if ( (ReferenceFastaFilefile.size()>0) and (pars.getDoMasking()) ) {
    set[0].contig[contigname].repeat.addMarks(ReferenceFastaFilefile, contigname,'N'); 
  }
}

//------------------------------------------------------------------------------
// extract read mapping info from BAM file - must be sorted by read name 
//------------------------------------------------------------------------------
void C_pairedfiles::loadBam( string  & file1, RunControlParameters & pars)
{
	
	// Bam file pattern - bam file shgould have bam extenstion
	string patternFile("^.*/(.+)(\\.BAM|\\.bam|\\.Bam)");
	// get allowed contig name regexp
	string AllowContigsRegex=pars.getAllowContigsRegex();
	// regex match
	string match;
	
	// Debug level 
	int dbg =pars.getDbg();  
	// flip relative to Illumina short PE FR convention
	MateMode=pars.getReadPairSenseConfig();
	// Scan only? 
	bool scan =pars.getDoScanOnly();
	if (scan) {
		SpannerMode=SPANNER_SCAN;
	} else {
		SpannerMode=SPANNER_BUILD;
	}

	vector<string> elements =pars.getMobileElements();
	bool moblist=elements.size()>0;
	if ( (!scan) & moblist ) {
		SpannerMode=SPANNER_MOBI;
	}
		
	// get Max Fragments
	MaxFragments  = pars.getMaxFragments();    
	// get Minimum Q mapping value for unique read
	Qmin  = pars.getQmin();
	// get Bam access mode (0=random access, 1=ZA tag, 2=sortedbyReadName)
	BamZ  = pars.getBamZA();
	
	// declare some ubiquitous stringy vars
	string s,s1,s2;
	// declare   pair read objects
	C_pairedread pair1, pair2;
	// default set name
	string name="BAM";
	if (RE2::FullMatch(file1.c_str(),patternFile.c_str(),&match) ) {
		name=match;
	} else {
		vector<string> parts;
		split(parts,file1,"/");
		name=parts[parts.size()-1];
	}
	
	
	// ReadGroups
	unsigned int ReadGroupCode = 0;

	//---------------------------------------------------------------------------
	// Check if this is a valid BAM alignment file (optional)
	//---------------------------------------------------------------------------
	if (checkBamFile(file1))  {	  
		BamFileNames.clear();
		BamFileNames.push_back(file1); 
	} else if (checkBamDirectory(file1))  {	  
		cout << "opening " << BamFileNames.size() << " Bam files for input" << endl;     
	} else  {
		cout << "Error opening " << file1 << " for input" << endl;     
		exit(101);
	}
	
	//---------------------------------------------------------------------------
	// open BAM  file(s)
	//---------------------------------------------------------------------------
	
	BamMultiReader ar1,ar2;	
	BamReader br1,br2;
	
	BamAlignment ba;
		
	ar1.Open(BamFileNames,BamZ<1,false);
	
	if (ar1.GetReferenceCount()<1) {
		cout << "ERROR: Unable to open the BAM file (" << file1.c_str()<< ")." << endl;
		exit(102);
		
	}  else {
		if (!ar1.GetNextAlignment(ba)) {
			cerr << " ERROR: Unable to load record BAM file  (" << file1.c_str() << ")." << endl;
			exit(103);
		}
		ar1.Rewind(); 
	}
	
	if ( (BamZ<1) ) {
		
		ar2.Open(BamFileNames,true,false);
		if (ar2.GetReferenceCount()<1) {
			cout << "ERROR: Unable to open the BAM file twice (" << file1.c_str() << ")." << endl;
			exit(104);
		}  else {
			if (!ar2.GetNextAlignment(ba)) {
				cerr << " ERROR: Unable to load record BAM file  (" << file1.c_str() << ")." << endl;
				exit(105);
			}
			ar2.Rewind(); 
		}
		
  }
	
	// retrieve the header text BAM 
	string samHeader = ar1.GetHeaderText();	
		
	// fill anchorinfo from Bam info  (don't use NT_ anchors)
	C_anchorinfo anchor1(ar1); 

	// override anchors with anchorfile? 
	string anchorfile1=pars.getAnchorfile();
	if (anchorfile1.size()>0) {
		C_anchorinfo anchorf(anchorfile1);
		anchor1=anchorf;
	}
	if (scan) {
		AllowContigsRegex = "doScanOnlyDoesNotUseContigs!";		
	}
	
	// only process specified contigs
	anchor1.anchorlimit(AllowContigsRegex);
	// bookeep mobile element anchors default
	anchor1.anchorElements(elements);
	// set this object anchors
	anchors = anchor1;
	
	//---------------------------------------------------------------------------
	// fetch library info from BAM alignments file
	//---------------------------------------------------------------------------
	C_libraries libs(ar1);
	if (libs.libmap.size()<1) {
		cerr << "no ReadGroups" << endl;
		exit(-1);
	}
	libraries = libs;

	//---------------------------------------------------------------------------
	// check first library for platform info. set flip if 454 or AB
	//---------------------------------------------------------------------------
	vector<char> seqtech = libraries.getSequencingTechnology();
	char seqtech1 = seqtech[0];
	switch (seqtech1) {
		case ST_454:
			MateMode=MATEMODE_454;
			break;
		case ST_SOLID:
			MateMode=MATEMODE_SOLID;
			break;
		default:
			MateMode=MATEMODE_ILLUMINA;
			break;
	}
	
	pars.setReadPairSenseConfig(MateMode);
	
	//---------------------------------------------------------------------------
	// make multiple sets only for Scan, one set for build, detect *(or 454 scan) 
	//---------------------------------------------------------------------------
	int Nset=0;
	if ((scan)&&(seqtech1!=ST_454)) { 
		Nset=makeSetsFromLibs(name,pars,anchors);
	} else {
		Nset=makeOneSetFromLibs(name,pars,anchors);
	}
	
	if (Nset<1) {
		cerr << " problem getting libraries " << endl;
		exit(-1);
	}
	
	//----------------------------------------------------------------------------
	// read mapping histograms
	//----------------------------------------------------------------------------
	int iset;
	for (iset=0; iset<Nset; iset++) {
		set[iset].initContigs();
	}
	
	
	cout << "-----------------------------------------------------------" <<   endl;
	cout << " load " << file1 <<   endl;
	cout << "-----------------------------------------------------------" <<   endl;
	
	// counter for uu frags
	int Nuu[2]={0,0};
	
	// fragment counter
	int Nfrag = 0;
	
	// clear doneFrag map
	doneFrag.clear();

	if ( BamFileNames.size()==1 ) {
		br1.Open(file1);
		br2.Open(file1,"",true);
	}
		
	bool next = true;
	//---------------------------------------------------------------------------
	// loop through file to load non-proper pairs 
	//---------------------------------------------------------------------------
	while(next) { 
		//if (BamFileNames.size()>1) { 
		//	next = nextBamAlignmentPair(ar1,ar2,mr)&&(Nfrag<=MaxFragments);
		//} else if (BamFileNames.size()==1) {
			next= nextBamAlignmentPair(br1,br2,pair1)&&(Nfrag<=MaxFragments);
		//}
		if (!next) continue;
		
		Nfrag++;
		
		
		// convert mosaik struct to Spanner pair object
		//pair1 = Mosaik2pair(mr);
		
		// ReadGroupCode
		ReadGroupCode = pair1.ReadGroupCode;
		iset=ReadGroupCode2set[ReadGroupCode];
		
		set[iset].Nfrag++;
		
		//--------------------------------------------------------------------------
		// occasional status report at intervals of dbg
		//--------------------------------------------------------------------------
		if ( (dbg!=0)&&(elapsedtime()>float(dbg))) {
			time(&tprev);
			cout << Nfrag << "\t";
			for (int iset1=0; iset1<Nset; iset1++) {
				cout << set[iset1];
			}
			cout << endl;
			
		}      
		//------------------------------------------------------------------------
		// this is a complete pair
		//------------------------------------------------------------------------
		/*
		if (scan) {
			pair2=pair1;
		} else {
			// worry about handling multiple sets here XXXXXXXXXXXXX ???? XXXXXXX
			pair2 = set[iset].resolvePairConstraint(pair1);  
		}
		*/
		
		
		set[iset].repeatStats.Fill1(int(pair1.read[0].Nalign));
		set[iset].repeatStats.Fill1(int(pair1.read[1].Nalign));
		
		bool bothends = (pair1.read[0].Nalign*pair1.read[1].Nalign)>0;
		
		bool uu = (pair1.read[0].Nalign*pair1.read[1].Nalign)==1;
		char constrain = 0; 
		if (bothends) {
			
			set[iset].Npair++;
			
			// mapping quality
			set[iset].qStats.Fill1(pair1.read[0].align[0].q);         
			set[iset].qStats.Fill1(pair1.read[1].align[0].q);         
			
			if (uu &&  (pair1.read[0].align[0].q>=Qmin)&&(pair1.read[1].align[0].q>=Qmin) ) { 
				//----------------------------------------------------------------------
				// make some stats with this unique mapped fragment  
				//----------------------------------------------------------------------
				int lm = set[iset].Fraglength(pair1);  
				int sl = set[iset].spanLength(pair1);
				if (lm!=-100000) {
					// make frag dist with high quality mapped pairs
					set[iset].fragStats.Fill1(lm);         
					set[iset].spanStats.Fill1(sl);         
					Nuu[1]++;  
				} else {  // abberant pairs
					Nuu[0]++;
				}
				
				int model = set[iset].pairModel(pair1);   
				set[iset].pairModelStats.Fill1(model);         
				set[iset].lengthStats.Fill1(pair1.read[0].align[0].len);         
				set[iset].lengthStats.Fill1(pair1.read[1].align[0].len);         
				set[iset].refStats.Fill1(int(pair1.read[0].align[0].anchor));         
				set[iset].refStats.Fill1(int(pair1.read[1].align[0].anchor));         
			}
			//------------------------------------------------------------------------
			// process pair for SV info
			//------------------------------------------------------------------------
			if (!scan) {
				set[iset].processPair(pair1,constrain);
			}
			//------------------------------------------------------------------------
			// count paired read alignments for set[0] pairstats
			//------------------------------------------------------------------------
			int N1[2]={int(pair1.read[0].align[0].nmap), int(pair1.read[1].align[0].nmap)};
			if (N1[0]*N1[1]==1) set[iset].Npair_11o++; 
			int n2[2]= {N1[0],N1[1]};
			if (n2[1]<n2[0]) { n2[1]=N1[0]; n2[0]=N1[1];}
			if (n2[0]>2) n2[0]=2;
			if (n2[1]>2) n2[1]=2;
			int npc = n2[0]*3 + n2[1];
			set[iset].pairCountStats.Fill1(npc);   // e is paired            
		} else { 
			//------------------------------------------------------------------------
			// not pair - dangling end
			//------------------------------------------------------------------------
			int np1 = pair1.read[0].Nalign;
			int np2 = pair1.read[1].Nalign;
			int np = (np1==0? np2: np1);
			if (np>2) { np=2; }
			set[iset].pairCountStats.Fill1(np);  
			set[iset].refStats.Fill1(int(pair1.read[0].align[0].anchor));         
			if (!scan) {
				// dangling missing ends
				set[iset].processPair(pair1,false);
			}      
		}
	}
	
	// finalize histos
	for (iset=0; iset<Nset; iset++) {
		set[iset].pairCountStats.Finalize();         
		set[iset].pairModelStats.Finalize();         
		set[iset].qStats.Finalize();         
		set[iset].fragStats.Finalize();         
		set[iset].spanStats.Finalize();         
		set[iset].lengthStats.Finalize();
		set[iset].repeatStats.Finalize();
		set[iset].refStats.Finalize();
		
		if (set[iset].libraries.libmap.size()!=1) {
			cerr << "==============================" << endl;
			cerr << " expected one library per set " << endl;
			cerr << set[iset].getSetName()  << " has " << set[iset].libraries.libmap.size() << endl;
			cerr << "==============================" << endl;
		}
		
		if (scan) {
			
			// finalize library: fragment length
			HistObj h=set[iset].fragStats.h;
			
			double lowf = (1.0-pars.getFragmentLengthWindow()/100.0)/2.0;
			double highf = (1.0-lowf); // half tail on high side  
			
			// store fraglength info from pars to library (override this data set)
			C_librarymap::iterator it;
			
			
			for ( it=set[iset].libraries.libmap.begin() ; it != set[iset].libraries.libmap.end(); it++ ) {
				ReadGroupCode  = (*it).first; 
				C_libraryinfo lib1 = (*it).second; 
				
				set[iset].libraries.libmap[ReadGroupCode].LM=int(h.mode);
				set[iset].libraries.libmap[ReadGroupCode].tailcut = pars.getFragmentLengthWindow();
				set[iset].libraries.libmap[ReadGroupCode].LMlow=int(h.p2xTrim(lowf)); // int(pars.getFragmentLengthLo());
				set[iset].libraries.libmap[ReadGroupCode].LMhigh=int(h.p2xTrim(highf)); //int(pars.getFragmentLengthHi()); 
				// store fraglength histogram from this data set to library
				set[iset].libraries.libmap[ReadGroupCode].fragHist=h;
				
				
				// Read length
				h=set[iset].lengthStats.h;
				set[iset].libraries.libmap[ReadGroupCode].readLengthHist=h;
				set[iset].libraries.libmap[ReadGroupCode].LR=int(set[0].lengthStats.h.mode);
				set[iset].libraries.libmap[ReadGroupCode].LRmax=int(h.xc[h.xc.size()-1]);
				int nb=h.n.size();
				int xmin= set[iset].libraries.libmap[ReadGroupCode].LRmax;
				for (int ib=nb; ib>0; ib--) {
					if ( (h.n[ib-1]>0)&&(h.xc[ib-1]<xmin) ) xmin=h.xc[ib-1];
				}
				set[iset].libraries.libmap[ReadGroupCode].LRmin=xmin;	  
				// Redundancy	  
				set[iset].calcLibraryRedundancy(set[iset].libraries.libmap[ReadGroupCode]);	  
			}
		}
	}
	cout << "-----------------------------------------------------------------\n";
	cout << "final stats: "<< file1 << endl;
	summaryLog();
	cout << "-----------------------------------------------------------------\n";
	
}

//------------------------------------------------------------------------------
// extract read mapping info from BAM file - must be sorted by read name 
//------------------------------------------------------------------------------
void C_pairedfiles::testMultiMapBam( string  & file1, RunControlParameters & pars)
{
	// Bam file pattern
	string patternFile("^.*/(.+)(\\.BAM|\\.bam|\\.Bam)");
	// get allowed contig name regexp
	string AllowContigsRegex=pars.getAllowContigsRegex();
	// regex match
	string match;
	// Debug level 
	int dbg =pars.getDbg();
	// Scan only? 
	bool scan =pars.getDoScanOnly();
	// get Max Fragments
	MaxFragments  = pars.getMaxFragments();
	
	// declare some ubiquitous stringy vars
	string s,s1,s2;
	// declare   pair read objects
	C_pairedread pair1;
	// set name
	string name="BAM";
	if (RE2::FullMatch(file1.c_str(),patternFile.c_str(),&match) ) {
		name=match;
	}	
	
	// ReadGroups
	//unsigned int ReadGroupCode = 0;
	
	//---------------------------------------------------------------------------
	// Check if this is a valid BAM alignment file (optional)
	//---------------------------------------------------------------------------
	
	if (checkBamFile(file1))  {	  
		BamFileNames.clear();
		BamFileNames.push_back(file1); 
	} else if (checkBamDirectory(file1))  {	  
		cout << "opening " << BamFileNames.size() << " Bam files for input" << endl;     
	} else  {
		cout << "Error opening " << file1 << " for input" << endl;     
		exit(106);
	}
	
	
	//---------------------------------------------------------------------------
	// open BAM  file
	//---------------------------------------------------------------------------
	BamMultiReader ar1;
	BamAlignment ba;
	
	// check reader1
	//ar1.Open(BamFileNames[0].c_str(),"",true,true);
	ar1.Open(BamFileNames,"",true,true);
	if (ar1.GetReferenceCount()<1) {
		cout << "ERROR: Unable to open the BAM file (" << file1.c_str() << ")." << endl;
		exit(106);
	}  else {
		if (!ar1.GetNextAlignment(ba)) {
			cerr << " ERROR: Unable to load record BAM file  (" << file1.c_str() << ")." << endl;
			exit(108);
		}
		ar1.Rewind(); 
	}
				
	// retrieve the header text BAM 
	string samHeader = ar1.GetHeaderText();
	
	//----------------------------------------------------------------------------
	// Bam files with multiple memory locations per fragment
	//----------------------------------------------------------------------------
	// pars.setDoRepeatCheck(true);
	
	// fill anchorinfo from Bam info  (don't use NT_ anchors)
	C_anchorinfo anchor1(ar1); 
	// override anchors with anchorfile? 
	string anchorfile1=pars.getAnchorfile();
	if (anchorfile1.size()>0) {
		C_anchorinfo anchorf(anchorfile1);
		anchor1=anchorf;
	}
	if (scan) {
		AllowContigsRegex = "doScanOnlyDoesNotUseContigs!";
	}
	
	// only process specified contigs
	anchor1.anchorlimit(AllowContigsRegex);
	// bookeep mobile element anchors default
	anchor1.anchorElements();

	anchors=anchor1;
	
	//---------------------------------------------------------------------------
	// fetch library info from BAM alignments file
	//---------------------------------------------------------------------------
	C_libraries libs(ar1);
	if (libs.libmap.size()<1) {
		cerr << "no ReadGroups" << endl;
		exit(-1);
	}
	libraries = libs;
	
	int Nset=makeSetsFromLibs(name,pars,anchors);
	if (Nset<1) {
		cerr << " problem getting libraries " << endl;
		exit(-1);
	}
	
	//----------------------------------------------------------------------------
	// read mapping histograms
	//----------------------------------------------------------------------------
	for (int iset=0; iset<Nset; iset++) {
		set[iset].initContigs();
	}
	
	//---------------------------------------------------------------------------
	// the number of reads in the archive
	//---------------------------------------------------------------------------
	//MaxFragments =900000000;
	
	cout << "-----------------------------------------------------------" <<   endl;
	cout << " load " << file1 <<   endl;
	cout << "-----------------------------------------------------------" <<   endl;
	
	// Bam structure
	BamAlignment ba1;
	
	string Name0="";
	//int Nuu[2]={0,0};
	
	// fragment counter
	int Nfrag = 0;
	
	//bool checkName=false;
	//---------------------------------------------------------------------------
	// loop through file...
	//---------------------------------------------------------------------------
	
	while ( ar1.GetNextAlignment(ba1) && (Nfrag<MaxFragments) ) {
       
		Nfrag++;
		set[0].Nfrag++;
		
		unsigned int p = ba1.Position;
		unsigned short a = ba1.RefID;
		string name = set[0].anchors.names[a];
		if (this->set[0].contig[name].Length>0) {
			set[0].contig[name].repeat.n[p]++;
		}
		
		//--------------------------------------------------------------------------
		// occasional status report 
		//--------------------------------------------------------------------------
		if ( (dbg!=0)&&(elapsedtime()>float(dbg))) {
			time(&tprev);
			cout << Nfrag << "\t" << ba1.RefID <<  "\t" << ba1.Position <<  "\n";
		}      		
		
    }
	
		
			
	// finalize histos
	for (int iset=0; iset<Nset; iset++) {
		
		C_librarymap::iterator imap;
		if (set[iset].libraries.libmap.size()!=1) {
			cerr << "==============================" << endl;
			cerr << " expected one library per set " << endl;
			cerr << set[iset].getSetName()  << " has " << set[iset].libraries.libmap.size() << endl;
			cerr << "==============================" << endl;
		}		
		
	}
	cout << "-----------------------------------------------------------------\n";
	cout << "final stats: "<< file1 << endl;
	summaryLog();
	cout << "-----------------------------------------------------------------\n";
	
}

//------------------------------------------------------------------------------
// make multiple sets for each read group in bam for scan   
//------------------------------------------------------------------------------
int C_pairedfiles::makeSetsFromLibs(string & name, RunControlParameters & pars, C_anchorinfo & anchors) {
  
  C_librarymap::iterator it;
  int N = libraries.libmap.size();
  
  ReadGroupID2set.clear();
	
  setNames.clear();
  set.clear();

  cout << " Number of Libraries:\t " <<  N<< endl;
  
  int iset=0;
  for ( it=libraries.libmap.begin() ; it != libraries.libmap.end(); it++ )
  {

    // cout << (*it).first << " => " << (*it).second << endl;
    //unsigned int ReadGroupCode=(*it).first;
    C_libraryinfo lib1 = (*it).second; 
    
    string setName = lib1.Info.ReadGroupID;
    C_set set1(setName, pars); 
    setNames.push_back(setName);
    set.push_back(set1);
    
    // hash from ReadGroupID to code
    ReadGroupID2Code[setName]=lib1.Info.ReadGroupCode;  

    // index from ReadGroupID to set 
    ReadGroupID2set[setName]=iset;

    // index from ReadGroupCode to set 
    ReadGroupCode2set[lib1.Info.ReadGroupCode]=iset;
    
    set[iset].anchors = anchors;
    set[iset].elementMinAnchor=anchors.anchorMinElement(); 
    set[iset].setFileName(name);

    // keep only lib in set libraries=libs; ????	  
	  
	  C_libraries libs1=libraries;
    libs1.libmap.clear();
    libs1.libmap[lib1.Info.ReadGroupCode]=lib1;

    set[iset].libraries = libs1; 
    // anchorinfo -> libraries
    set[iset].libraries.anchorinfo=anchors;

    iset++;

  }

  return iset;
}

//------------------------------------------------------------------------------
// make one set for all read group in bam for build  
//------------------------------------------------------------------------------
int C_pairedfiles::makeOneSetFromLibs(string & name,
									RunControlParameters & pars, C_anchorinfo & anchors) {
	
	C_librarymap::iterator it;
	int N = libraries.libmap.size();
	
	ReadGroupID2set.clear();
	
	C_HistoGroups H= pars.getHistoGroups();

	cout << " Number of Libraries:\t " <<  N<< endl;

	setNames.clear();
	set.clear();
	
	C_set set1(name, pars); 
	setNames.push_back(name);
	set.push_back(set1);
	
	int iset=0;
	
	set[iset].anchors = anchors;
	set[iset].elementMinAnchor=anchors.anchorMinElement(); 
	set[iset].setFileName(name);

	
	for ( it=libraries.libmap.begin() ; it != libraries.libmap.end(); it++ )
	{
		
		// cout << (*it).first << " => " << (*it).second << endl;
		//unsigned int ReadGroupCode=(*it).first;
		C_libraryinfo lib1 = (*it).second; 
				
		// hash from ReadGroupID to code
		ReadGroupID2Code[lib1.Info.ReadGroupID]=lib1.Info.ReadGroupCode;  
		
		// index from ReadGroupID to set 
		ReadGroupID2set[name]=iset;
		
		// index from ReadGroupCode to set 
		ReadGroupCode2set[lib1.Info.ReadGroupCode]=iset;
		
		
		// fragment length 
		int HGI = H.ReadGroupIndex[lib1.Info.ReadGroupID];		
		libraries.libmap[lib1.Info.ReadGroupCode].fragHist = H.Groups[HGI].h["LF"];

		// read length 
		libraries.libmap[lib1.Info.ReadGroupCode].readLengthHist = H.Groups[HGI].h["LR"];
		
			
	}
	
	// set threshold for each library
	double tailcut = pars.getFragmentLengthWindow();	
	libraries.resetFragLimits(tailcut);

	set[iset].libraries = libraries; 

	return 1;
}

/*
//------------------------------------------------------------------------------
// samtools compare function for query name sorting
//------------------------------------------------------------------------------
int C_pairedfiles::strnum_cmp(const string & a0, const string & b0)
{
  const char *a = a0.c_str();
  const char *b = b0.c_str();
	char *pa, *pb;
	pa = (char*)a; pb = (char*)b;
	while (*pa && *pb) {
		if (isdigit(*pa) && isdigit(*pb)) {
			long ai, bi;
			ai = strtol(pa, &pa, 10);
			bi = strtol(pb, &pb, 10);
			if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
		} else {
			if (*pa != *pb) break;
			++pa; ++pb;
		}
	}
	if (*pa == *pb)
		return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
	return *pa<*pb? -1 : *pa>*pb? 1 : 0;
}
*/




/*
//------------------------------------------------------------------------------
// extract read mapping info from MosaikAligner file - put into PairData objects
//------------------------------------------------------------------------------
void C_pairedfiles::loadMosaik( string  & file1, RunControlParameters & pars)
{
  // bam file pattern
  string patternFile("^.* /(.+)(\\.align|\\.dat|\\.mka)");
  // get allowed contig name regexp
  string AllowContigsRegex=pars.getAllowContigsRegex();
  // Read info ID
  string IDRegex=pars.getReadIDRegex();
  string patternID(IDRegex);
  // Read info Number
  string NumberRegex=pars.getReadNumberRegex();
  string patternNumber(NumberRegex);
  // regex match
  string match;
  // Debug level 
  int dbg =pars.getDbg();
  // flip relative to Illumina short PE FR convention
  MateMode =pars.getReadPairSenseConfig();
  // Scan only? 
  bool scan =pars.getDoScanOnly();
  // get Max Fragments
  int MaxFragments  = pars.getMaxFragments();

  // get Minimum Q mapping value for unique read
  Qmin  = pars.getQmin();

  // declare some relevant stringy vars
  string s,s1,s2;
  // read identifiers 
  string readName=" "; 
  // declare   pair read objects
  C_pairedread pair1, pair2;
  // set name
  string name="MSK";
  if  (RE2::FullMatch(file1.c_str(),patternFile.c_str(),&match) )  {
      name=match;
  }


  // ReadGroups
  unsigned int ReadGroupCode = 0;
  unsigned int ReadGroupCode0 = 0;

  // data set
  C_set set1(name, pars); 
  setNames.push_back(name);
  // local set class
  set.push_back(set1);
  //

  //---------------------------------------------------------------------------
  // extract set name from MOSAIK alignment filename
  //---------------------------------------------------------------------------
  if (RE2::FullMatch(file1.c_str(),patternFile.c_str(),&match) ) {
    string f = match;
    set[0].setFileName(f);
  }
  
  //---------------------------------------------------------------------------
  // Check if this is a valid MOSAIK alignment file (optional)
	//---------------------------------------------------------------------------
  if (!Mosaik::CAlignmentReader::CheckFile(file1, true))  {
    cout << "Error opening " << file1 << " for input" << endl;
    exit(-1);
  }
  
  
  //---------------------------------------------------------------------------
  // open the MOSAIK alignments file
  //---------------------------------------------------------------------------
  Mosaik::CAlignmentReader ar;
  ar.Open(file1);

  //---------------------------------------------------------------------------
  // the number of reads in the archive
  //---------------------------------------------------------------------------
  int Ntotal =int(ar.GetNumReads());
  // revise MaxFragments if warranted
  if (Ntotal<MaxFragments) { MaxFragments = Ntotal;} 
  cout << "-----------------------------------------------------------" <<   endl;
  cout << " load " << file1 <<   endl;
  cout << " Pairs " << Ntotal <<  "\t Max " << MaxFragments <<   endl;  

  //---------------------------------------------------------------------------
  // anchorinfo 
  //---------------------------------------------------------------------------
  C_anchorinfo anchor1(ar); 
  anchors=anchor1;
  // override anchors with anchorfile? 
  string anchorfile1=pars.getAnchorfile();
  if (anchorfile1.size()>0) {
    C_anchorinfo anchorf(anchorfile1);
    anchors=anchorf;
  }
  if (scan) {
     AllowContigsRegex = "doScanOnlyDoesNotUseContigs!";
  }
  anchors.anchorlimit(AllowContigsRegex);
  // the only set is set[0]
  set[0].anchors = anchors;
  set[0].elementMinAnchor=anchors.anchorMinElement(); 

  //---------------------------------------------------------------------------
	// fetch library info from MOSAIK alignments file
  //---------------------------------------------------------------------------
  C_libraries libs(ar);
  if (libs.libmap.size()>1) {
    cerr << "Multiple ReadGroups" << endl;
    cerr << libs << endl;
    exit(-1);
  }
  set[0].libraries = libs; 
  // anchorinfo -> libraries
  set[0].libraries.anchorinfo=anchors;
  // ReadGroupCode from header
  ReadGroupCode0=(*set[0].libraries.libmap.begin()).first;
  //----------------------------------------------------------------------------
  // read mapping histograms
  //----------------------------------------------------------------------------
  set[0].initContigs();
	
  //---------------------------------------------------------------------------
  // declare two Mosaik read objects
  //---------------------------------------------------------------------------

  Mosaik::AlignedRead mr;
  // debug
  int Nuu[2]={0,0};
  //---------------------------------------------------------------------------
  // loop through file...
  //---------------------------------------------------------------------------
  while(ar.LoadNextRead(mr)&&(int(set[0].Nfrag)<=MaxFragments)) {
    set[0].Nfrag++;
    // skip if read 
    if (mr.Name.size()==0) { 
        cout << " bad record " << set[0].Nfrag << " in " << file1 << endl;
        exit(-1);
    }      
    // unique read id part of read name
    readName="";
    if (RE2::FullMatch(mr.Name.c_str(),patternID.c_str(),&match) ) {
          readName=match;
    }
    pair1 = Mosaik2pair(mr);

    // ReadGroupCode
    ReadGroupCode = mr.ReadGroupCode;
    if ( (set[0].Nfrag==1) &(ReadGroupCode!=ReadGroupCode0) ) {  // hacked fix for sim
      set[0].libraries.libmap[ReadGroupCode]=set[0].libraries.libmap[ReadGroupCode0];
      set[0].libraries.libmap[ReadGroupCode].Info.ReadGroupCode=ReadGroupCode;
      set[0].libraries.libmap.erase(ReadGroupCode0);   
      ReadGroupCode0=ReadGroupCode;
    }
    if ( ReadGroupCode!=ReadGroupCode0) {
      cerr << "Changing ReadGroup\t " << ReadGroupCode0 << "\t to \t"<< ReadGroupCode << endl;
      exit(-1);
    }
    ReadGroupCode0=ReadGroupCode;
    
    int e;
    for (e = 0; e<2; e++) {
      int Nmap = pair1.read[e].align.size();
      set[0].repeatStats.Fill1(Nmap);   
      //--------------------------------------------------------------------------
      // repeat map
      //--------------------------------------------------------------------------
      if ( (Nmap>1) && (!scan) ) {
         set[0].processRepeat(pair1.read[e]); 
      }
    }
    
    //--------------------------------------------------------------------------
    // occasional status report 
    //--------------------------------------------------------------------------
    if ( (dbg!=0)&&(elapsedtime()>float(dbg))) {
      time(&tprev);
      if (scan) {
        printf("\t Reads %15d \t uufragments %12d\t normal %12d\n",
          set[0].Nfrag,Nuu[0]+Nuu[1],Nuu[1]);
      } else {  
        cout << set[0];
      }
    }      
    
    //------------------------------------------------------------------------
    // this is a complete pair
    //------------------------------------------------------------------------
    if (scan) {
      pair2=pair1;
    } else {
      // pair2 is the same as pair1 except the closest constrained pair is first
      pair2 = set[0].resolvePairConstraint(pair1);
      
      //--------------------------------------------------------------------------
      // debug
      //--------------------------------------------------------------------------
      int N0=pair2.read[0].align.size();
      int N1=pair2.read[1].align.size();      
      bool um=( ( (N0>1)&&(N1==1) ) || ( (N1>1)&&(N0==1) ) ); 
      if (um && (dbg >10000) && ((pair2.read[0].element==1)|| (pair2.read[1].element==1)) ) {
				int lmt=set[0].Fraglength(pair2);
				if (abs(lmt-140)>100) {
					cerr << " unconstrained alu " << mr.Name << "\t " << N0 
					<< "\t " << pair2.read[0].align[0].anchor 
					<< "\t " << pair2.read[0].align[0].pos+1 << "\t " << N1
					<< "\t " << pair2.read[1].align[0].anchor  
					<< "\t " << pair2.read[1].align[0].pos+1 << "\t " << lmt << endl;
				}
      }
      //--------------------------------------------------------------------------
 

    }
        
    bool bothends = (pair2.read[0].align.size()*pair2.read[1].align.size())>0;
    bool uu = (pair2.read[0].align.size()*pair2.read[1].align.size())==1;
    char constrain = 0; 
    if (bothends) {
        if (uu) {
          //constrain = (pair1.read[0].align.size()>1) + 2*(pair1.read[1].align.size()>1);          
          //----------------------------------------------------------------------
          // make some stats with this unique mapped fragment  
          //----------------------------------------------------------------------
          int lm = set[0].Fraglength(pair2);   
          if (lm!=-100000) {
             // make frag dist with high quality mapped pairs (q>30)
            if ( (pair1.read[0].align[0].q>30)&&(pair1.read[1].align[0].q>Qmin) ) { 
               set[0].fragStats.Fill1(lm);         
            }
            Nuu[1]++;  
          } else {
            Nuu[0]++;
          }
          set[0].lengthStats.Fill1(pair2.read[0].align[0].len);         
          set[0].lengthStats.Fill1(pair2.read[1].align[0].len);   
          set[0].qStats.Fill1(pair2.read[0].align[0].q);         
          set[0].qStats.Fill1(pair2.read[1].align[0].q);   
        }
        //------------------------------------------------------------------------
        // process pair for SV info
        //------------------------------------------------------------------------
        if (!scan) {
          set[0].processPair(pair2,constrain);
        }
        //------------------------------------------------------------------------
        // count paired read alignments for set[0] pairstats
        //------------------------------------------------------------------------
        int N1[2]={int(pair1.read[0].align.size()), int(pair1.read[1].align.size())};
        if (N1[0]*N1[1]==1) set[0].Npair_11o++; 
        int N2[2]={int(pair2.read[0].align.size()), int(pair2.read[1].align.size())};
        int n2[2]= {N2[0],N2[1]};
        if (n2[1]<n2[0]) { n2[1]=N2[0]; n2[0]=N2[1];}
        if (n2[0]>2) n2[0]=2;
        if (n2[1]>2) n2[1]=2;
        int npc = n2[0]*3 + n2[1];
        set[0].pairCountStats.Fill1(npc);   // e is paired            
    } else { 
        //------------------------------------------------------------------------
        // not pair - dangling end
        //------------------------------------------------------------------------
        int np1 = pair1.read[0].align.size();
        int np2 = pair1.read[1].align.size();
        int np = (np1==0? np2: np1);
        if (np>2) { np=2; }
        set[0].pairCountStats.Fill1(np);  
        if (!scan) {
          // dangling missing ends
            set[0].processPair(pair1,false);
        }      
    }
  }
  // finalize histos
  set[0].pairCountStats.Finalize();         
  set[0].qStats.Finalize();         
  set[0].fragStats.Finalize();         
  set[0].lengthStats.Finalize();
  set[0].repeatStats.Finalize();
  cout << "-----------------------------------------------------------------\n";
  cout << "final stats: "<< file1 << endl;
  summaryLog();
  cout << "-----------------------------------------------------------------\n";
  
  // finalize library: fragment length
  //double lowf = (1.0-pars.getFragmentLengthWindow()/100.0)/2.0;
  //double highf = (1.0-lowf); // half tail on high side  
  HistObj h=set[0].fragStats.h;
  // store fraglength info from pars to library (override this data set)
  set[0].libraries.libmap[ReadGroupCode].LM=int(pars.getFragmentLength()); //int(h.mode);
  set[0].libraries.libmap[ReadGroupCode].tailcut = pars.getFragmentLengthWindow();
  set[0].libraries.libmap[ReadGroupCode].LMlow=int(pars.getFragmentLengthLo()); ;
  set[0].libraries.libmap[ReadGroupCode].LMhigh=int(pars.getFragmentLengthHi()); //int(h.p2xTrim(highf));
  // store fraglength histogram from this data set to library
  set[0].libraries.libmap[ReadGroupCode].fragHist=h;

    
  // Read length
  h=set[0].lengthStats.h;
  set[0].libraries.libmap[ReadGroupCode].readLengthHist=h;
  set[0].libraries.libmap[ReadGroupCode].LR=int(set[0].lengthStats.h.mode);
  set[0].libraries.libmap[ReadGroupCode].LRmax=int(h.xc[h.xc.size()-1]);
  int nb=h.n.size();
  int xmin= set[0].libraries.libmap[ReadGroupCode].LRmax;
  for (int ib=nb; ib>0; ib--) {
     if ( (h.n[ib-1]>0)&&(h.xc[ib-1]<xmin) ) xmin=h.xc[ib-1];
  }
  set[0].libraries.libmap[ReadGroupCode].LRmin=xmin;

  // Redundancy
  
  set[0].calcLibraryRedundancy(set[0].libraries.libmap[ReadGroupCode]);
  
}

/*
//--------------------------------------------------------------------------
// utility function to extract Mosaik alignment into Spanner pairedread form 
//--------------------------------------------------------------------------
C_pairedread   C_pairedfiles::Mosaik2pair(  Mosaik::AlignedRead& mr) {

  C_pairedread pair1a;	
  
  unsigned int numAlignments1 = mr.Mate1Alignments.size();
  unsigned int numAlignments2 = mr.Mate2Alignments.size();
  
  if (numAlignments1>50000) {
    cerr << "Nalign " << numAlignments1 << " " << mr.Name << endl;
  }
  if (numAlignments2>50000) {
    cerr << "Nalign " << numAlignments2 << " " << mr.Name << endl;
  }
  
  C_readmaps read1;	
  bool mobhit=false;
	
  for(unsigned int i = 0; i < numAlignments1; i++) {
    C_readmap r1;       
    // 
    r1.q = mr.Mate1Alignments[i].Quality;
		
    // *** skip low quality alignments ***
    if ( (int(r1.q)< Qmin)&&(!scan) )   continue;
    // *** skip false quality> 100 alignments ***
    //if ( (int(r1.q)< Qmin) | (int(r1.q)>100)  ) continue;
		
    //--------------------------------------------------------------------------
    r1.pos=mr.Mate1Alignments[i].ReferenceBegin;
		// r1.len=1+mr.Mate1Alignments[i].QueryEnd-mr.Mate1Alignments[i].QueryBegin;
		r1.len=1+mr.Mate1Alignments[i].ReferenceEnd-mr.Mate1Alignments[i].ReferenceBegin;
		r1.anchor=mr.Mate1Alignments[i].ReferenceIndex;
    r1.sense=(mr.Mate1Alignments[i].IsReverseStrand ? 'R' : 'F');
    r1.mm = mr.Mate1Alignments[i].NumMismatches;
    //--------------------------------------------------------------------------
    // transform 454 long PE to Illumina short PE convention 
    //--------------------------------------------------------------------------
    if (MateMode==MATEMODE_454) r1.sense=(mr.Mate1Alignments[i].IsReverseStrand ? 'F' : 'R');
		
    if (r1.pos>=anchors.L[anchors.names[r1.anchor]]) {
			// cerr << " bad coordinate " << mr.Name << endl;
			continue;
    }
		
    
		// elements
    if (r1.anchor>=(anchors.element.size())) {
			cerr << " anchor mixup in Mosik2pair " << r1.anchor << " > " << anchors.element.size() << endl;
			cerr << " read " << 	mr.Name  << " 1" << endl;
			cerr << " Number of Alignments " << 	numAlignments1  << endl;
			cerr << " Alignment anchor pos len sense " << 	r1.anchor  << "\t " << r1.pos << "\t " << r1.len << "\t " << r1.sense << endl;       
			continue;
    }
    // check for mobhits
    //mobhit=mobhit|(r1.anchor>24);
    mobhit=mobhit|(anchors.element[r1.anchor]>0);
		
    //--------------------------------------------------------------------------
    // skip if mobhit has more than 10 bases in poly-A tail (L1HS p>6030) -kludge 9/16/2009
    //--------------------------------------------------------------------------
    if ( mobhit && ((r1.pos+r1.len)>6040) && (r1.anchor>56) )  {
      // cerr << "\t skip " << r1.anchor << "\t " << r1.pos << endl;
      continue;
    }
    
    read1.align.push_back(r1);
    // mobile element
    int me = anchors.element[r1.anchor];
    if (me>=0) {
      read1.element = read1.element | (1<<me);
    }
  }
  read1.Nalign=numAlignments1;
  pair1a.read[0]=read1;
  read1.align.clear();      
  read1.element=0; 
	
  for (unsigned int i = 0; i < numAlignments2; i++) {
    C_readmap r2;       
		
		r2.q = mr.Mate2Alignments[i].Quality; 
		// *** skip low quality alignments ***
		if ( (int(r2.q)< Qmin) && (!scan)   ) continue;
		// *** skip false quality> 100 alignments ***
		//if ( (int(r2.q)< Qmin) | (int(r2.q)>100)  ) continue;
		
    r2.pos=mr.Mate2Alignments[i].ReferenceBegin; 
    //r2.len=1+mr.Mate2Alignments[i].QueryEnd-mr.Mate2Alignments[i].QueryBegin;
		r2.len=1+mr.Mate2Alignments[i].ReferenceEnd-mr.Mate2Alignments[i].ReferenceBegin;
    r2.anchor=mr.Mate2Alignments[i].ReferenceIndex;
    r2.sense=(mr.Mate2Alignments[i].IsReverseStrand ? 'R' : 'F');
    r2.mm = mr.Mate2Alignments[i].NumMismatches;
    //--------------------------------------------------------------------------
    // transform AB SOLiD long PE to Illumina short PE convention 
    //--------------------------------------------------------------------------
    if (MateMode==MATEMODE_SOLID) r2.sense=(mr.Mate2Alignments[i].IsReverseStrand ? 'F' : 'R');
		
    if (r2.pos>=anchors.L[anchors.names[r2.anchor]]) {
			// cerr << " bad coordinate " << mr.Name << endl;
			continue;
    }
    
    if (r2.anchor>=(anchors.element.size())) {
			cerr << " anchor mixup in Mosik2pair " << r2.anchor << " > " << anchors.element.size() << endl;
			cerr << " read " << 	mr.Name  << " 2" << endl;
			cerr << " Number of Alignments " << 	numAlignments2  << endl;
			cerr << " Alignment anchor pos len sense " << 	r2.anchor  << "\t " << r2.pos << "\t " << r2.len << "\t " << r2.sense << endl;       
			continue;
    }
		
    //mobhit=mobhit|(r2.anchor>24);
    mobhit=mobhit|(anchors.element[r2.anchor]>0);
		
    //--------------------------------------------------------------------------
    // skip if mobhit is in poly-A tail (L1HS p>6030) -kludge 9/14/2009
    //--------------------------------------------------------------------------
    if ( mobhit && (r2.pos>6025) && (r2.anchor>56) )  {
      // cerr << "\t skip " << r2.anchor << "\t " << r2.pos << endl;
      continue;
    }
    
    read1.align.push_back(r2);
    // mobile element
    int me = anchors.element[r2.anchor];
    if (me>=0) {
      read1.element = read1.element | (1<<me);
    }
		
  }
  
  read1.Nalign=numAlignments2;
  pair1a.read[1]=read1;
  pair1a.ReadGroupCode=mr.ReadGroupCode;
  return pair1a;
}
*/

//------------------------------------------------------------------------------
// convert Bam read record with ZA tag to Mosaik aligned pair record 
// returns:
//         -1 = error 
//          0 = no reads
//          1 = one read
//          2 = both ends
//------------------------------------------------------------------------------
int  C_pairedfiles::BamZA2PairedRead(BamAlignment & ba1, C_pairedread & pr1) 
{      
	
	C_readmap r1,r2;
	C_readmaps rr1,rr2;
	
	// ZA tag pattern
	RE2 patternZA("<(.);(\\d*?);(\\d*?);(.*?);(\\d*?);(.*?);(.*?)>");
	// <@;Q1;Q2;Mob;# mappings;;> <=;Q1;Q2;Mob;# mappings;CIGAR;MD>
	// example: <@;41;1;L1;226;;><=;35;0;;1;100M;28T5T27G37>
	
	int NZA=0;
	string tagZA = "ZAs";
	string tagNM = "NMs";
	string tagMD = "MDs";
	string ZA,mob;
	if (!ba1.GetTag(tagZA,ZA) ) { 
		return -1;
	} 			  
	int q1,q2,nmap1;
	string this1,mob1,cig1,md1;
	string ReadGroupID1;
	
	StringPiece ZASP(ZA);    // Wrap a StringPiece around it
	int lenQ,lenR, mm;
	
	while (RE2::FindAndConsume(&ZASP,patternZA,&this1,&q1,&q2,&mob1,&nmap1,&cig1,&md1) ) {
		
		char t1 = this1[0];
		
		switch (t1) {
				
			case '@':
				
				NZA++;
							
				// 
				r1.anchor=ba1.RefID;
				// define start of query at start of mapped part
				// length of read mapping in query coordinates
				r1.len=BamCigarData2Len(ba1.CigarData,1);				
				r1.pos=ba1.Position;
				r1.sense=(ba1.IsReverseStrand()? 'R': 'F');
				r1.q=q1;
				r1.q2=q2;
				r1.nmap=nmap1;
				r1.mob=mob1;
				
				
				// mismatches NM
				if (ba1.GetTag(tagNM,mm)) {
					r1.mm=(mm>=0? mm: 0);
				} else { 
					int mm=BamCigarData2mm(ba1.CigarData);
					string MD;
					if (ba1.GetTag(tagMD,MD)) {						
						mm+=getMDMismatchCount(MD);
					}
					r1.mm=mm;
				}
				rr1.align.push_back(r1);
				rr1.Nalign=r1.nmap;
				
				// mark reads hitting elements
				rr1.element=r1.mob.size();				
				
				break;
				
			// mate
			case '&':
			case '=':  
				
				NZA++;
				
				if (!getCigarLengths(cig1,lenQ,lenR,mm)) continue;
				
				r2.anchor=ba1.MateRefID;
				r2.len=lenR;
				r2.pos=ba1.MatePosition;
				r2.sense=(ba1.IsMateReverseStrand()? 'R': 'F');
				r2.q=q1;
				r2.q2=q2;
				r2.nmap=nmap1;
				r2.mob=mob1;
				
				// fix this with MD
				mm+=getMDMismatchCount(md1);
				r2.mm=mm;
				
				rr2.align.push_back(r2);
				rr2.Nalign=r2.nmap;
				
				// mark reads hitting elements
				rr2.element=r2.mob.size();				
				
				
		}
	}

	// convert ReadGroupID to readGroupCode (Mosaik/Spanner) 				
	ba1.GetReadGroup(ReadGroupID1);							
	// all libraries stored in every set for this map to work 
	pr1.ReadGroupCode=ReadGroupID2Code[ReadGroupID1];
	
	if (NZA>1) { 
		// 454 & SOLiD pair orientation gymnastics
		if (MateMode==MATEMODE_454) { // 454
			if (ba1.IsFirstMate()) {
				rr1.align[0].sense=!rr1.align[0].sense;
			} else {
				rr2.align[0].sense=!rr2.align[0].sense;
			}
		} else if (MateMode==MATEMODE_SOLID) { // SOLiD
			if (!ba1.IsFirstMate()) {
				rr1.align[0].sense=!rr1.align[0].sense;
			} else {
				rr2.align[0].sense=!rr2.align[0].sense;
			}
		}
  }
	
	// convert ReadGroupID to readGroupCode (Mosaik/Spanner) 				
	ba1.GetReadGroup(ReadGroupID1);							
	// all libraries stored in every set for this map to work 
	pr1.ReadGroupCode=ReadGroupID2Code[ReadGroupID1];
	
	pr1.read[0]=rr1;
	pr1.read[1]=rr2;
	

	if (NZA < 1)  {
		cerr << "no ZA tag info " << endl;
	}
	 
	return(NZA);
}

//------------------------------------------------------------------------------
// convert Bam read record with ZA tag to Mosaik aligned pair record 
// returns:
//         -1 = error 
//          0 = no reads
//          1 = one read
//          2 = both ends
//------------------------------------------------------------------------------
int  C_pairedfiles::BamBam2PairedRead(BamAlignment & ba1, BamAlignment & ba2, C_pairedread & pr1)
{      
	
	C_readmap r1,r2;
	C_readmaps rr1,rr2;
		
	string tagNM = "NMs";
	string tagMD = "MDs";
	string ReadGroupID1;
	
	
	
	int Nmap=1,mm=0;
				
	// 
	r1.anchor=ba1.RefID;
	r1.len=BamCigarData2Len(ba1.CigarData,0);				
	r1.pos=ba1.Position;
	r1.sense=(ba1.IsReverseStrand()? 'R': 'F');
  r1.q=ba1.MapQuality;
	r1.q2=0;
	r1.nmap=(r1.q>0? 1: 2);
				
	// mismatches NM
	if (ba1.GetTag(tagNM,mm)) {
		r1.mm=(mm>=0? mm: 0);
	} else { 
		int mm=BamCigarData2mm(ba1.CigarData);
		string MD;
		if (ba1.GetTag(tagMD,MD)) {						
			mm+=getMDMismatchCount(MD);
		}
		r1.mm=mm;
	}

	rr1.align.push_back(r1);
	rr1.Nalign=r1.nmap;
	
	if (ba2.IsMapped()) { 
		
		Nmap=2;
		r2.anchor=ba2.RefID;
		r2.len=BamCigarData2Len(ba2.CigarData,0);		
		r2.pos=ba2.Position;
		r2.sense=(ba2.IsReverseStrand()? 'R': 'F');
		r2.q=ba2.MapQuality;
		r1.q2=0;
		r2.nmap=(r2.q>0? 1: 2);
		
		// mismatches NM
		if (ba2.GetTag(tagNM,mm)) {
			r2.mm=(mm>=0? mm: 0);
		} else { 
			int mm=BamCigarData2mm(ba1.CigarData);
			string MD;
			if (ba2.GetTag(tagMD,MD)) {						
				mm+=getMDMismatchCount(MD);
			}
			r2.mm=mm;
		}
		rr2.align.push_back(r2);
		rr2.Nalign=r2.nmap;
		
	}

	
	if (Nmap>1) { 
		// 454 & SOLiD pair orientation gymnastics
		if (MateMode==MATEMODE_454) { // 454
			if (ba1.IsFirstMate()) {
				rr1.align[0].sense=!rr1.align[0].sense;
			} else {
				rr2.align[0].sense=!rr2.align[0].sense;
			}
		} else if (MateMode==MATEMODE_SOLID) { // SOLiD
			if (!ba1.IsFirstMate()) {
				rr1.align[0].sense=!rr1.align[0].sense;
			} else {
				rr2.align[0].sense=!rr2.align[0].sense;
			}
		}
	}
	
	// convert ReadGroupID to readGroupCode (Mosaik/Spanner) 				
	ba1.GetReadGroup(ReadGroupID1);							
	// all libraries stored in every set for this map to work 
	pr1.ReadGroupCode=ReadGroupID2Code[ReadGroupID1];
	
	pr1.read[0]=rr1;
	pr1.read[1]=rr2;
		
	return Nmap;
}
		
		


//------------------------------------------------------------------------------
// convert Bam read record with ZA tag to Mosaik aligned pair record 
// returns:
//         -1 = error 
//          0 = no reads
//          1 = one read
//          2 = both ends
//------------------------------------------------------------------------------
int  C_pairedfiles::BamSpecial2PairedRead(BamAlignment & ba1, BamAlignment & ba2, C_pairedread & pr1)
{      
	
	C_readmap r1,r2;
	C_readmaps rm1,rm2,ra1,ra2;
	C_pairedread prm,pra;
	
	string tagNM = "NMs";
	string tagMD = "MDs";
	string tagZA = "ZAs";
	int NZA=0;
	string ZAa,ZAm;
	string ReadGroupID1;
	int q1,q2,nmap1,me1;
	string this1,mob1,cig1,md1;
	int lenQ,lenR, mm;
	
	
	// ZA tag pattern
	RE2 patternZA("<(.);(\\d*?);(\\d*?);(.*?);(\\d*?);(.*?);(.*?)>");
	if (!ba1.GetTag(tagZA,ZAa) ) { 
		return -1;
	} 			  
	if (!ba2.GetTag(tagZA,ZAm) ) { 
		return -1;
	} 			  

	StringPiece ZASPm(ZAm);    // Wrap a StringPiece around the moblist record ZAm 
	
	while (RE2::FindAndConsume(&ZASPm,patternZA,&this1,&q1,&q2,&mob1,&nmap1,&cig1,&md1) ) {
		
		char t1 = this1[0];
		
		switch (t1) {
				
			case '@':
				
				NZA++;
				
				r1.anchor=ba2.RefID;
				// define start of query at start of mapped part
				// length of read mapping in query coordinates
				r1.len=BamCigarData2Len(ba2.CigarData,1);				
				r1.pos=ba2.Position;
				r1.sense=(ba2.IsReverseStrand()? 'R': 'F');
				r1.q=q1;
				r1.q2=q2;
				r1.nmap=nmap1;
				r1.mob=mob1;
				
				
				// mismatches NM
				if (ba1.GetTag(tagNM,mm)) {
					r1.mm=(mm>=0? mm: 0);
				} else { 
					int mm=BamCigarData2mm(ba2.CigarData);
					string MD;
					if (ba1.GetTag(tagMD,MD)) {						
						mm+=getMDMismatchCount(MD);
					}
					r1.mm=mm;
				}
				rm1.align.push_back(r1);
				rm1.Nalign=r1.nmap;
				
								
				// mobile element
				me1 = anchors.element[r1.anchor];
				if (me1>=0) {
					rm1.element = rm1.element | (1<<me1);
				}
				
				break;
				
				// mate
			case '&':  
			case '=':  
				
				NZA++;
				
				if (!getCigarLengths(cig1,lenQ,lenR,mm)) continue;
				
				r2.anchor=ba2.MateRefID;
				r2.len=lenR;
				r2.pos=ba2.MatePosition;
				r2.sense=(ba2.IsMateReverseStrand()? 'R': 'F');
				r2.q=q1;
				r2.q2=q2;
				r2.nmap=nmap1;
				r2.mob=mob1;
				
				// fix this with MD
				mm+=getMDMismatchCount(md1);
				r2.mm=mm;
				
				rm2.align.push_back(r2);
				rm2.Nalign=r2.nmap;
				
				
				me1 = anchors.element[r2.anchor];
				if (me1>=0) {
					rm2.element = rm2.element | (1<<me1);
				}												
		}
	}

	// flip strands for non-illumina modes
	if (NZA>1) { 
		// 454 & SOLiD pair orientation gymnastics
		if (MateMode==MATEMODE_454) { // 454
			if (ba2.IsFirstMate()) {
				rm1.align[0].sense=!rm1.align[0].sense;
			} else {
				rm2.align[0].sense=!rm2.align[0].sense;
			}
		} else if (MateMode==MATEMODE_SOLID) { // SOLiD
			if (!ba2.IsFirstMate()) {
				rm1.align[0].sense=!rm1.align[0].sense;
			} else {
				rm2.align[0].sense=!rm2.align[0].sense;
			}
		}
	}

	NZA=0;
	
	prm.read[0]=rm1;
	prm.read[1]=rm2;
	
	
	StringPiece ZASPa(ZAa);    // Wrap a StringPiece around the official reference part ZAa
			
	while (RE2::FindAndConsume(&ZASPa,patternZA,&this1,&q1,&q2,&mob1,&nmap1,&cig1,&md1) ) {
		
		char t1 = this1[0];
		
		switch (t1) {
				
			case '@':
				
				NZA++;
				
				r1.anchor=ba1.RefID;
				// define start of query at start of mapped part
				// length of read mapping in query coordinates
				r1.len=BamCigarData2Len(ba1.CigarData,1);				
				r1.pos=ba1.Position;
				r1.sense=(ba1.IsReverseStrand()? 'R': 'F');
				r1.q=q1;
				r1.q2=q2;
				r1.nmap=nmap1;
				r1.mob=mob1;
				
				
				// mismatches NM
				if (ba1.GetTag(tagNM,mm)) {
					r1.mm=(mm>=0? mm: 0);
				} else { 
					int mm=BamCigarData2mm(ba1.CigarData);
					string MD;
					if (ba1.GetTag(tagMD,MD)) {						
						mm+=getMDMismatchCount(MD);
					}
					r1.mm=mm;
				}

				
				// here - this is the mate that mapped to a moblist
				// check if properpair.. if so then make this alignment precede 
				// the moblist alignment
				
				ra1.align.push_back(r1);
				ra1.Nalign=r1.nmap;
				
				break;
				
			// mate
			case '&':  
			case '=':  
				
				NZA++;
				
				if (!getCigarLengths(cig1,lenQ,lenR,mm)) continue;
				
				r2.anchor=ba1.MateRefID;
				r2.len=lenR;
				r2.pos=ba1.MatePosition;
				r2.sense=(ba1.IsMateReverseStrand()? 'R': 'F');
				r2.q=q1;
				r2.q2=q2;
				r2.nmap=nmap1;
				r2.mob=mob1;
				
				// fix this with MD
				mm+=getMDMismatchCount(md1);
				r2.mm=mm;
				
				ra2.align.push_back(r2);
				ra2.Nalign=r2.nmap;
												
		}
	}
	
	// flip strands for non-illumina modes
	if (NZA>1) { 
		// 454 & SOLiD pair orientation gymnastics
		if (MateMode==MATEMODE_454) { // 454
			if (ba1.IsFirstMate()) {
				ra1.align[0].sense=!ra1.align[0].sense;
			} else {
				ra2.align[0].sense=!ra2.align[0].sense;
			}
		} else if (MateMode==MATEMODE_SOLID) { // SOLiD
			if (!ba1.IsFirstMate()) {
				ra1.align[0].sense=!ra1.align[0].sense;
			} else {
				ra2.align[0].sense=!ra2.align[0].sense;
			}
		}
	}
	

	
	pra.read[0]=ra1;
	pra.read[1]=ra2;
	
	int LF = set[0].Fraglength(pra);
	string rgid;
	ba1.GetReadGroup(rgid);
	unsigned int rgcode=libraries.ReadGroupID2Code[rgid];
	int LMlow = libraries.libmap[rgcode].LMlow;
	int LMhigh = libraries.libmap[rgcode].LMhigh;
	bool properOrientation=false;
	bool rev = ra1.align[0].sense=='R';
	bool revMate = ra2.align[0].sense=='R';
	if (ra1.align[0].anchor==ra2.align[0].anchor) {
		if (ra1.align[0].pos<ra2.align[0].pos) {				
			properOrientation = (!rev)&&revMate;
		}	else {
			properOrientation = rev&&(!revMate);
		}
	}
	
	bool properpair=false;
	if ( (LF>LMlow)&&(LF<LMhigh)&&properOrientation) {					
		properpair=true;
	}	
	
	// if 
	if (properpair) {
		pr1=pra;
	} else {
		pr1=prm;
	}
	pr1.ReadGroupCode=rgcode;

	
	return(NZA);
	
}





		

//------------------------------------------------------------------------------
// parse Bam tagData for tag:i:N
//------------------------------------------------------------------------------
int  C_pairedfiles::parseBamTagData(string & tagData, string & tag) {
      int N=-1;
      string str = tagData;
      size_t f1 = str.find(tag); //"H0");
      if (f1!=string::npos)
      {
         f1=f1+tag.length();
         char nm=str[f1];
         N=nm;
      }
       /* 
          size_t f2=str.find_first_not_of("0123456789",f1);
          if (f2!=string::npos)  {
            str = str.substr(f1,f2-f1);
            N= string2Int (str);
          }
       */
      return N;
}

//------------------------------------------------------------------------------
// parse Bam tagData for tag:i:Z
//------------------------------------------------------------------------------
string  C_pairedfiles::parseBamTagDataString(string & tagData, string & tag) {
	string S="";
	string str = tagData;
	size_t f1 = str.find(tag); //"H0");
	if (f1!=string::npos)
	{
		
		vector<string> tv;
		int nt = split(tv,tag,":");
		if (nt>2) {
			S=tv[2];
		}
	}
	return S;
}

//------------------------------------------------------------------------------
// Bam file sort order
//------------------------------------------------------------------------------
string  C_pairedfiles::parseBamHeader(string & samheader, const string & tag) const {
      string str = samheader;
      string order="    ";
      size_t f1 = str.find(tag); //"SO");
      if (f1!=string::npos)
      {
         f1=f1+tag.length();
         order = str.substr (f1,4);;
      }
      return order;
}

  
//------------------------------------------------------------------------------
// parse Bam read name into end (0 or 1) and ID (name)
//------------------------------------------------------------------------------
bool  C_pairedfiles::parseBamReadName(string & name0, int & e1, string & id) {
      
      size_t Lname0=name0.size();
      string pe=name0.substr(Lname0-2,2); 
      e1=0;
      
      // sometime read names end in .1 or .2 indicating pair order
      
      if (pe==".1") {  
         e1=0;                   // end 0
      } else if (pe==".2") {
         e1=1;                   // end 1
      } else {
         e1=0;                   // single end
         if (pe[0]=='.') {
            return false;
         }        
      }     
      // ID0
      if (pe[0]=='.') {  
        id=name0.substr(0,Lname0-2); 
      } else {
        id=name0;
      }
      return true;
}



//------------------------------------------------------------------------------
// convert a Bam read pair to one complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentPair( BamReader & ar1,BamReader & ar2, C_pairedread & pr1) 
{
	// reset pr1
	pr1.read[0].align.clear();
	pr1.read[1].align.clear();
	pr1.read[0].element=0;
	pr1.read[1].element=0;
	pr1.read[0].Nalign=0;
	pr1.read[1].Nalign=0;
	pr1.Name="";
	pr1.ReadGroupCode=0;

	if (BamZ==1) {  
		// use ZA tag to extract  all pair info from first end in bam 
		return (nextBamAlignmentZA( ar1, pr1));
	}
	
	if (BamZ==2) {  
		// scan all pairs from first end in bam 
		return (nextBamAlignmentPairSortedByName( ar1, pr1));
	}
	
	if (BamZ==3) {  
		// scan all pairs from first end in bam 
		return (nextBamAlignmentPairSpecial( ar1, pr1));
	}
	
	// scan all pairs jumping to mates in bam 
	return (nextBamAlignmentJump( ar1, ar2, pr1));
	
}


//------------------------------------------------------------------------------
// convert a Bam read pair to one semi-complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentZA( BamReader & ar1, C_pairedread & pr) 
{
	
	// Mosaik structures
	//C_pairedread pr1;
	
	// Bam structure
	BamAlignment ba1;
	
	
	int ndone = 0;
	bool ok=false;
	bool scan = SpannerMode==SPANNER_SCAN;
	// bool build = SpannerMode==SPANNER_BUILD;
	
	bool rev,revMate,skip;
	bool properOrientation = false;
	
	while ( ar1.GetNextAlignment(ba1) ) {
		
		if (doneFrag[ba1.Name]) {
			continue;
		}
				
		doneFrag[ba1.Name]=true;
		
		// skip unmapped reads (another program somewhere ...)
		if (ba1.RefID<0) {
			continue;
		}	
		if (!ba1.IsMapped() ) {
			continue;
		}	
		
		int LF = ba1.InsertSize;
		// SLX RP only
		if (ba1.IsReverseStrand()) LF=-LF; 
		
		// check if fragment done already
		ndone++;
		if (ba1.Name.size()<1) {			
			cerr << "empty read name: " << ba1.RefID+1 <<"\t"<< ba1.Position+1<<"\t"<< ba1.MateRefID+1<<"\t"<< ba1.MatePosition+1 <<"\t"<< ndone << endl;
			continue;
		}
		
		//-------------------------------------------------------------------------------
		// skip proper pairs right here to save span file space from  build
		//-------------------------------------------------------------------------------
		skip=false;
  	if (!scan) { 
			
			string rgid;
			ba1.GetReadGroup(rgid);
			unsigned int rgcode=libraries.ReadGroupID2Code[rgid];
			int LMlow = libraries.libmap[rgcode].LMlow;
			int LMhigh = libraries.libmap[rgcode].LMhigh;
			
			// orientation transform to Illumina FR			
			rev=ba1.IsReverseStrand();
			revMate=ba1.IsMateReverseStrand();
			
			//454 -> IL
			if ( MateMode==MATEMODE_454) {
				if (ba1.IsFirstMate()) {
					rev=!rev;			
				} else {
					revMate=!revMate;
				}
				// SOLiD -> IL
			} else if ( MateMode==MATEMODE_SOLID) {
				if (ba1.IsFirstMate()) {
					revMate=!revMate;
				} else {
					rev=!rev;			
				}
			}		
			
			
			
			if (ba1.RefID==ba1.MateRefID) {
				if (ba1.Position<ba1.MatePosition) {				
					properOrientation = (!rev)&&revMate;
				}	else {
					properOrientation = rev&&(!revMate);
				}
			}
			
			if (ba1.IsMateMapped()&&(LF>LMlow)&&(LF<LMhigh)&&properOrientation) {					
				skip=true;
			}	
			
		}
		
		if (skip) continue;		
		
		ok=true;
		break;
		
	}
	
	// check out when GetNextAlignment returns false
	if (!ok) { 
		cerr << ba1.Name.c_str() << ba1.RefID << endl;
		return(ok);
	}	
	
	// mark this one as done
	doneFrag[ba1.Name]=true;
	
	int pe = BamZA2PairedRead(ba1, pr);
	
	
	return (pe>=0);		
}


//------------------------------------------------------------------------------
// convert a Bam read pair to one complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentPairSortedByName( BamReader & ar1, C_pairedread & pr) 
{
	
	
	//------------------------------------------------------------------------------
	// select abberant pairs that meet Qmin criteria, include dangling (unmapped) ends
	//------------------------------------------------------------------------------
	
	// pair structures
	//C_pairedread pr1;	

	// Bam structure
	BamAlignment ba1,ba2;
	
	bool scan = SpannerMode==SPANNER_SCAN;
	//bool build = SpannerMode==SPANNER_BUILD;
		
	int LF=0;
	bool findmate=false;	
	int NmateFound=-1;
	string FR="FR";
	bool rev,revMate;
	bool properOrientation = false;	
	
	
	while ( ar1.GetNextAlignment(ba1) ) {
		
		if (doneFrag[ba1.Name]) {
			continue;
		}
		
		doneFrag[ba1.Name]=true;
		
		// skip unmapped reads (another function in a module somewhere ...)
		if (ba1.RefID<0) {
			continue;
		}	
		
		// single ends not processed here
		if (!ba1.IsPaired()) {
  		continue;
		}
		
		
		// dangling end is processed
		if (ba1.IsMateMapped() ) {
			ar1.GetNextAlignment(ba2);
		} else {
			NmateFound=1;  // set to mark dangling end
			break;
		}	
		
		if ( (!scan) && (ba1.MapQuality<Qmin) ) {
			continue;
		}	
		
		
		//if (ba1.Position>=1130515) { 
		// cerr << endl; }
		
		// skip proper pairs	=  abberant pair criteria 	
		LF = ba1.InsertSize;
		// SLX RP only  ??? 
		if (ba1.IsReverseStrand()) LF=-LF; 
		string rgid;
		ba1.GetReadGroup(rgid);
		unsigned int rgcode=libraries.ReadGroupID2Code[rgid];
		int LMlow = libraries.libmap[rgcode].LMlow;
		int LMhigh = libraries.libmap[rgcode].LMhigh;
		
		// orientation transform to Illumina FR
    
		rev=ba1.IsReverseStrand();
		revMate=ba1.IsMateReverseStrand();
		
		//454 -> IL
		//454 -> IL
		if ( MateMode==MATEMODE_454) {
			if (ba1.IsFirstMate()) {
				rev=!rev;			
			} else {
				revMate=!revMate;
			}
			// SOLiD -> IL
		} else if ( MateMode==MATEMODE_SOLID) {
			if (ba1.IsFirstMate()) {
				revMate=!revMate;
			} else {
				rev=!rev;			
			}
		}		
		
		// now in Illumina convention... check proper pair orientation
		if (ba1.RefID==ba1.MateRefID) {
			if (ba1.Position<ba1.MatePosition) {				
				properOrientation = (!ba1.IsReverseStrand())&&ba1.IsMateReverseStrand();
			}	else {
				properOrientation = ba1.IsReverseStrand()&&(!ba1.IsMateReverseStrand());
			}
		}
		
		//-------------------------------------------------------------------------------
		// skip proper pairs
		//-------------------------------------------------------------------------------
		if ( (LF>LMlow)&&(LF<LMhigh)&&properOrientation&&(!scan)) {					
			continue;
		}	
		
		// artifact check for mate mapped = this read (illumina problem)
		if  ( (ba1.RefID==ba1.MateRefID) && (ba1.Position==ba1.MatePosition) ) { 
			continue;
		}
		
		//------------------------------------------------------------------------------
		// go mate hunting
		//------------------------------------------------------------------------------		
		findmate=false;		
		
		if (ba2.RefID!=ba1.MateRefID) {
			// this set of bam files doesn't have mateRefID.... log this somewhere other than cerr, cout
			cerr << " Mates out of order :   " << ba1.Name <<"\t" << ba1.MateRefID << "\t" << ba2.RefID << "\n";  // oops
			break;
		}
		
		// check read name should be the same but not the same read (depend on First/Second mate flags). 
		if ( ( ba1.Name.compare(ba2.Name)==0 ) && (ba1.IsFirstMate()!=ba2.IsFirstMate()) )  { 
			findmate=true;			
			NmateFound++;
			if (scan) break;
			if (ba2.MapQuality>=Qmin )  break;
			// low quality reads 
			NmateFound=0;
		}										
		
	}	
	
  if(NmateFound<0) {
		return(false);
	}
	
	NmateFound = BamBam2PairedRead(ba1, ba2, pr);
		
	return (NmateFound>0);		
}



//------------------------------------------------------------------------------
// convert a Bam read pair to one complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentJump( BamReader & ar1, BamReader & ar2, C_pairedread & pr) 
{
	
	//------------------------------------------------------------------------------
	// select abberant pairs that meet Qmin criteria, include dangling (unmapped) ends
	//------------------------------------------------------------------------------
		
	// Bam structure
	BamAlignment ba1,ba2;
	
	bool scan = SpannerMode==SPANNER_SCAN;
	//bool build = SpannerMode==SPANNER_BUILD;

	int LF=0;
	bool findmate=false;	
	int NmateFound=-1;
	string FR="FR";
	bool rev,revMate;
	bool properOrientation = false;	
	
	while ( ar1.GetNextAlignment(ba1) ) {
		
		if (doneFrag[ba1.Name]) {
			continue;
		}
		
		doneFrag[ba1.Name]=true;
		
		// skip unmapped reads (another function in a module somewhere ...)
		if (ba1.RefID<0) {
			continue;
		}	
		
		// single ends not processed here
		if (!ba1.IsPaired()) { 
			continue;
		}
		
		if ( (!scan)&&(ba1.MapQuality<Qmin) ) {
			continue;
		}	
		
		// dangling end is processed
		if (!ba1.IsMateMapped() ) {
			NmateFound=1;  // set to mark dangling end
			break;
		}	
				
		// skip proper pairs	=  abberant pair criteria 	
		LF = ba1.InsertSize;
		// SLX RP only  ??? 
		if (ba1.IsReverseStrand()) LF=-LF; 
		string rgid;
		ba1.GetReadGroup(rgid);
		unsigned int rgcode=libraries.ReadGroupID2Code[rgid];
		int LMlow = libraries.libmap[rgcode].LMlow;
		int LMhigh = libraries.libmap[rgcode].LMhigh;
		
		// orientation transform to Illumina FR
    
		rev=ba1.IsReverseStrand();
		revMate=ba1.IsMateReverseStrand();
		
		//454 -> IL
		if ( MateMode==MATEMODE_454) {
			if (ba1.IsFirstMate()) {
				rev=!rev;			
			} else {
				revMate=!revMate;
			}
			// SOLiD -> IL
		} else if ( MateMode==MATEMODE_SOLID) {
			if (ba1.IsFirstMate()) {
				revMate=!revMate;
			} else {
				rev=!rev;			
			}
		}		
		
		
		if (ba1.RefID==ba1.MateRefID) {
			if (ba1.Position<ba1.MatePosition) {				
				properOrientation = (!ba1.IsReverseStrand())&&ba1.IsMateReverseStrand();
			}	else {
				properOrientation = ba1.IsReverseStrand()&&(!ba1.IsMateReverseStrand());
			}
		}
		
		//-------------------------------------------------------------------------------
		// skip proper pairs
		//-------------------------------------------------------------------------------
		if (ba1.IsMateMapped()&&(LF>LMlow)&&(LF<LMhigh)&&properOrientation&&(!scan)) {					
			continue;
		}	
		
		// artifact check for mate mapped = this read (illumina problem)
		if  ( (ba1.RefID==ba1.MateRefID) && (ba1.Position==ba1.MatePosition) ) { 
			continue;
		}
		
		//------------------------------------------------------------------------------
		// go mate hunting
		//------------------------------------------------------------------------------		
		findmate=false;
		
		int matePosGuess = ba1.Position+ba1.InsertSize-1;
		
		if (ba1.RefID!=ba1.MateRefID)			{
			matePosGuess = ba1.MatePosition+5;
		}
		
		// jump to first alignment that overlaps mate
		if ( ar2.Jump(ba1.MateRefID, matePosGuess) ) {
			
			int nshot=0;
			
			while ( ar2.GetNextAlignment(ba2) ) {
				
				nshot++;
				
				if (ba2.RefID!=ba1.MateRefID) {
					Shots2Mate.push_back(nshot);
					// this set of bam files doesn't have mateRefID.... log this somewhere other than cerr, cout
					//cerr << "no mate RefID in bam:  want " << ba1.MateRefID << "\tget" << ba2.RefID << "\n";  // oops
					break;
				}
				
				// check read name should be the same but not the same read (depend on First/Second mate flags). 
				if ( ( ba1.Name.compare(ba2.Name)==0 ) && (ba1.IsFirstMate()!=ba2.IsFirstMate()) )  { 
					findmate=true;					
					NmateFound=2;
					Shots2Mate.push_back(nshot);
					//if (nskip>500)  { 
					//	cerr << "bam find after large skip " << nskip << "\n";  // oops
					//	cerr << ba1.Name << "\t " << matePosGuess << "\t " << ba2.Position << "\t " << ba2.Length << endl;
					//}
					break;
				}
				
				if (ba2.Position>ba1.MatePosition) {
					cerr << "bam jump too far " << ba1.MatePosition << "\t" << ba2.Position << "\n";  // oops
					// jump to first alignment that overlaps mate
					if ( ar2.Jump(ba1.MateRefID, matePosGuess-5) ) {						
						continue;
					}
				}
				
				if( nshot>100) { // pileup hot spot for bad alignments - bail on this pair
					if (ba1.Position>(NbadPos+1000)) { 
						cerr << "bam jump skip " << ba1.Name << "\t " << matePosGuess << "\t " << ba2.Position << "\t " << ba2.Length << endl;
						NbadPos=ba1.Position;
					}
					break;
				}
				
			}
			
		}
		
		if (findmate && (scan||(ba2.MapQuality>=Qmin) )) {
			break;
		}	
		
		if (findmate && (!scan) && (ba2.MapQuality<Qmin )) {
			// back up to next first alignment
			findmate=false;
			NmateFound=0;
			continue;
		}	
		
		
	}	
	
  if(NmateFound<0) {
		return(false);
	}
	
	NmateFound = BamBam2PairedRead(ba1, ba2, pr);
		
	// check consistent LF
	if ( (ba1.RefID==ba2.MateRefID) & properOrientation) {
		int LF1= set[0].Fraglength(pr);
		if (LF1!=LF)  {
			char i=char(1-ba1.IsReverseStrand()*1);
			char s1=FR[i];
			i=char(1-ba2.IsReverseStrand()*1);
			char s2=FR[i];
			//s+= FR[char(1-ma2.IsReverseStrand*1)];
			cerr	 << " fraglength problem " << LF1 << "\t" << LF << "\t" << s1 << s2 << "\t" << ba1.Name << endl;
		}
	}	
		
  return (NmateFound>0);		

}



//------------------------------------------------------------------------------
// convert a Bam read pair to one complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentPairSpecial( BamReader & ar1, C_pairedread & pr) 
{
	
	
	//------------------------------------------------------------------------------
	// select abberant pairs that meet Qmin criteria, include dangling (unmapped) ends
	//------------------------------------------------------------------------------
	
	// pair structures
	//C_pairedread pr1;	
	
	// Bam structure
	BamAlignment ba1,ba2;
	
	bool scan = SpannerMode==SPANNER_SCAN;
	//bool build = SpannerMode==SPANNER_BUILD;
	
	int QminSpecial=0;
	int LF=0;
	bool findmate=false;	
	int NmateFound=-1;
	string FR="FR";
	bool rev,revMate;
	bool properOrientation = false;	
	
	
	while ( ar1.GetNextAlignment(ba1) ) {
		
		ar1.GetNextAlignment(ba2);
		
		if (doneFrag[ba1.Name]) {
			continue;
		}
		
		doneFrag[ba1.Name]=true;
		
		// skip unmapped reads (another function in a module somewhere ...)
		if (ba1.RefID<0) {
			continue;
		}	
		
		// single ends not processed here
		if (!ba1.IsPaired()) {
  		continue;
		}
		
		
		if ( (!scan) && (ba1.MapQuality<QminSpecial) ) {
			continue;
		}	
		
				
		// skip proper pairs	=  abberant pair criteria 	
		LF = ba1.InsertSize;
		// SLX RP only  ??? 
		if (ba1.IsReverseStrand()) LF=-LF; 
		string rgid;
		ba1.GetReadGroup(rgid);
		unsigned int rgcode=libraries.ReadGroupID2Code[rgid];
		int LMlow = libraries.libmap[rgcode].LMlow;
		int LMhigh = libraries.libmap[rgcode].LMhigh;
		
		// orientation transform to Illumina FR
    
		rev=ba1.IsReverseStrand();
		revMate=ba1.IsMateReverseStrand();
		
		//454 -> IL
		//454 -> IL
		if ( MateMode==MATEMODE_454) {
			if (ba1.IsFirstMate()) {
				rev=!rev;			
			} else {
				revMate=!revMate;
			}
			// SOLiD -> IL
		} else if ( MateMode==MATEMODE_SOLID) {
			if (ba1.IsFirstMate()) {
				revMate=!revMate;
			} else {
				rev=!rev;			
			}
		}		
		
		// now in Illumina convention... check proper pair orientation
		if (ba1.RefID==ba1.MateRefID) {
			if (ba1.Position<ba1.MatePosition) {				
				properOrientation = (!ba1.IsReverseStrand())&&ba1.IsMateReverseStrand();
			}	else {
				properOrientation = ba1.IsReverseStrand()&&(!ba1.IsMateReverseStrand());
			}
		}
		
		//-------------------------------------------------------------------------------
		// skip proper pairs
		//-------------------------------------------------------------------------------
		if ( (LF>LMlow)&&(LF<LMhigh)&&properOrientation&&(!scan)) {					
			continue;
		}	
		
		// artifact check for mate mapped = this read (illumina problem)
		if  ( (ba1.RefID==ba1.MateRefID) && (ba1.Position==ba1.MatePosition) ) { 
			continue;
		}
		
		//------------------------------------------------------------------------------
		// go mate hunting
		//------------------------------------------------------------------------------		
		findmate=false;		
		
		
		// check read name should be the same but not the same read (depend on First/Second mate flags). 
		if ( ( ba1.Name.compare(ba2.Name)==0 )  )  { 
			findmate=true;			
			NmateFound++;
			if (scan) break;
			if (ba2.MapQuality>=QminSpecial )  break;
			// low quality reads 
			NmateFound=0;
		}										
		
	}	
	
  if(NmateFound<0) {
		return(false);
	}
	
	NmateFound = BamSpecial2PairedRead(ba1, ba2, pr);
	
	return (NmateFound>0);		
}


/*
 
 //------------------------------------------------------------------------------
 // convert pair of Bam read records to Mosaik aligned pair record 
 // returns:
 //			0 = only end mapped
 //          1 = first of pair 
 //          2 = second in pair
 //          <1 = no mapping
 //------------------------------------------------------------------------------
 int  C_pairedfiles::Bam2MosaikPair(BamAlignment & ba1, BamAlignment & ba2, Mosaik::Alignment & ma1, 
 Mosaik::Alignment & ma2) 
 {   
 
 if ((BamZ==1)&&(ba1.MateRefID>=0))  {
 int pe = BamZA2Mosaik(ba1,ma1,ma2);
 return pe;
 }
 
 ma1.ReferenceBegin=0;
 ma1.ReferenceEnd=0;
 ma1.ReferenceIndex=0;
 ma1.Quality=0;
 ma1.QueryBegin=0;
 ma1.QueryEnd=0;
 ma1.IsReverseStrand=false;
 ma1.IsPairedEnd=false;
 ma1.IsResolvedAsPair=false;
 ma1.NumMismatches=0;
 ma1.IsFirstMate=false;
 ma2.ReferenceBegin=0;
 ma2.ReferenceEnd=0;
 ma1.ReferenceIndex=0;
 ma2.Quality=0;
 ma2.QueryBegin=0;
 ma2.QueryEnd=0;
 ma2.IsReverseStrand=false;
 ma2.IsPairedEnd=false;
 ma2.IsResolvedAsPair=false;
 ma2.NumMismatches=0;
 ma2.IsFirstMate=false;
 
 int pe=0;
 
 if ( ba1.IsMapped() ) {
 
 ma1.ReferenceBegin=ba1.Position;
 int len=BamCigarData2Len(ba1.CigarData,false);
 ma1.ReferenceEnd=ma1.ReferenceBegin+len; //ba1.Length-1;
 ma1.ReferenceIndex=ba1.RefID;
 ma1.QueryBegin=1;
 ma1.QueryEnd=BamCigarData2Len(ba1.CigarData,true); //ba1.QueryBases.size();
 ma1.Quality=ba1.MapQuality;
 ma1.IsReverseStrand=ba1.IsReverseStrand();
 ma1.IsPairedEnd=ba1.IsPaired();
 ma2.IsPairedEnd=ba1.IsPaired();
 ma1.IsResolvedAsPair=ba1.IsProperPair();
 ma1.IsFirstMate=ba1.IsFirstMate();
 
 // convert ReadGroupID to readGroupCode (Mosaik/Spanner) 
 string ReadGroupID;
 ba1.GetReadGroup(ReadGroupID);
 
 // all libraries stored in every set for this map to work 
 ma1.ReadGroupCode=ReadGroupID2Code[ReadGroupID];
 ma2.ReadGroupCode=ma1.ReadGroupCode;
 
 // mismatches NM
 string tag="NMs";
 int itag;
 if (ba1.GetTag(tag,itag)) {
 ma1.NumMismatches=(itag>=0? itag: 0);
 } else { 
 int mm=BamCigarData2mm(ba1.CigarData);
 tag="MDs";
 string MD;
 if (ba1.GetTag(tag,MD)) {						
 mm+=getMDMismatchCount(MD);
 }
 ma1.NumMismatches=mm;
 }
 
 }
 
 if ( ba1.IsMapped() & ba2.IsMapped() ) {
 
 ma2.ReferenceBegin=ba2.Position;
 int len=BamCigarData2Len(ba2.CigarData,false);
 ma2.ReferenceEnd=ma2.ReferenceBegin+len; //ba1.Length-1;
 //ma2.ReferenceEnd=ma2.ReferenceBegin+ba2.AlignedBases.size()-1; //ba1.Length-1;
 ma2.ReferenceIndex=ba2.RefID;
 ma2.QueryBegin=1;
 ma2.QueryEnd=BamCigarData2Len(ba2.CigarData,true); //ba2.QueryBases.size();
 ma2.Quality=ba2.MapQuality;
 ma2.IsReverseStrand=ba2.IsReverseStrand();
 ma2.IsFirstMate=ba2.IsFirstMate();
 
 // consistency check
 if (ba1.MateRefID != ba2.RefID) {
 cerr << " bad bam mate: " << ba1.Name << "\t" << ba1.MateRefID << "\t" << ba2.RefID << "\n";
 return(-1);
 }
 if (ba1.MatePosition != ba2.Position) {
 cerr << " bad bam mate: " << ba1.Name << "\t" << ba1.MatePosition << "\t" << ba2.Position << "\n";
 return(-2);
 }
 if (ba1.IsMateReverseStrand() != ba2.IsReverseStrand() )  {
 cerr << " bad bam mate: " << ba1.Name << "\t" << ba1.IsMateReverseStrand() << "\t" << ba2.IsReverseStrand() << "\n";
 return(-3);
 }
 
 // convert ReadGroupID to readGroupCode (Mosaik/Spanner) 
 string ReadGroupID;
 ba1.GetReadGroup(ReadGroupID);
 
 // mismatches NM
 string tag="NMs";
 int itag;
 if (ba2.GetTag(tag,itag)) {
 ma2.NumMismatches=(itag>=0? itag: 0);
 } else { 
 int mm=BamCigarData2mm(ba2.CigarData);
 tag="MDs";
 string MD;
 if (ba2.GetTag(tag,MD)) {						
 mm+=getMDMismatchCount(MD);
 }
 ma2.NumMismatches=mm;
 }
 }	
 
 // return index to first end in machine order
 if (ba1.IsFirstMate() && ba1.IsMapped() && ba2.IsMapped()) {
 pe=1;
 } else 	if (ba2.IsFirstMate() && ba1.IsMapped() && ba2.IsMapped()) {
 pe=2;
 }
 
 return pe;
 }
 
 
 
 //------------------------------------------------------------------------------
 // convert one Bam record to Mosaik aligned pair record from position sorted bam
 // returns:
 //			0 = no mapped reads (error)
 //          1 = single end  
 //          2 = pair
 //------------------------------------------------------------------------------
 int  C_pairedfiles::Bam2MosaikScan(BamAlignment & ba1,  Mosaik::Alignment & ma1, 
 Mosaik::Alignment & ma2) 
 {   
 
 if ((BamZ==1)&&(ba1.MateRefID>=0))  {
 int pe = BamZA2Mosaik(ba1,ma1,ma2);
 return pe;
 }
 
 ma1.ReferenceBegin=0;
 ma1.ReferenceEnd=0;
 ma1.ReferenceIndex=0;
 ma1.Quality=0;
 ma1.QueryBegin=0;
 ma1.QueryEnd=0;
 ma1.IsReverseStrand=false;
 ma1.IsPairedEnd=false;
 ma1.IsResolvedAsPair=false;
 ma1.NumMismatches=0;
 ma1.IsFirstMate=false;
 ma2.ReferenceBegin=0;
 ma2.ReferenceEnd=0;
 ma1.ReferenceIndex=0;
 ma2.Quality=0;
 ma2.QueryBegin=0;
 ma2.QueryEnd=0;
 ma2.IsReverseStrand=false;
 ma2.IsPairedEnd=false;
 ma2.IsResolvedAsPair=false;
 ma2.NumMismatches=0;
 ma2.IsFirstMate=false;
 
 int pe=0;
 
 if ( ba1.IsMapped() ) {
 
 pe=1;
 
 ma1.ReferenceBegin=ba1.Position;
 int len=BamCigarData2Len(ba1.CigarData,false);
 ma1.ReferenceEnd=ma1.ReferenceBegin+len; //ba1.Length-1;
 //ma1.ReferenceEnd=ma1.ReferenceBegin+ba1.AlignedBases.size()-1; //ba1.Length-1;
 ma1.ReferenceIndex=ba1.RefID;
 ma1.QueryBegin=1;
 len=BamCigarData2Len(ba1.CigarData,true);
 ma1.QueryEnd=len; //ba1.QueryBases.size();
 ma1.Quality=ba1.MapQuality;
 if (ba1.MapQuality>100) ma1.Quality=100; 
 ma1.IsReverseStrand=ba1.IsReverseStrand();
 ma1.IsPairedEnd=ba1.IsPaired();
 ma1.IsResolvedAsPair=ba1.IsProperPair();
 ma1.IsFirstMate=ba1.IsFirstMate();
 
 // convert ReadGroupID to readGroupCode (Mosaik/Spanner) 
 string ReadGroupID;
 ba1.GetReadGroup(ReadGroupID);
 // all libraries stored in every set for this map to work 
 ma1.ReadGroupCode=ReadGroupID2Code[ReadGroupID];
 
 // mismatches NM
 string tag="NMs";
 int itag;
 if (ba1.GetTag(tag,itag)) {
 ma1.NumMismatches=(itag>=0? itag: 0);
 } else { 
 int mm=BamCigarData2mm(ba1.CigarData);
 tag="MDs";
 string MD;
 if (ba1.GetTag(tag,MD)) {						
 mm+=getMDMismatchCount(MD);
 }
 ma1.NumMismatches=mm;
 }
 
 
 if (ma1.IsPairedEnd && ba1.IsMateMapped()) {
 
 pe=2;
 
 ma2.ReadGroupCode=ma1.ReadGroupCode;
 ma2.IsPairedEnd=ba1.IsPaired();
 ma2.IsFirstMate=!ba1.IsFirstMate();
 ma2.IsReverseStrand=ba1.IsMateReverseStrand();
 ma2.ReferenceIndex=ba1.MateRefID;
 ma2.ReferenceBegin=ba1.MatePosition;
 
 
 // awful hack for scan  
 ma2.Quality=ma1.Quality;
 ma2.NumMismatches=ma1.NumMismatches;
 
 ma2.QueryBegin=1;
 // no info
 ma2.QueryEnd=0; // ma2.ReferenceEnd-ma2.ReferenceStart;
 //ma2.Quality=0; //ma1.Quality;
 
 // 454 & SOLiD pair orientation gymnastics
 if (MateMode==MATEMODE_454) {  // 454
 if (ma1.IsFirstMate) {
 ma1.IsReverseStrand=!ma1.IsReverseStrand;
 } else {
 ma2.IsReverseStrand=!ma2.IsReverseStrand;
 }
 } else if (MateMode==MATEMODE_SOLID) { // SOLiD
 if (!ma1.IsFirstMate) {
 ma1.IsReverseStrand=!ma1.IsReverseStrand;
 } else {
 ma2.IsReverseStrand=!ma2.IsReverseStrand;
 }
 }
 
 // FR orientation ... not for 454, SOLiD....
 bool FR = (!ma1.IsReverseStrand)&&(ma2.IsReverseStrand)&&(ba1.InsertSize>0);
 if  (FR) {
 ma2.ReferenceEnd=ma1.ReferenceBegin+ba1.InsertSize;
 } else {
 // hack - mate2 is like mate1 for scan....
 ma2.ReferenceEnd=ma2.ReferenceBegin+ma1.ReferenceEnd-ma1.ReferenceBegin;
 }
 
 }
 
 }	
 
 return pe;
 }
 
 
 //------------------------------------------------------------------------------
 // convert pair of Bam read records to Mosaik aligned pair record 
 // returns:
 //			0 = only end mapped
 //          1 = first of pair 
 //          2 = second in pair
 //          <1 = no mapping
 //------------------------------------------------------------------------------
 int  C_pairedfiles::Bam2MosaikRead(BamAlignment & ba1,  Mosaik::Alignment & ma1) 
 {   
 ma1.ReferenceBegin=0;
 ma1.ReferenceEnd=0;
 ma1.ReferenceIndex=0;
 ma1.Quality=0;
 ma1.QueryBegin=0;
 ma1.QueryEnd=0;
 ma1.IsReverseStrand=false;
 ma1.IsPairedEnd=false;
 ma1.IsResolvedAsPair=false;
 ma1.NumMismatches=0;
 ma1.IsFirstMate=false;
 
 int pe=-1;
 
 if ( ba1.IsMapped() ) {
 
 ma1.ReferenceBegin=ba1.Position;
 ma1.ReferenceIndex=ba1.RefID;
 ma1.QueryBegin=1;
 if (ba1.AlignedBases.size()>0) {  
 ma1.ReferenceEnd=ma1.ReferenceBegin+ba1.AlignedBases.size()-1; //ba1.Length-1;
 ma1.QueryEnd=ba1.QueryBases.size();
 }
 ma1.Quality=ba1.MapQuality;
 ma1.IsReverseStrand=ba1.IsReverseStrand();
 ma1.IsPairedEnd=ba1.IsPaired();
 ma1.IsFirstMate=ba1.IsFirstMate();
 if(ba1.IsMateMapped()&&ma1.IsPairedEnd) {
 // cerr << " Problem in Bam2MosaikRead: " << ba1.Name << "\n";
 //return(pe);
 }
 
 //if (opt>0)  { 
 // convert ReadGroupID to readGroupCode (Mosaik/Spanner) 
 string ReadGroupID;
 ba1.GetReadGroup(ReadGroupID);
 
 // all libraries stored in every set for this map to work 
 ma1.ReadGroupCode=ReadGroupID2Code[ReadGroupID];
 //}
 
 // mismatches NM
 string tag="NMs";
 int itag;
 if (ba1.GetTag(tag,itag)) {
 ma1.NumMismatches=(itag>=0? itag: 0);
 } else { 
 int mm=BamCigarData2mm(ba1.CigarData);
 tag="MDs";
 string MD;
 if (ba1.GetTag(tag,MD)) {						
 mm+=getMDMismatchCount(MD);
 }
 ma1.NumMismatches=mm;
 }
 
 pe=0;
 
 }
 
 
 return pe;
 }
 
 
 //------------------------------------------------------------------------------
 // convert Bam read records to Mosaik aligned pair record 
 // returns:
 //          0 = single unique mapping
 //          1 = first of pair in unique mapping
 //          2 = second of pair in unique mapping
 //         -1 = no mapping
 //          0+10*NM = single-end mapped read with NM mappings
 //          1+10*NM = first of pair with NM mappings
 //          2+10*NM = second of pair with NM mappings
 //------------------------------------------------------------------------------
 int  C_pairedfiles::Bam2Mosaik(BamAlignment & ba1, Mosaik::Alignment & ma1, 
 Mosaik::Alignment & ma2) 
 {      
 if ((BamZ==1)&&(ba1.IsPaired())&&(ba1.MateRefID>=0) ) {
 int pe = BamZA2Mosaik(ba1,ma1,ma2);
 return pe;
 }
 
 if ( ba1.IsMapped() ) {
 // check all this with Michael
 ma1.ReferenceBegin=ba1.Position;
 ma1.ReferenceEnd=ma1.ReferenceBegin+ba1.AlignedBases.size()-1; //ba1.Length-1;
 ma1.ReferenceIndex=ba1.RefID;
 ma1.QueryBegin=1;
 ma1.QueryEnd=ba1.QueryBases.size();
 ma1.Quality=ba1.MapQuality;
 ma1.IsReverseStrand=ba1.IsReverseStrand();
 ma1.IsPairedEnd=ba1.IsPaired();
 ma2.IsPairedEnd=ba1.IsPaired();
 ma1.IsResolvedAsPair=ba1.IsProperPair();
 ma2.IsResolvedAsPair=ba1.IsProperPair();
 
 // convert ReadGroupID to readGroupCode (Mosaik/Spanner) 
 string ReadGroupID;
 ba1.GetReadGroup(ReadGroupID);
 
 // all libraries stored in every set for this map to work 
 ma1.ReadGroupCode=ReadGroupID2Code[ReadGroupID];
 ma2.ReadGroupCode=ma1.ReadGroupCode;
 
 // mismatches NM
 string tag="NMs";
 int itag;
 if (ba1.GetTag(tag,itag)) {
 ma1.NumMismatches=(itag>=0? itag: 0);
 } else { 
 int mm=BamCigarData2mm(ba1.CigarData);
 tag="MDs";
 string MD;
 if (ba1.GetTag(tag,MD)) {						
 mm+=getMDMismatchCount(MD);
 }
 ma1.NumMismatches=mm;
 }
 
 // non-unique alignments - H0+H1+H2 
 int nmap=0;
 if (!ma1.IsResolvedAsPair) {		  
 tag = "H0s";
 if (ba1.GetTag(tag,itag) ) { 
 if (itag>=0) nmap+=itag;
 }
 tag = "H1s";
 if (ba1.GetTag(tag,itag) ) { 
 if (itag>=0) nmap+=itag;
 }
 tag = "H2s";
 if (ba1.GetTag(tag,itag) ) { 
 if (itag>=0) nmap+=itag;
 }
 }
 
 // mate mapping Q in pair 
 tag = "MQs";
 int Q2 = 0;
 if (ba1.GetTag(tag,itag) ) { 
 if (itag>=0) Q2=itag;
 }
 
 // set missing mate mismatches to 99 to flag as Bam mated end  
 // no optional flag for this
 ma2.NumMismatches=99;
 
 if ( ba1.IsMateMapped() & (ba1.MateRefID>=0) & (ba1.MatePosition>=0)) {
 
 
 ma2.ReferenceIndex=ba1.MateRefID;      // ID for reference sequence that mate was aligned to
 ma2.ReferenceBegin=ba1.MatePosition;   // position that mate was aligned to
 ma2.IsReverseStrand=ba1.IsMateReverseStrand();
 ma2.QueryBegin=1;	
 
 
 }
 
 // single end
 return 10*nmap;
 
 } else {
 
 // no ends
 ma1.ReferenceBegin=0;
 ma1.ReferenceEnd=0;
 ma1.ReferenceIndex=0;
 ma1.Quality=0;
 ma1.QueryBegin=0;
 ma1.QueryEnd=0;
 ma1.IsReverseStrand=false;
 ma1.IsPairedEnd=false;
 ma1.IsResolvedAsPair=false;
 ma1.NumMismatches=0;
 ma2.ReferenceBegin=0;
 ma2.ReferenceIndex=0;
 ma2.Quality=0;
 }
 return -1;
 }
 
 //------------------------------------------------------------------------------
 // convert Bam read record with ZA tag to Mosaik aligned pair record 
 // returns:
 //         -1 = error 
 //          0 = no reads
 //          1 = one read
 //          2 = both ends
 //------------------------------------------------------------------------------
 int  C_pairedfiles::BamZA2Mosaik(BamAlignment & ba1, Mosaik::Alignment & ma1, 
 Mosaik::Alignment & ma2) 
 {      
 // ZA tag pattern
 
 
 //RE2 patternZA("<(.);(\\d?);(\\d?);(.?);(\\d*);(.*?)>");
 RE2 patternZA("<(.);(\\d*?);(\\d*?);(.*?);(\\d*?);(.*?)>");
 
 
 // <@;Q1;Q2;Mob;# mappings;> <=;Q1;Q2;Mob;# mappings;CIGAR>
 // example: <@;0;0;;518;><=;0;0;;518;45M1I30M>
 
 
 // clear Mosaik structures
 ma1.ReferenceBegin=0;
 ma1.ReferenceEnd=0;
 ma1.ReferenceIndex=0;
 ma1.QueryBegin=0;
 ma1.QueryEnd=0;
 ma1.Quality=0;
 ma1.NumMismatches=0;
 ma2.ReferenceBegin=0;
 ma2.ReferenceEnd=0;
 ma2.ReferenceIndex=0;
 ma2.QueryBegin=0;
 ma2.QueryEnd=0;
 ma2.Quality=0;
 ma2.NumMismatches=0;
 
 int NZA=0;
 string tagZA = "ZAs";
 string tagNM = "NMs";
 string tagMD = "MDs";
 string ZA;
 if (!ba1.GetTag(tagZA,ZA) ) { 
 return -1;
 } 			  
 int q1,q2,nmap1;
 string this1,mob1,cig1,md1;
 string ReadGroupID1;
 
 StringPiece ZASP(ZA);    // Wrap a StringPiece around it
 int lenQ,lenR, mm;
 while (RE2::FindAndConsume(&ZASP,patternZA,&this1,&q1,&q2,&mob1,&nmap1,&cig1) ) {
 
 char t1 = this1[0];
 switch (t1) {
 
 case '@':
 
 NZA++;
 
 ma1.ReferenceIndex=ba1.RefID;
 // define start of query at start of mapped part
 ma1.QueryBegin=1;  
 // length of read mapping in query coordinates
 lenQ=BamCigarData2Len(ba1.CigarData,0);
 ma1.QueryEnd=ma1.QueryBegin+lenQ-1;
 
 ma1.ReferenceBegin=ba1.Position;
 // length of read mapping in reference coordinates
 lenR=BamCigarData2Len(ba1.CigarData,1);
 ma1.ReferenceEnd=ma1.ReferenceBegin+lenR-1; //ba1.Length-1;
 
 
 ma1.Quality=ba1.MapQuality;
 ma1.IsReverseStrand=ba1.IsReverseStrand();
 ma1.IsPairedEnd=ba1.IsPaired();
 ma1.IsFirstMate=ba1.IsFirstMate();
 
 // convert ReadGroupID to readGroupCode (Mosaik/Spanner) 
 
 ba1.GetReadGroup(ReadGroupID1);							
 // all libraries stored in every set for this map to work 
 ma1.ReadGroupCode=ReadGroupID2Code[ReadGroupID1];
 
 // mismatches NM
 if (ba1.GetTag(tagNM,mm)) {
 ma1.NumMismatches=(mm>=0? mm: 0);
 } else { 
 int mm=BamCigarData2mm(ba1.CigarData);
 string MD;
 if (ba1.GetTag(tagMD,MD)) {						
 mm+=getMDMismatchCount(MD);
 }
 ma1.NumMismatches=mm;
 }
 
 break;
 
 // mate
 case '=':  
 
 NZA++;
 
 ma2.ReadGroupCode=ma1.ReadGroupCode;
 ma2.IsPairedEnd=ba1.IsPaired();
 ma2.IsFirstMate=!ba1.IsFirstMate();
 ma2.IsReverseStrand=ba1.IsMateReverseStrand();
 ma2.ReferenceIndex=ba1.MateRefID;
 ma2.ReferenceBegin=ba1.MatePosition;	
 
 if (!getCigarLengths(cig1,lenQ,lenR,mm)) continue;
 
 ma2.Quality=q1;
 int mm=getCigarMismatchCount(cig1);
 // fix this with MD
 md1="";
 mm+=getMDMismatchCount(md1);
 ma2.NumMismatches=mm;
 
 // mate read starts with first mapped base
 ma2.QueryBegin=1;
 ma2.QueryEnd=ma2.QueryBegin+lenQ-1;
 
 ma2.ReferenceBegin=ba1.MatePosition;	
 ma2.ReferenceEnd=ma2.ReferenceBegin+lenR-1;
 
 // 454 & SOLiD pair orientation gymnastics
 if (MateMode==MATEMODE_454) { // 454
 if (ma1.IsFirstMate) {
 ma1.IsReverseStrand=!ma1.IsReverseStrand;
 } else {
 ma2.IsReverseStrand=!ma2.IsReverseStrand;
 }
 } else if (MateMode==MATEMODE_SOLID) { // SOLiD
 if (!ma1.IsFirstMate) {
 ma1.IsReverseStrand=!ma1.IsReverseStrand;
 } else {
 ma2.IsReverseStrand=!ma2.IsReverseStrand;
 }
 }
 
 }
 }
 
 if (NZA<1) {
 cerr << "no ZA tag info \n";
 }
 
 return NZA;
 }
 

 
//------------------------------------------------------------------------------
// convert a Bam read pair to one complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentPair( BamMultiReader & ar1,BamMultiReader & ar2, Mosaik::AlignedRead & mr1) 
{
	
	if (BamZ==1) {  
		// use ZA tag to extract  all pair info from first end in bam 
		return (nextBamAlignmentScan( ar1, mr1));
	}
	
	if (BamZ==2) {  
		// scan all pairs from first end in bam 
		return (false) ; //nextBamAlignmentPairSortedByName( ar1, mr1));
	}
	
	if (scan) {  
		// scan all pairs from first end in bam 
		return (nextBamAlignmentScan( ar1, mr1));
	}
	
	
	
	//------------------------------------------------------------------------------
	// select abberant pairs that meet Qmin criteria, include dangling (unmapped) ends
	//------------------------------------------------------------------------------
	
	// Mosaik structures
	Mosaik::Alignment ma1,ma2;
	Mosaik::Alignment mz1,mz2;
	
	// Bam structure
	BamAlignment ba1,ba2;
	
	// clear Mosaik structures
	ma1.ReferenceBegin=0;
	ma1.ReferenceEnd=0;
	ma1.ReferenceIndex=0;
	ma1.QueryBegin=0;
	ma1.QueryEnd=0;
	ma1.Quality=0;
	ma1.IsPairedEnd=false;
	ma1.IsReverseStrand=false;
	ma2.ReferenceBegin=0;
	ma2.ReferenceEnd=0;
	ma2.ReferenceIndex=0;
	ma2.QueryBegin=0;
	ma2.QueryEnd=0;
	ma2.Quality=0;
	ma2.IsPairedEnd=false;
	ma2.IsReverseStrand=false;
	mr1.Mate1Alignments.clear();
	mr1.Mate2Alignments.clear();
	
	
	int LF=0;
	bool findmate=false;	
	int NmateFound=-1;
	string FR="FR";
	bool rev,revMate;
	bool properOrientation = false;	
	
	while ( ar1.GetNextAlignment(ba1) ) {
		
		if (doneFrag[ba1.Name]) {
			continue;
		}
		
		doneFrag[ba1.Name]=true;
		
		// skip unmapped reads (another function in a module somewhere ...)
		if (ba1.RefID<0) {
			continue;
		}	
		
		// single ends not processed here
		if (!ba1.IsPaired()) { 
			continue;
		}
		
		if ( (ba1.MapQuality<Qmin) && (!scan) ) {
			continue;
		}	
		
		// dangling end is processed
		if (!ba1.IsMateMapped() ) {
			NmateFound=1;  // set to mark dangling end
			break;
		}	
		
		//if (ba1.Position>=1130515) { 
		// cerr << endl; }
		
		// skip proper pairs	=  abberant pair criteria 	
		LF = ba1.InsertSize;
		// SLX RP only  ??? 
		if (ba1.IsReverseStrand()) LF=-LF; 
		string rgid;
		ba1.GetReadGroup(rgid);
		unsigned int rgcode=libraries.ReadGroupID2Code[rgid];
		int LMlow = libraries.libmap[rgcode].LMlow;
		int LMhigh = libraries.libmap[rgcode].LMhigh;
		
		// orientation transform to Illumina FR
		
		rev=ba1.IsReverseStrand();
		revMate=ba1.IsMateReverseStrand();
		
		//454 -> IL
		if ( MateMode==MATEMODE_454) {
			if (ba1.IsFirstMate()) {
				rev=!rev;			
			} else {
				revMate=!revMate;
			}
			// SOLiD -> IL
		} else if ( MateMode==MATEMODE_SOLID) {
			if (ba1.IsFirstMate()) {
				revMate=!revMate;
			} else {
				rev=!rev;			
			}
		}		
		
		
		if (ba1.RefID==ba1.MateRefID) {
			if (ba1.Position<ba1.MatePosition) {				
				properOrientation = (!ba1.IsReverseStrand())&&ba1.IsMateReverseStrand();
			}	else {
				properOrientation = ba1.IsReverseStrand()&&(!ba1.IsMateReverseStrand());
			}
		}
		
		//-------------------------------------------------------------------------------
		// skip proper pairs
		//-------------------------------------------------------------------------------
		if (ba1.IsMateMapped()&&(LF>LMlow)&&(LF<LMhigh)&&properOrientation) {					
			continue;
		}	
		
		// artifact check for mate mapped = this read (illumina problem)
		if  ( (ba1.RefID==ba1.MateRefID) && (ba1.Position==ba1.MatePosition) ) { 
			continue;
		}
		
		//------------------------------------------------------------------------------
		// go mate hunting
		//------------------------------------------------------------------------------		
		findmate=false;
		
		int matePosGuess = ba1.Position+ba1.InsertSize-1;
		
		if (ba1.RefID!=ba1.MateRefID)			{
			matePosGuess = ba1.MatePosition+5;
		}
		
		// jump to first alignment that overlaps mate
		if ( ar2.Jump(ba1.MateRefID, matePosGuess) ) {
			
			int nshot=0;
			
			while ( ar2.GetNextAlignment(ba2) ) {
				
				nshot++;
				
				if (ba2.RefID!=ba1.MateRefID) {
					Shots2Mate.push_back(nshot);
					// this set of bam files doesn't have mateRefID.... log this somewhere other than cerr, cout
					//cerr << "no mate RefID in bam:  want " << ba1.MateRefID << "\tget" << ba2.RefID << "\n";  // oops
					break;
				}
				
				// check read name should be the same but not the same read (depend on First/Second mate flags). 
				if ( ( ba1.Name.compare(ba2.Name)==0 ) && (ba1.IsFirstMate()!=ba2.IsFirstMate()) )  { 
					findmate=true;					
					NmateFound=2;
					Shots2Mate.push_back(nshot);
					//if (nskip>500)  { 
					//	cerr << "bam find after large skip " << nskip << "\n";  // oops
					//	cerr << ba1.Name << "\t " << matePosGuess << "\t " << ba2.Position << "\t " << ba2.Length << endl;
					//}
					break;
				}
				
				if (ba2.Position>ba1.MatePosition) {
					cerr << "bam jump too far " << ba1.MatePosition << "\t" << ba2.Position << "\n";  // oops
					// jump to first alignment that overlaps mate
					if ( ar2.Jump(ba1.MateRefID, matePosGuess-5) ) {						
						continue;
					}
				}
				
				if( nshot>10000) { // pileup hot spot for bad alignments - bail on this pair
					if (ba1.Position>(NbadPos+1000)) { 
						cerr << "bam jump skip " << ba1.Name << "\t " << matePosGuess << "\t " << ba2.Position << "\t " << ba2.Length << endl;
						NbadPos=ba1.Position;
					}
					break;
				}
				
			}
			
		}
		
		if (findmate && (scan||(ba2.MapQuality>=Qmin ) )) {
			break;
		}	
		
		if (findmate && (!scan)&& (ba2.MapQuality<Qmin )) {
			// back up to next first alignment
			findmate=false;
			NmateFound=0;
			continue;
		}	
		
		
	}	
	
	if(NmateFound<0) {
		return(false);
	}
	
	if (findmate ) {
		
		NmateFound = Bam2MosaikPair(ba1, ba2, ma1,ma2);
		
		// check consistent LF
		if ( (ba1.RefID==ba2.MateRefID) & properOrientation) {		
			int LF1=ma2.ReferenceEnd-ma1.ReferenceBegin-1;
			if (LF1!=LF)  {
				char i=char(1-ma1.IsReverseStrand*1);
				char s1=FR[i];
				i=char(1-ma2.IsReverseStrand*1);
				char s2=FR[i];
				//s+= FR[char(1-ma2.IsReverseStrand*1)];
				cerr	 << " fraglength problem " << LF1 << "\t" << LF << "\t" << s1 << s2 << "\t" << ba1.Name << endl;
			}
		}	
		
	} else {
		NmateFound = Bam2MosaikRead(ba1, ma1);
	}
	
	if (NmateFound>=0)  {
		
		mr1.Name=ba1.Name;
		mr1.ReadGroupCode=ma1.ReadGroupCode;
		
		// shortest read length
		// shortestReadLength=(shortestReadLength<ma1.QueryEnd) ? shortestReadLength: ma1.QueryEnd;
		
		switch (NmateFound) {
				
			case 1:
				mr1.Mate1Alignments.push_back(ma1);
				mr1.Mate2Alignments.push_back(ma2);
				mr1.IsPairedEnd=true;
				break;
				
			case 2:
				mr1.Mate1Alignments.push_back(ma2);
				mr1.Mate2Alignments.push_back(ma1);
				mr1.IsPairedEnd=true;
				break;
				
			case 0:
				mr1.Mate1Alignments.push_back(ma1);
				
		} 
		
		
	}  //else {
	//cerr << " problem finding mate \t" << ba1.Name << " in chrom " << ba1.MateRefID << "\n"; 
	//}
	
	return (true);		
}


//------------------------------------------------------------------------------
// convert a Bam read pair to one complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentPair( BamReader & ar1,BamReader & ar2, Mosaik::AlignedRead & mr1) 
{
	
	if (BamZ==1) {  
		// use ZA tag to extract  all pair info from first end in bam 
		return (nextBamAlignmentZA( ar1, mr1));
	}
	
	if (BamZ==2) {  
		// scan all pairs from first end in bam 
		return (nextBamAlignmentPairSortedByName( ar1, mr1));
	}
	
	if (scan) {  
		// scan all pairs from first end in bam 
		return (nextBamAlignmentScan( ar1, mr1));
	}
	
	
	
	//------------------------------------------------------------------------------
	// select abberant pairs that meet Qmin criteria, include dangling (unmapped) ends
	//------------------------------------------------------------------------------
	
	// Mosaik structures
	Mosaik::Alignment ma1,ma2;
	Mosaik::Alignment mz1,mz2;
	
	// Bam structure
	BamAlignment ba1,ba2;
	
	// clear Mosaik structures
	ma1.ReferenceBegin=0;
	ma1.ReferenceEnd=0;
	ma1.ReferenceIndex=0;
	ma1.QueryBegin=0;
	ma1.QueryEnd=0;
	ma1.Quality=0;
	ma1.IsPairedEnd=false;
	ma1.IsReverseStrand=false;
	ma2.ReferenceBegin=0;
	ma2.ReferenceEnd=0;
	ma2.ReferenceIndex=0;
	ma2.QueryBegin=0;
	ma2.QueryEnd=0;
	ma2.Quality=0;
	ma2.IsPairedEnd=false;
	ma2.IsReverseStrand=false;
	mr1.Mate1Alignments.clear();
	mr1.Mate2Alignments.clear();
	
	
	int LF=0;
	bool findmate=false;	
	int NmateFound=-1;
	string FR="FR";
	bool rev,revMate;
	bool properOrientation = false;	
	
	while ( ar1.GetNextAlignment(ba1) ) {
		
		if (doneFrag[ba1.Name]) {
			continue;
		}
		
		doneFrag[ba1.Name]=true;
		
		// skip unmapped reads (another function in a module somewhere ...)
		if (ba1.RefID<0) {
			continue;
		}	
		
		// single ends not processed here
		if (!ba1.IsPaired()) { 
			continue;
		}
		
		if ( (!scan)&&(ba1.MapQuality<Qmin) ) {
			continue;
		}	
		
		// dangling end is processed
		if (!ba1.IsMateMapped() ) {
			NmateFound=1;  // set to mark dangling end
			break;
		}	
		
		//if (ba1.Position>=1130515) { 
		// cerr << endl; }
		
		// skip proper pairs	=  abberant pair criteria 	
		LF = ba1.InsertSize;
		// SLX RP only  ??? 
		if (ba1.IsReverseStrand()) LF=-LF; 
		string rgid;
		ba1.GetReadGroup(rgid);
		unsigned int rgcode=libraries.ReadGroupID2Code[rgid];
		int LMlow = libraries.libmap[rgcode].LMlow;
		int LMhigh = libraries.libmap[rgcode].LMhigh;
		
		// orientation transform to Illumina FR
		
		rev=ba1.IsReverseStrand();
		revMate=ba1.IsMateReverseStrand();
		
		//454 -> IL
		if ( MateMode==MATEMODE_454) {
			if (ba1.IsFirstMate()) {
				rev=!rev;			
			} else {
				revMate=!revMate;
			}
			// SOLiD -> IL
		} else if ( MateMode==MATEMODE_SOLID) {
			if (ba1.IsFirstMate()) {
				revMate=!revMate;
			} else {
				rev=!rev;			
			}
		}		
		
		
		if (ba1.RefID==ba1.MateRefID) {
			if (ba1.Position<ba1.MatePosition) {				
				properOrientation = (!ba1.IsReverseStrand())&&ba1.IsMateReverseStrand();
			}	else {
				properOrientation = ba1.IsReverseStrand()&&(!ba1.IsMateReverseStrand());
			}
		}
		
		//-------------------------------------------------------------------------------
		// skip proper pairs
		//-------------------------------------------------------------------------------
		if (ba1.IsMateMapped()&&(LF>LMlow)&&(LF<LMhigh)&&properOrientation) {					
			continue;
		}	
		
		// artifact check for mate mapped = this read (illumina problem)
		if  ( (ba1.RefID==ba1.MateRefID) && (ba1.Position==ba1.MatePosition) ) { 
			continue;
		}
		
		//------------------------------------------------------------------------------
		// go mate hunting
		//------------------------------------------------------------------------------		
		findmate=false;
		
		int matePosGuess = ba1.Position+ba1.InsertSize-1;
		
		if (ba1.RefID!=ba1.MateRefID)			{
			matePosGuess = ba1.MatePosition+5;
		}
		
		// jump to first alignment that overlaps mate
		if ( ar2.Jump(ba1.MateRefID, matePosGuess) ) {
			
			int nshot=0;
			
			while ( ar2.GetNextAlignment(ba2) ) {
				
				nshot++;
				
				if (ba2.RefID!=ba1.MateRefID) {
					Shots2Mate.push_back(nshot);
					// this set of bam files doesn't have mateRefID.... log this somewhere other than cerr, cout
					//cerr << "no mate RefID in bam:  want " << ba1.MateRefID << "\tget" << ba2.RefID << "\n";  // oops
					break;
				}
				
				// check read name should be the same but not the same read (depend on First/Second mate flags). 
				if ( ( ba1.Name.compare(ba2.Name)==0 ) && (ba1.IsFirstMate()!=ba2.IsFirstMate()) )  { 
					findmate=true;					
					NmateFound=2;
					Shots2Mate.push_back(nshot);
					//if (nskip>500)  { 
					//	cerr << "bam find after large skip " << nskip << "\n";  // oops
					//	cerr << ba1.Name << "\t " << matePosGuess << "\t " << ba2.Position << "\t " << ba2.Length << endl;
					//}
					break;
				}
				
				if (ba2.Position>ba1.MatePosition) {
					cerr << "bam jump too far " << ba1.MatePosition << "\t" << ba2.Position << "\n";  // oops
					// jump to first alignment that overlaps mate
					if ( ar2.Jump(ba1.MateRefID, matePosGuess-5) ) {						
						continue;
					}
				}
				
				if( nshot>10000) { // pileup hot spot for bad alignments - bail on this pair
					if (ba1.Position>(NbadPos+1000)) { 
						cerr << "bam jump skip " << ba1.Name << "\t " << matePosGuess << "\t " << ba2.Position << "\t " << ba2.Length << endl;
						NbadPos=ba1.Position;
					}
					break;
				}
				
			}
			
		}
		
		if (findmate && (scan||(ba2.MapQuality>=Qmin) )) {
			break;
		}	
		
		if (findmate && (!scan) && (ba2.MapQuality<Qmin )) {
			// back up to next first alignment
			findmate=false;
			NmateFound=0;
			continue;
		}	
		
		
	}	
	
	if(NmateFound<0) {
		return(false);
	}
	
	if (findmate ) {
		
		NmateFound = Bam2MosaikPair(ba1, ba2, ma1,ma2);
		
		// check consistent LF
		if ( (ba1.RefID==ba2.MateRefID) & properOrientation) {		
			int LF1=ma2.ReferenceEnd-ma1.ReferenceBegin-1;
			if (LF1!=LF)  {
				char i=char(1-ma1.IsReverseStrand*1);
				char s1=FR[i];
				i=char(1-ma2.IsReverseStrand*1);
				char s2=FR[i];
				//s+= FR[char(1-ma2.IsReverseStrand*1)];
				cerr	 << " fraglength problem " << LF1 << "\t" << LF << "\t" << s1 << s2 << "\t" << ba1.Name << endl;
			}
		}	
		
	} else {
		NmateFound = Bam2MosaikRead(ba1, ma1);
	}
	
	if (NmateFound>=0)  {
		
		mr1.Name=ba1.Name;
		mr1.ReadGroupCode=ma1.ReadGroupCode;
		
		// shortest read length
		// shortestReadLength=(shortestReadLength<ma1.QueryEnd) ? shortestReadLength: ma1.QueryEnd;
		
		switch (NmateFound) {
				
			case 1:
				mr1.Mate1Alignments.push_back(ma1);
				mr1.Mate2Alignments.push_back(ma2);
				mr1.IsPairedEnd=true;
				break;
				
			case 2:
				mr1.Mate1Alignments.push_back(ma2);
				mr1.Mate2Alignments.push_back(ma1);
				mr1.IsPairedEnd=true;
				break;
				
			case 0:
				mr1.Mate1Alignments.push_back(ma1);
				
		} 
		
		
	}  //else {
	//cerr << " problem finding mate \t" << ba1.Name << " in chrom " << ba1.MateRefID << "\n"; 
	//}
	
	return (true);		
}




//------------------------------------------------------------------------------
// convert a Bam read pair to one semi-complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
//bool  C_pairedfiles::nextBamAlignmentScan( BamReader & ar1, Mosaik::AlignedRead & mr1) 
bool  C_pairedfiles::nextBamAlignmentScan( BamMultiReader & ar1, Mosaik::AlignedRead & mr1) 
{
	
	// Mosaik structures
	Mosaik::Alignment ma1,ma2;
	
	// Bam structure
	BamAlignment ba1;
	
	// clear Mosaik structures
	ma1.ReferenceBegin=0;
	ma1.ReferenceEnd=0;
	ma1.ReferenceIndex=0;
	ma1.QueryBegin=0;
	ma1.QueryEnd=0;
	ma1.Quality=0;
	ma1.IsPairedEnd=false;
	ma1.IsReverseStrand=false;
	ma2.ReferenceBegin=0;
	ma2.ReferenceEnd=0;
	ma2.ReferenceIndex=0;
	ma2.QueryBegin=0;
	ma2.QueryEnd=0;
	ma2.Quality=0;
 	ma2.IsPairedEnd=false;
	ma2.IsReverseStrand=false;
	mr1.Mate1Alignments.clear();
	mr1.Mate2Alignments.clear();
	
	int pe=-1;
	
	int ndone = 0;
	bool ok=false;
	while ( ar1.GetNextAlignment(ba1) ) {
		
		// skip unmapped reads (another program somewhere ...)
		if (ba1.RefID<0) {
			continue;
		}	
		if (!ba1.IsMapped() ) {
			continue;
		}	
		
		int LF = ba1.InsertSize;
		// SLX RP only
		if (ba1.IsReverseStrand()) LF=-LF; 
		
		// check if fragment done already
		ndone++;
		if (ba1.Name.size()<1) {			
			cerr << "empty read name: " << ba1.RefID+1 <<"\t"<< ba1.Position+1<<"\t"<< ba1.MateRefID+1<<"\t"<< ba1.MatePosition+1 <<"\t"<< ndone << endl;
			continue;
		}
		
		if (!doneFrag[ba1.Name]) {
			ok=true;
			break;
		}
		
	}
	
	// check out when GetNextAlignment returns false
	if (!ok) { 
		cerr << ba1.Name.c_str() << ba1.RefID << endl;
		return(ok);
	}
	
	if (ndone>1000) {
		cerr << "problem in nextBamAlignmentPair: " << ndone << "\n";
	}
	
	// mark this one as done
	doneFrag[ba1.Name]=true;
		
	pe = Bam2MosaikScan(ba1, ma1,ma2);
	
	if (pe>=0)  {
		
		mr1.Name=ba1.Name;
		mr1.ReadGroupCode=ma1.ReadGroupCode;
		
		// shortest read length
		//shortestReadLength=(shortestReadLength<ma1.QueryEnd) ? shortestReadLength: ma1.QueryEnd;
		
		switch (pe) {
				
			case 2:
				mr1.IsPairedEnd=true;
				if  (!ma1.IsReverseStrand) { 
					mr1.Mate1Alignments.push_back(ma1);
					mr1.Mate2Alignments.push_back(ma2);
				} else {
					mr1.Mate1Alignments.push_back(ma2);
					mr1.Mate2Alignments.push_back(ma1);
				}
				break;
				
			case 1:
				mr1.Mate1Alignments.push_back(ma1);
				
		} 
		
		
	} else {
		cerr << ba1.Name.c_str() << ba1.RefID << endl;
	}
	
	// trap this 
	if (mr1.Name.size()==0) {
		cerr << ba1.Name.c_str() << ba1.RefID << endl;
	}		
	
	return (pe>=0);		
}


//------------------------------------------------------------------------------
// convert a Bam read pair to one semi-complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentScan( BamReader & ar1, Mosaik::AlignedRead & mr1) 
{
	
	// Mosaik structures
	Mosaik::Alignment ma1,ma2;
	
	// Bam structure
	BamAlignment ba1;
	
	// clear Mosaik structures
	ma1.ReferenceBegin=0;
	ma1.ReferenceEnd=0;
	ma1.ReferenceIndex=0;
	ma1.QueryBegin=0;
	ma1.QueryEnd=0;
	ma1.Quality=0;
	ma1.IsPairedEnd=false;
	ma1.IsReverseStrand=false;
	ma2.ReferenceBegin=0;
	ma2.ReferenceEnd=0;
	ma2.ReferenceIndex=0;
	ma2.QueryBegin=0;
	ma2.QueryEnd=0;
	ma2.Quality=0;
 	ma2.IsPairedEnd=false;
	ma2.IsReverseStrand=false;
	mr1.Mate1Alignments.clear();
	mr1.Mate2Alignments.clear();
	
	int pe=-1;
	
	int ndone = 0;
	bool ok=false;
	
	bool rev,revMate,skip;
	bool properOrientation = false;
	
	while ( ar1.GetNextAlignment(ba1) ) {
		
		// skip unmapped reads (another program somewhere ...)
		if (ba1.RefID<0) {
			continue;
		}	
		if (!ba1.IsMapped() ) {
			continue;
		}	
		
		
		int LF = ba1.InsertSize;
		// SLX RP only
		if (ba1.IsReverseStrand()) LF=-LF; 
		
		// check if fragment done already
		ndone++;
		if (ba1.Name.size()<1) {			
			cerr << "empty read name: " << ba1.RefID+1 <<"\t"<< ba1.Position+1<<"\t"<< ba1.MateRefID+1<<"\t"<< ba1.MatePosition+1 <<"\t"<< ndone << endl;
			continue;
		}

		//-------------------------------------------------------------------------------
		// skip proper pairs to save span file space from  build
		//-------------------------------------------------------------------------------
		skip=false;
  	if (!scan) { 
	
			string rgid;
			ba1.GetReadGroup(rgid);
			unsigned int rgcode=libraries.ReadGroupID2Code[rgid];
			int LMlow = libraries.libmap[rgcode].LMlow;
			int LMhigh = libraries.libmap[rgcode].LMhigh;
			
			// orientation transform to Illumina FR			
			rev=ba1.IsReverseStrand();
			revMate=ba1.IsMateReverseStrand();
			
			//454 -> IL
			if ( MateMode==MATEMODE_454) {
				if (ba1.IsFirstMate()) {
					rev=!rev;			
				} else {
					revMate=!revMate;
				}
				// SOLiD -> IL
			} else if ( MateMode==MATEMODE_SOLID) {
				if (ba1.IsFirstMate()) {
					revMate=!revMate;
				} else {
					rev=!rev;			
				}
			}		
			
			if (ba1.RefID==ba1.MateRefID) {
				if (ba1.Position<ba1.MatePosition) {				
					properOrientation = (!rev)&&revMate;
				}	else {
					properOrientation = rev&&(!revMate);
				}
			}
			
			if (ba1.IsMateMapped()&&(LF>LMlow)&&(LF<LMhigh)&&properOrientation) {					
				skip=true;
			}	
			
		}

		if (skip) continue;
		
		if (!doneFrag[ba1.Name]) {
			ok=true;
			break;
		}
						
	}
	
	// check out when GetNextAlignment returns false
	if (!ok) { 
		cerr << ba1.Name.c_str() << ba1.RefID << endl;
		return(ok);
	}
	
	if (ndone>1000) {
		cerr << "problem in nextBamAlignmentPair: " << ndone << "\n";
	}
	
	// mark this one as done
	doneFrag[ba1.Name]=true;
	
	pe = Bam2MosaikScan(ba1, ma1,ma2);
	
	if (pe>=0)  {
		
		mr1.Name=ba1.Name;
		mr1.ReadGroupCode=ma1.ReadGroupCode;
				
		switch (pe) {
				
			case 2:
				mr1.IsPairedEnd=true;
				if  (!ma1.IsReverseStrand) { 
					mr1.Mate1Alignments.push_back(ma1);
					mr1.Mate2Alignments.push_back(ma2);
				} else {
					mr1.Mate1Alignments.push_back(ma2);
					mr1.Mate2Alignments.push_back(ma1);
				}
				break;
				
			case 1:
				mr1.Mate1Alignments.push_back(ma1);
				
		} 
		
		
	} else {
		cerr << "bad read \t";
		cerr << ba1.Name.c_str() << ba1.RefID << endl;
	}
	
	// trap this 
	if (mr1.Name.size()==0) {
		cerr << "baf read name " << ba1.Name.c_str() << ba1.RefID << endl;
	}		
	
	return (pe>=0);		
}



//------------------------------------------------------------------------------
// convert a Bam read pair to one semi-complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentZA( BamReader & ar1, Mosaik::AlignedRead & mr1) 
{
	
	// Mosaik structures
	Mosaik::Alignment ma1,ma2;
	
	// Bam structure
	BamAlignment ba1;
	
	// clear Mosaik structures
	ma1.ReferenceBegin=0;
	ma1.ReferenceEnd=0;
	ma1.ReferenceIndex=0;
	ma1.QueryBegin=0;
	ma1.QueryEnd=0;
	ma1.Quality=0;
	ma1.IsPairedEnd=false;
	ma1.IsReverseStrand=false;
	ma2.ReferenceBegin=0;
	ma2.ReferenceEnd=0;
	ma2.ReferenceIndex=0;
	ma2.QueryBegin=0;
	ma2.QueryEnd=0;
	ma2.Quality=0;
 	ma2.IsPairedEnd=false;
	ma2.IsReverseStrand=false;
	mr1.Mate1Alignments.clear();
	mr1.Mate2Alignments.clear();
	
	int pe=-1;
	
	int ndone = 0;
	bool ok=false;
	
	bool rev,revMate,skip;
	bool properOrientation = false;
	
	while ( ar1.GetNextAlignment(ba1) ) {
		
		// skip unmapped reads (another program somewhere ...)
		if (ba1.RefID<0) {
			continue;
		}	
		if (!ba1.IsMapped() ) {
			continue;
		}	
		
		
		int LF = ba1.InsertSize;
		// SLX RP only
		if (ba1.IsReverseStrand()) LF=-LF; 
		
		// check if fragment done already
		ndone++;
		if (ba1.Name.size()<1) {			
			cerr << "empty read name: " << ba1.RefID+1 <<"\t"<< ba1.Position+1<<"\t"<< ba1.MateRefID+1<<"\t"<< ba1.MatePosition+1 <<"\t"<< ndone << endl;
			continue;
		}
		
		//-------------------------------------------------------------------------------
		// skip proper pairs to save span file space from  build
		//-------------------------------------------------------------------------------
		skip=false;
  	if (!scan) { 
			
			string rgid;
			ba1.GetReadGroup(rgid);
			unsigned int rgcode=libraries.ReadGroupID2Code[rgid];
			int LMlow = libraries.libmap[rgcode].LMlow;
			int LMhigh = libraries.libmap[rgcode].LMhigh;
			
			// orientation transform to Illumina FR			
			rev=ba1.IsReverseStrand();
			revMate=ba1.IsMateReverseStrand();
			
			//454 -> IL
			if ( MateMode==MATEMODE_454) {
				if (ba1.IsFirstMate()) {
					rev=!rev;			
				} else {
					revMate=!revMate;
				}
				// SOLiD -> IL
			} else if ( MateMode==MATEMODE_SOLID) {
				if (ba1.IsFirstMate()) {
					revMate=!revMate;
				} else {
					rev=!rev;			
				}
			}					
			
			if (ba1.RefID==ba1.MateRefID) {
				if (ba1.Position<ba1.MatePosition) {				
					properOrientation = (!rev)&&revMate;
				}	else {
					properOrientation = rev&&(!revMate);
				}
			}
			
			if (ba1.IsMateMapped()&&(LF>LMlow)&&(LF<LMhigh)&&properOrientation) {					
				skip=true;
			}	
			
		}
		
		if (skip) continue;
		
		if (!doneFrag[ba1.Name]) {
			ok=true;
			break;
		}
		
	}
	
	// check out when GetNextAlignment returns false
	if (!ok) { 
		cerr << ba1.Name.c_str() << ba1.RefID << endl;
		return(ok);
	}
	
	if (ndone>1000) {
		cerr << "problem in nextBamAlignmentPair: " << ndone << "\n";
	}
	
	// mark this one as done
	doneFrag[ba1.Name]=true;
	
	pe = BamZA2Mosaik(ba1, ma1,ma2);
	
	if (pe>=0)  {
		
		mr1.Name=ba1.Name;
		mr1.ReadGroupCode=ma1.ReadGroupCode;
		
		switch (pe) {
				
			case 2:
				mr1.IsPairedEnd=true;
				if  (!ma1.IsFirstMate) { 
					mr1.Mate1Alignments.push_back(ma1);
					mr1.Mate2Alignments.push_back(ma2);
				} else {
					mr1.Mate1Alignments.push_back(ma2);
					mr1.Mate2Alignments.push_back(ma1);
				}
				break;
				
			case 1:
				mr1.Mate1Alignments.push_back(ma1);
				
		} 
		
		
	} else {
		cerr << "bad read \t";
		cerr << ba1.Name.c_str() << ba1.RefID << endl;
	}
	
	// trap this 
	if (mr1.Name.size()==0) {
		cerr << "baf read name " << ba1.Name.c_str() << ba1.RefID << endl;
	}		
	
	return (pe>=0);		
}






//------------------------------------------------------------------------------
// convert a Bam read pair to one complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentPairSortedByName( BamReader & ar1, Mosaik::AlignedRead & mr1) 
{
	
	
	//------------------------------------------------------------------------------
	// select abberant pairs that meet Qmin criteria, include dangling (unmapped) ends
	//------------------------------------------------------------------------------
	
	// Mosaik structures
	Mosaik::Alignment ma1,ma2;
	Mosaik::Alignment mz1,mz2;
	
	// Bam structure
	BamAlignment ba1,ba2;
	
	// clear Mosaik structures
	ma1.ReferenceBegin=0;
	ma1.ReferenceEnd=0;
	ma1.ReferenceIndex=0;
	ma1.QueryBegin=0;
	ma1.QueryEnd=0;
	ma1.Quality=0;
	ma1.IsPairedEnd=false;
	ma1.IsReverseStrand=false;
	ma2.ReferenceBegin=0;
	ma2.ReferenceEnd=0;
	ma2.ReferenceIndex=0;
	ma2.QueryBegin=0;
	ma2.QueryEnd=0;
	ma2.Quality=0;
 	ma2.IsPairedEnd=false;
	ma2.IsReverseStrand=false;
	mr1.Mate1Alignments.clear();
	mr1.Mate2Alignments.clear();
	
	
	int LF=0;
	bool findmate=false;	
	int NmateFound=-1;
	string FR="FR";
	bool rev,revMate;
	bool properOrientation = false;	
	
	
	while ( ar1.GetNextAlignment(ba1) ) {
		
		if (doneFrag[ba1.Name]) {
			continue;
		}
		
		doneFrag[ba1.Name]=true;
		
		// skip unmapped reads (another function in a module somewhere ...)
		if (ba1.RefID<0) {
			continue;
		}	
		
		// single ends not processed here
		if (!ba1.IsPaired()) {
  		continue;
		}
		
		
		// dangling end is processed
		if (ba1.IsMateMapped() ) {
			ar1.GetNextAlignment(ba2);
		} else {
			NmateFound=1;  // set to mark dangling end
			break;
		}	

		if ( (!scan) && (ba1.MapQuality<Qmin) ) {
			continue;
		}	
		
		
		//if (ba1.Position>=1130515) { 
		// cerr << endl; }
		
		// skip proper pairs	=  abberant pair criteria 	
		LF = ba1.InsertSize;
		// SLX RP only  ??? 
		if (ba1.IsReverseStrand()) LF=-LF; 
		string rgid;
		ba1.GetReadGroup(rgid);
		unsigned int rgcode=libraries.ReadGroupID2Code[rgid];
		int LMlow = libraries.libmap[rgcode].LMlow;
		int LMhigh = libraries.libmap[rgcode].LMhigh;
		
		// orientation transform to Illumina FR
    
		rev=ba1.IsReverseStrand();
		revMate=ba1.IsMateReverseStrand();
		
		//454 -> IL
		//454 -> IL
		if ( MateMode==MATEMODE_454) {
			if (ba1.IsFirstMate()) {
				rev=!rev;			
			} else {
				revMate=!revMate;
			}
			// SOLiD -> IL
		} else if ( MateMode==MATEMODE_SOLID) {
			if (ba1.IsFirstMate()) {
				revMate=!revMate;
			} else {
				rev=!rev;			
			}
		}		
		
		// now in Illumina convention... check proper pair orientation
		if (ba1.RefID==ba1.MateRefID) {
			if (ba1.Position<ba1.MatePosition) {				
				properOrientation = (!ba1.IsReverseStrand())&&ba1.IsMateReverseStrand();
			}	else {
				properOrientation = ba1.IsReverseStrand()&&(!ba1.IsMateReverseStrand());
			}
		}
		
		//-------------------------------------------------------------------------------
		// skip proper pairs
		//-------------------------------------------------------------------------------
		if ( (LF>LMlow)&&(LF<LMhigh)&&properOrientation&&(!scan)) {					
			continue;
		}	
		
		// artifact check for mate mapped = this read (illumina problem)
		if  ( (ba1.RefID==ba1.MateRefID) && (ba1.Position==ba1.MatePosition) ) { 
			continue;
		}
		
		//------------------------------------------------------------------------------
		// go mate hunting
		//------------------------------------------------------------------------------		
		findmate=false;		
		
		if (ba2.RefID!=ba1.MateRefID) {
			// this set of bam files doesn't have mateRefID.... log this somewhere other than cerr, cout
			cerr << " Mates out of order :   " << ba1.Name <<"\t" << ba1.MateRefID << "\t" << ba2.RefID << "\n";  // oops
			break;
		}
		
		// check read name should be the same but not the same read (depend on First/Second mate flags). 
		if ( ( ba1.Name.compare(ba2.Name)==0 ) && (ba1.IsFirstMate()!=ba2.IsFirstMate()) )  { 
			findmate=true;			
			NmateFound++;
			if (scan) break;
			if (ba2.MapQuality>=Qmin )  break;
			// low quality reads 
			NmateFound=0;
		}										
		
	}	
	
  if(NmateFound<0) {
		return(false);
	}
	
	if (findmate ) {
		
		NmateFound = Bam2MosaikPair(ba1, ba2, ma1,ma2);
		
		// check consistent LF
    if ( (ba1.RefID==ba2.MateRefID) && (ba1.RefID==ba1.MateRefID) && (ma1.IsReverseStrand!=ma2.IsReverseStrand)  ) {		
			int LF1=int(ma2.ReferenceEnd)-int(ma1.ReferenceBegin)-1;
			if (ma1.IsReverseStrand) { 
				LF1=int(ma1.ReferenceEnd)-int(ma2.ReferenceBegin)-1;
			}
			string s=FR;
			if (LF1!=LF)  {				
				if (ma1.ReferenceBegin>ma2.ReferenceBegin) { 
					reverse(s.begin(),s.end());
				}
				cerr	 << " fraglength problem " << LF1 << "\t" << LF << "\t" << s << "\t" << ba1.Name << endl;
			}
		}	
		
	} else {
		NmateFound = Bam2MosaikRead(ba1, ma1);
	}
	
	if (NmateFound>=0)  {
		
		mr1.Name=ba1.Name;
		mr1.ReadGroupCode=ma1.ReadGroupCode;
		
		// shortest read length
		// shortestReadLength=(shortestReadLength<ma1.QueryEnd) ? shortestReadLength: ma1.QueryEnd;
		
		switch (NmateFound) {
				
			case 1:
				mr1.Mate1Alignments.push_back(ma1);
				mr1.Mate2Alignments.push_back(ma2);
				mr1.IsPairedEnd=true;
				break;
				
			case 2:
				mr1.Mate1Alignments.push_back(ma2);
				mr1.Mate2Alignments.push_back(ma1);
				mr1.IsPairedEnd=true;
				break;
				
			case 0:
				mr1.Mate1Alignments.push_back(ma1);
				
		} 
		
		
	}  //else {
	//cerr << " problem finding mate \t" << ba1.Name << " in chrom " << ba1.MateRefID << "\n"; 
	//}
	
	return (true);		
}

//------------------------------------------------------------------------------
// convert a Bam read pair to one complete Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentPairSortedByName( BamMultiReader & ar1, Mosaik::AlignedRead & mr1) 
{
	
	
	//------------------------------------------------------------------------------
	// select abberant pairs that meet Qmin criteria, include dangling (unmapped) ends
	//------------------------------------------------------------------------------
	
	// Mosaik structures
	Mosaik::Alignment ma1,ma2;
	Mosaik::Alignment mz1,mz2;
	
	// Bam structure
	BamAlignment ba1,ba2;
	
	// clear Mosaik structures
	ma1.ReferenceBegin=0;
	ma1.ReferenceEnd=0;
	ma1.ReferenceIndex=0;
	ma1.QueryBegin=0;
	ma1.QueryEnd=0;
	ma1.Quality=0;
	ma1.IsPairedEnd=false;
	ma1.IsReverseStrand=false;
	ma2.ReferenceBegin=0;
	ma2.ReferenceEnd=0;
	ma2.ReferenceIndex=0;
	ma2.QueryBegin=0;
	ma2.QueryEnd=0;
	ma2.Quality=0;
 	ma2.IsPairedEnd=false;
	ma2.IsReverseStrand=false;
	mr1.Mate1Alignments.clear();
	mr1.Mate2Alignments.clear();
	
	
	int LF=0;
	bool findmate=false;	
	int NmateFound=-1;
	string FR="FR";
	bool rev,revMate;
	bool properOrientation = false;	
	
	
	while ( ar1.GetNextAlignment(ba1) ) {
		
		if (doneFrag[ba1.Name]) {
			continue;
		}
		
		doneFrag[ba1.Name]=true;
		
		// skip unmapped reads (another function in a module somewhere ...)
		if (ba1.RefID<0) {
			continue;
		}	
		
		// single ends not processed here
		if (!ba1.IsPaired()) { 
			continue;
		}
		
		if ( (!scan)&&(ba1.MapQuality<Qmin) ) {
			continue;
		}	
		
		// dangling end is processed
		if (!ba1.IsMateMapped() ) {
			NmateFound=1;  // set to mark dangling end
			break;
		}	
		
		//if (ba1.Position>=1130515) { 
		// cerr << endl; }
		
		// skip proper pairs	=  abberant pair criteria 	
		LF = ba1.InsertSize;
		// SLX RP only  ??? 
		if (ba1.IsReverseStrand()) LF=-LF; 
		string rgid;
		ba1.GetReadGroup(rgid);
		unsigned int rgcode=libraries.ReadGroupID2Code[rgid];
		int LMlow = libraries.libmap[rgcode].LMlow;
		int LMhigh = libraries.libmap[rgcode].LMhigh;
		
		// orientation transform to Illumina FR
    
		rev=ba1.IsReverseStrand();
		revMate=ba1.IsMateReverseStrand();
		
		//454 -> IL
		if ( MateMode==MATEMODE_454) {
			if (ba1.IsFirstMate()) {
				rev=!rev;			
			} else {
				revMate=!revMate;
			}
			// SOLiD -> IL
		} else if ( MateMode==MATEMODE_SOLID) {
			if (ba1.IsFirstMate()) {
				revMate=!revMate;
			} else {
				rev=!rev;			
			}
		}		
		
		
		if (ba1.RefID==ba1.MateRefID) {
			if (ba1.Position<ba1.MatePosition) {				
				properOrientation = (!ba1.IsReverseStrand())&&ba1.IsMateReverseStrand();
			}	else {
				properOrientation = ba1.IsReverseStrand()&&(!ba1.IsMateReverseStrand());
			}
		}
		
		//-------------------------------------------------------------------------------
		// skip proper pairs
		//-------------------------------------------------------------------------------
		if (ba1.IsMateMapped()&&(LF>LMlow)&&(LF<LMhigh)&&properOrientation&&(!scan)) {					
			continue;
		}	
		
		// artifact check for mate mapped = this read (illumina problem)
		if  ( (ba1.RefID==ba1.MateRefID) && (ba1.Position==ba1.MatePosition) ) { 
			continue;
		}
		
		//------------------------------------------------------------------------------
		// go mate hunting
		//------------------------------------------------------------------------------		
		findmate=false;		
		
		if ( ar1.GetNextAlignment(ba2) ) {
			
			if (ba2.RefID!=ba1.MateRefID) {
				// this set of bam files doesn't have mateRefID.... log this somewhere other than cerr, cout
				cerr << " Mates out of order :   " << ba1.Name <<"\t" << ba1.MateRefID << "\t" << ba2.RefID << "\n";  // oops
				break;
			}
			
			// check read name should be the same but not the same read (depend on First/Second mate flags). 
			if ( ( ba1.Name.compare(ba2.Name)==0 ) && (ba1.IsFirstMate()!=ba2.IsFirstMate()) )  { 
				findmate=true;			
				NmateFound++;
				break;
			}								
			
		}
		
		
		if (findmate && (scan || (ba2.MapQuality>=Qmin) )) {
			break;
		}	
		
		if (findmate && (!scan) && (ba2.MapQuality<Qmin )) {
			// back up to next first alignment
			findmate=false;
			NmateFound=0;
			continue;
		}	
		
		
	}	
	
  if(NmateFound<0) {
		return(false);
	}
	
	if (findmate ) {
		
		NmateFound = Bam2MosaikPair(ba1, ba2, ma1,ma2);
		
		// check consistent LF
    if ( (ba1.RefID==ba2.MateRefID) & properOrientation) {		
			int LF1=ma2.ReferenceEnd-ma1.ReferenceBegin-1;
			if (ma1.ReferenceBegin>ma2.ReferenceBegin) 
				LF1=ma1.ReferenceEnd-ma2.ReferenceBegin-1;
			if (LF1!=LF)  {
				char i=char(1-ma1.IsReverseStrand*1);
				char s1=FR[i];
				i=char(1-ma2.IsReverseStrand*1);
				char s2=FR[i];
				//s+= FR[char(1-ma2.IsReverseStrand*1)];
				cerr	 << " fraglength problem " << LF1 << "\t" << LF << "\t" << s1 << s2 << "\t" << ba1.Name << endl;
		  }
		}	
		
	} else {
		NmateFound = Bam2MosaikRead(ba1, ma1);
	}
	
	if (NmateFound>=0)  {
		
		mr1.Name=ba1.Name;
		mr1.ReadGroupCode=ma1.ReadGroupCode;
		
		// shortest read length
		// shortestReadLength=(shortestReadLength<ma1.QueryEnd) ? shortestReadLength: ma1.QueryEnd;
		
		switch (NmateFound) {
				
			case 1:
				mr1.Mate1Alignments.push_back(ma1);
				mr1.Mate2Alignments.push_back(ma2);
				mr1.IsPairedEnd=true;
				break;
				
			case 2:
				mr1.Mate1Alignments.push_back(ma2);
				mr1.Mate2Alignments.push_back(ma1);
				mr1.IsPairedEnd=true;
				break;
				
			case 0:
				mr1.Mate1Alignments.push_back(ma1);
				
		} 
		
		
	}  //else {
	//cerr << " problem finding mate \t" << ba1.Name << " in chrom " << ba1.MateRefID << "\n"; 
	//}
	
	return (true);		
}




//------------------------------------------------------------------------------
// convert a Bam read to one complete Mosaik aligned read record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentRead( BamMultiReader & ar1,Mosaik::AlignedRead & mr1, bool checkName=false) 
{
	
	// Mosaik structures
	Mosaik::Alignment ma1;
	
	// Bam structure
	BamAlignment ba1;
	
	// clear Mosaik structures
	ma1.ReferenceBegin=0;
	ma1.ReferenceEnd=0;
	ma1.ReferenceIndex=0;
	ma1.QueryBegin=0;
	ma1.QueryEnd=0;
	ma1.Quality=0;
	ma1.IsPairedEnd=false;
	ma1.IsReverseStrand=false;
	
	int pe1=-1;
	
	int ndone = 0;
	bool ok=false;
	while ( ar1.GetNextAlignment(ba1) ) {
		
		if (!checkName) {
			ok = true;
			break; 
		}
		
		// check if fragment done already
		ndone++;
		if (!doneFrag[ba1.Name]) {
			ok=true;
			break;
		}
		
	}
	
	// check out when GetNextAlignment returns false
	if (!ok) return(ok);
	
	if (ndone>1000) {
		cerr << "problem in nextBamAlignmentPair: " << ndone << "\n";
	}
	
	// mark this one as done
	if(checkName) doneFrag[ba1.Name]=true;
		
	pe1 = Bam2MosaikRead(ba1, ma1);
	
	if (pe1>=0)  {
		
		mr1.Name=ba1.Name;
		mr1.ReadGroupCode=ma1.ReadGroupCode;
		
		mr1.Mate1Alignments.push_back(ma1);
		
		return (true);		
		
	} 
	return (false);		
}

*/

int C_pairedfiles::BamCigarData2Len(vector<CigarOp> CigarData, bool Query) {
	// iterate over CIGAR operations to calculate length 
	// Query true : query length
	// Query false : reference length
	int len=0;
	const int numCigarOps = (const int)CigarData.size();
	for (int i = 0; i < numCigarOps; ++i ) { 
		const CigarOp& op = CigarData.at(i);
		
		// if op is MATCH
		if ( op.Type == 'M' ) {
			len+=op.Length;
		}

		// if op is INS
		if ( (!Query)&&( op.Type == 'D' )) {
			len+=op.Length;
		}

		// if op is DEL
		if ( (Query)&&( op.Type == 'I' )) {
			len+=op.Length;
		}
	}
	return(len);
}


int C_pairedfiles::BamCigarData2mm(vector<CigarOp> CigarData) {
	// iterate over CIGAR operations to calculate number of indel bases
	int indel=0;
	const int numCigarOps = (const int)CigarData.size();
	for (int i = 0; i < numCigarOps; ++i ) { 
		const CigarOp& op = CigarData.at(i);
		
		// if op is INS
		if ( op.Type == 'D' ) {
			indel+=op.Length;
		}
		
		// if op is DEL
		if ( op.Type == 'I' ) {
			indel+=op.Length;
		}
	}
	return(indel);
}

				
/*			
				
//------------------------------------------------------------------------------
// convert one Bam read record to one (incomplete) Mosaik aligned pair record 
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentPair1( BamReader & ar1, Mosaik::AlignedRead & mr1) 
{

    // Mosaik structures
    Mosaik::Alignment ma1,ma2;
    
    // Bam structure
    BamAlignment ba1;

    // clear Mosaik structures
    ma1.ReferenceBegin=0;
    ma1.ReferenceEnd=0;
    ma1.ReferenceIndex=0;
    ma1.QueryBegin=0;
    ma1.QueryEnd=0;
    ma1.Quality=0;
    ma2.ReferenceBegin=0;
    ma2.ReferenceEnd=0;
    ma2.ReferenceIndex=0;
    ma2.QueryBegin=0;
    ma2.QueryEnd=0;
    ma2.Quality=0;
  	mr1.Mate1Alignments.clear();
	mr1.Mate2Alignments.clear();
    
	int pe1=-1;
	
    while ( ar1.GetNextAlignment(ba1) ) {

        // fill first mr1 alignment record
        mr1.Name=ba0.Name;
        // stash alignment info into one end or the other
		
		int optZA=1;
        pe1 = Bam2Mosaik(ba1, ma1,ma2,optZA);
		if (pe1>=0) break;
	}
    
    if (pe1>=0)  {
		mr1.Name=ba1.Name;
		mr1.ReadGroupCode=ma1.ReadGroupCode;

        switch (pe1) {
        case 1:
          mr1.Mate1Alignments.push_back(ma1);
          mr1.Mate2Alignments.push_back(ma2);
          mr1.IsPairedEnd=true;
          break;
        case 2:
          mr1.Mate1Alignments.push_back(ma2);
          mr1.Mate2Alignments.push_back(ma1);
          mr1.IsPairedEnd=true;
          break;
        case 0:
          mr1.Mate1Alignments.push_back(ma1);
        } 
        
		return (true);		
  
    } 
    return (false);		
}
*/


/*

//------------------------------------------------------------------------------
// convert Bam read records to Mosaik aligned pair record 
// This function will fetch the info the the last BAM record loaded, and continue
// to load BAM records until the read name changes
// Name0 is the read name from the previous BAM record
//------------------------------------------------------------------------------
bool  C_pairedfiles::nextBamAlignmentPairSort( BamReader & ar1, Mosaik::AlignedRead & mr1) 
{
    // pair end
    int pe1;

    // alignment counter
		int alignmentCount = 0;
    int pairCount=0;

    // previous read name - 
    BamAlignment ba1;
    ba1=ba0;
    string name0=ba1.Name;
    size_t Lname=name0.size();

    // Mosaik structures
    Mosaik::Alignment ma1,ma2;

    // clear Mosaik structures
    ma1.ReferenceBegin=0;
    ma1.ReferenceEnd=0;
    ma1.ReferenceIndex=0;
    ma1.QueryBegin=0;
    ma1.QueryEnd=0;
    ma1.Quality=0;
    ma2.ReferenceBegin=0;
    ma2.ReferenceEnd=0;
    ma2.ReferenceIndex=0;
    ma2.QueryBegin=0;
    ma2.QueryEnd=0;
    ma2.Quality=0;
  	mr1.Mate1Alignments.clear();
		mr1.Mate2Alignments.clear();

    // if prev read exists (with a name) parse for pe0 and ID0
    if (Lname>0) {
       //if (parseBamReadName(name0, pe1, id0) ) {      
  
       // fill first mr1 alignment record
       mr1.Name=name0;

        // stash alignment info into one end (ma1) and/or the other end (ma2)
       pe1 =Bam2Mosaik(ba1, ma1,ma2,0);
       if (pe1>0) {
          mr1.Mate1Alignments.push_back(ma1);
          pairCount++;
          alignmentCount = 1;
       }      
    }
        
    // loop over reads    
    while ( ar1.GetNextAlignment(ba1) ) {

        ba0=ba1;

        //if (parseBamReadName(ba1.Name, pe1, id1) ){      
  
        if ((ba1.Name!=mr1.Name)&&(Lname>0)) {
           return true;
        }  
        // fill first mr1 alignment record
        mr1.Name=ba1.Name;
        Lname=mr1.Name.size();
        // stash alignment info into one end or the other
        pe1 = Bam2Mosaik(ba1, ma1,ma2,0);
        if (pe1>0) {
          pairCount++;
          
          switch (pairCount) {
            case 1:
              mr1.Mate1Alignments.push_back(ma1);
              mr1.Name=ba1.Name;
              mr1.ReadGroupCode=ma1.ReadGroupCode;
              break;
            case 2:
              mr1.Mate2Alignments.push_back(ma1);
              mr1.Name=ba1.Name;
              mr1.ReadGroupCode=ma1.ReadGroupCode;
              mr1.IsPairedEnd=true;
              break;
            default:
              cerr << "not paired-end data" << endl;
              exit(1);
          }

        }   

				++alignmentCount;

		/ *
			  cout << "----------------------------" << endl;
				cout << "Alignment " << alignmentCount << endl;
				cout << ba1.Name << endl;
				cout << ba1.AlignedBases << endl;
     		cout << "Aligned to " <<  set[0].anchors.names[ba1.RefID] << ":" << ba1.Position << endl;
        * /

    } 
    // no more alignments - reset ba0     
    this->ba0.Name.clear();
    return (alignmentCount>0);		
}
*/


// function to check for existing Set element of vector
char C_pairedfiles::checkSetName(string & setname1)
{
  // loop over existing contigsets
  vector<string>::const_iterator cii;
  int ii;
  for(ii=0; ii < int(this->setNames.size()); ii++)
  {
    // return index to existing contigset
    if (this->setNames[ii]==setname1) {
      return ii;
    }
  }
  // no set found, make a new one
  C_set set1(set1);
  this->setNames.push_back(set1.getSetName());
  this->set.push_back(set1); 
  this->Nset=this->setNames.size();
  return Nset-1;
}

//----------------------------------------------------------------------------
// check Bam file
//----------------------------------------------------------------------------

bool C_pairedfiles::checkBamFile(string & filename1) {

  const char* bamFilename1 = filename1.c_str();

  // check if this is a file...
  DIR  *d;
  d = opendir(bamFilename1);
  if (d) {
	  return false;
  }

  BamReader file1;
	// file1.SetFilename(bamFilename1);

  file1.Open(bamFilename1,"",true);
  if (file1.GetReferenceCount()>0) { 
   	//cerr << endl;
	//	cerr << "BAM" << endl; 
    return true;
	}
  file1.Close();
  return false;
}


bool C_pairedfiles::checkSpannerFiles(string & fs1) {
  //C_headers h;
  vector<string> ff;
  //C_libraries libs;
  int nf = selectfiles(ff,fs1);
  if (nf<7) return false;
  //----------------------------------------------------------------------------
  // load anchor info  
  //----------------------------------------------------------------------------
  //C_anchorinfo a1();
  for (int i=0; i<nf; i++) {
     string f1=ff[i];
     size_t found = f1.find(".anchors.txt");
     if (found!=string::npos) {
        C_anchorinfo a(f1);
        if (anchors.L.size()>0)  {
          if (!(anchors == a)) {
            cerr << "\t mismatched anchors :\t"<< f1 << endl;
            exit(1);
          }
        }
        anchors = a;
     }
     found = f1.find(".library.span");
     if (found!=string::npos) {
        C_libraries libs1(f1);
        if (libraries.libmap.size()==0) {
          libraries=libs1;
        }  else {
          // merge this set of libraries 
          C_librarymap::iterator imap;
          for(imap=libs1.libmap.begin(); imap != libs1.libmap.end(); ++imap) {
            unsigned int ReadGroupCode = imap->first;
            C_libraryinfo lib1=imap->second;
            if (libraries.libmap.count(ReadGroupCode)>0) {
              // problem with merging ReadGroupCodes for this set of libraries 
               cerr << "\t redundant ReadGroupCode:\t"<< ReadGroupCode << endl;
               cerr << "\t from:\t"<< f1 << endl;
               exit(1);
            }
            libraries.libmap[ReadGroupCode]=lib1;
          }
        }
     }
  }
  if (anchors.L.size()<1) return false;
  if (libraries.libmap.size()<1) return false;
   
  C_headerSpan h1;
  vector<string> x;
  split(x, fs1,"/");
  // reconstruct path - everything but last thing after last "/"
  string path="/";
  for (int i=0; i<int(x.size()-1); i++) {
    path = path + x[i]+"/";
  }
  // filename stub
  string stub = x[x.size()-1];
  // span file list
  //vector<string> s1;
  //s1.clear();
  int Nspanx = int(h1.spanext.size());
  for (int i=0; i<Nspanx; i++) {
     string f1 = fs1+h1.spanext[i];
     fstream input(f1.c_str(), ios::in|ios::binary);
     if (!input) {
        cerr << "Unable to open input file: " << f1 << endl;
        headers.clear();
        return false;
     }
     C_headerSpan hf(input);
     // set type here
     hf.type=i;
     headers[f1]=hf;
     input.close();
  }
  return true;  
}

bool C_pairedfiles::checkSpannerDirectory(string & fs1) {
  // ---------------------------------------------------------------------------
  // check all files in directory for *.pair.span
  // ---------------------------------------------------------------------------
  string patternReadSpan("(.+)\\.pair\\.span$");
  string match;
  string f;
  DIR  *d;
  struct dirent *dir;
  SpannerFileNames.clear();
  d = opendir(fs1.c_str());
  if (d)
  {
    while ((dir = readdir(d)) != NULL)
    {
      string filename = dir->d_name;
      if (RE2::FullMatch(filename.c_str(),patternReadSpan.c_str(),&match)) {      
         f = match;
         printf("%s\n", f.c_str());
         SpannerFileNames.push_back(f);
      }
    }
    closedir(d);
  }
  return (SpannerFileNames.size()>0);  
}


bool C_pairedfiles::checkBamDirectory(string & fs1) {
	// ---------------------------------------------------------------------------
	// check all files in directory for *.pair.span
	// ---------------------------------------------------------------------------
	string patternReadBam("(.+)\\.bam$");
	string match;
	string fd,ff;
	DIR  *d;
	struct dirent *dir;
	BamFileNames.clear();
	d = opendir(fs1.c_str());
	if (d)
	{
		// trim trailing "/" from directory name
		fd=fs1;		
		if (fd.find_last_of("/")==fd.size()) {
			fd.resize(fd.size()-1);
		}
		
		while ((dir = readdir(d)) != NULL)
		{
			string filename = dir->d_name;
			if (RE2::FullMatch(filename.c_str(),patternReadBam.c_str(),&match)) {      
				ff = fd+"/"+filename;				
				if (checkBamFile(ff)) {
					printf("%s\n", ff.c_str());
					BamFileNames.push_back(ff);
				}
			}
		}
		closedir(d);
	}
	return (BamFileNames.size()>0);  
}