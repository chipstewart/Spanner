/*
 *  headerInfo.cpp
 *  SpanDet
 *
 *  Created by Chip Stewart on 8/19/08.
 *  Copyright 2008 Boston College. All rights reserved.
 *
 */

#include "headerSpan.h"

C_headerSpan::C_headerSpan() {               // constructor
    V=207; 
    contigName="";                           
    setName="";                              
    typeName="";
    reclen=0;
    light=0;
    type = 0;
    //--------------------------------------------------------------------------
    // Spanner build file set extensions
    //--------------------------------------------------------------------------
    string s[] = {".pair.span",".cross.span",".dangle.span",".multi.span"};
    vector<string> s1(s, s + 4);
    spanext = s1;
    //--------------------------------------------------------------------------
    // Spanner optional intermediate (clustering...) file set extensions
    //--------------------------------------------------------------------------
    string sc[] = {".short.cls.span",".long.cls.span",".invert5.cls.span",
    ".invert3.cls.span",".dangle5.cls.span",".dangle3.cls.span",
    ".cross5.cls.span",".cross3.cls.span"};
    vector<string> s1c(sc, sc + 8);
    spanextc = s1c;
    /*
    vector<string> s1e(s + 6, s + 7);
    spanext1 = s1e;
    */
}

ostream &operator<<(ostream &output,  C_headerSpan & h)
{
    output << "header\t  Version:" << h.V << "\t " << h.typeName << endl ;      
    output << " set:\t " << h.setName <<  "\t contig:\t " << h.contigName ;      
    if (h.light!=0) {
      output << " Light :\t " << h.light << endl ;      
    } else {
      output << endl ;      
    }
    output << " Number :\t " << h.N <<  "\t reclen:\t " << h.reclen << endl ;      
    return output;
}

//------------------------------------------------------------------------------
// Constructor from file
//------------------------------------------------------------------------------
C_headerSpan::C_headerSpan(fstream & input) 
{
  input.seekg (0, ios::beg);
  if (!input.is_open()) { 
      cerr << "file not open in loadHeader" << endl;
  }
  // grab version number from latest default constructor
  C_headerSpan h;
  int V0 = h.V;
  // buffer for loading strings
  char buff[512];
  // read version 
  input.read(reinterpret_cast < char * > (&V), sizeof(int));
  if (V!=V0) {
      cerr << "Spanner file version "<< V << " doesn't match expected version " << V0 << endl;
  }  
  // length of c-style string
  int nc = 0;
  input.read(reinterpret_cast < char * > (&nc), sizeof(int));
  //  contigName
  input.read(buff, nc);
  // convert c-style string buff to string contigName 
  buff[nc]=0;
  contigName = buff;
  // length of setName string
  int ns = 0;
  input.read(reinterpret_cast < char * > (&ns), sizeof(int));
  //  contigName
  input.read(buff, ns);
  // convert c-style string buff to string contigName 
  buff[ns]=0;
  setName = buff;
  // length of c-style string
  int nt = 0;
  input.read(reinterpret_cast < char * > (&nt), sizeof(int));
  //  contigName
  input.read(buff, nt);
  // convert c-style string buff to string contigName 
  buff[nt]=0;
  typeName = buff;
  // read record size
  input.read(reinterpret_cast < char * > (&reclen), sizeof(int));
  // read info long long light if version > 204
  if (h.V>204) { 
    // read info 
    input.read(reinterpret_cast < char * > (&light), sizeof(light));
  }
  // read number of fragments 
  input.read(reinterpret_cast < char * > (&N), sizeof(int));
}


// I/O function write
bool C_headerSpan::write(fstream & output) 
{
  // open output binary file. bomb if unable to open
  //fstream output(outfilename.c_str(), ios::out | ios::binary);
  if (!output.is_open()) {
      cerr << "Unable to write header into output file " << endl;
      return false;
  }
  // write version 
  output.write(reinterpret_cast<const char *>(&V), sizeof(int));
  // convert contigName to c-style string
  const char * cC = this->contigName.c_str();    
  // length of c-style string
  int cCL = strlen(cC)+1;
  // write length of contigName 
  output.write(reinterpret_cast<const char *>(&cCL), sizeof(int));
  // write contigName
  output.write((const char *)(cC), cCL * sizeof(char));
  // convert setName to c-style string
  const char * sC = this->setName.c_str();    
  // length of c-style string
  int sCL = strlen(sC)+1;
  // write length of setName 
  output.write(reinterpret_cast<const char *>(&sCL), sizeof(int));
  // write setName
  output.write((const char *)(sC), sCL * sizeof(char));
  // convert typeName to c-style string
  const char * sT = this->typeName.c_str();    
  // length of c-style string
  int sTL = strlen(sT)+1;
  // write length of setName 
  output.write(reinterpret_cast<const char *>(&sTL), sizeof(int));
  // write setName
  output.write((const char *)(sT), sTL * sizeof(char));
  // write record size
  output.write(reinterpret_cast<const char *>(&reclen), sizeof(int));
  // write info
  output.write(reinterpret_cast<const char *>(&light), sizeof(light));
  // write number of records
  output.write(reinterpret_cast<const char *>(&N), sizeof(int));
  return true;
}

