/*
 *  Histo.cpp
 *  cnv
 *
 *  Created by Chip Stewart on 11/3/07.
 *  Copyright 2007 Boston College. All rights reserved.
 *
 */
#include <iostream>
using std::cin;
using std::cout;
using std::clog;
using std::cerr;
using std::endl;
#include <iomanip>
using std::setw;
using std::setprecision;
#include <math.h>
#include "Function-Generic.h"
#include "Histo.h"

//==========================================================
// Histogram class constructor 
//==========================================================

HistObj::HistObj() {
    Nbin = 0;
    Ntot = 0;
    Nin = 0;
    title = "";
    xlabel = "x";
	collapsed=false;
}

//==========================================================
// constructor with arg list
//==========================================================

HistObj::HistObj(const vector<double> & x, const int Nb, const double Xlow, const double Xhigh) {
    this->Fill(x, Nb, Xlow, Xhigh);
}

//============================================================
// Histogram class constructor from stored file 
//============================================================

HistObj::HistObj(string & filename) {
    ifstream file1;
    file1.open(filename.c_str(), ios::in);
    if (!file1) {
        cerr << "HistObj unable to open file: " << filename << endl;
        exit(121);
    }
    string line;
    // title
    // output <<setw(10)<<"TOT"<< setw(10)<<"MEAN"<< setw(10) <<"STD"
    // <<setw(10)<<"IN"<< setw(10)<<"OVER"<< setw(10) <<"UNDER" << endl;
    getline(file1, line);
    // scan file until "TOT" and "STD" are found in the line...
    while ((line.find("TOT") == string::npos) & (line.find("STD") == string::npos)) {
        title = line;
        getline(file1, line);
    }
    //numbers
    getline(file1, line);
    //output <<setw(10)<<H1.Ntot<< setw(10)<<H1.mean<< setw(10) <<H1.std
    //  <<setw(10)<<H1.Nin<< setw(10)<<H1.Nover<< setw(10) <<H1.Nunder<< endl;
    vector<string> s;
    int i = split(s, line, " \t");
    if (i < 6) {
        cerr << "HistObj unable to parse file header: " << filename << endl;
        exit(122);
    }
    Ntot = string2Double(s[0]);
    mean = string2Double(s[1]);
    std = string2Double(s[2]);
    Nin = string2Double(s[3]);
    Nover = string2Double(s[4]);
    Nunder = string2Double(s[5]);
    getline(file1, line);
    //output <<setw(10)<<"bin"<< setw(10)<<"x"<< setw(10)
    //       <<"n" << setw(10) <<"cum" << endl;
    i = split(s, line, " \t");
    if (i == 4) {
        xlabel = s[1];
    }
    Nbin = 0;
    //vector<double> xc1;
    //vector<int> n1;
    //vector<int> cn1;
    mode = 0;
    mode1 = 0;
    median = -1e20;
    bool donemedian = false;
    double nmax = -1;
    double nmax1 = -1;
    while (getline(file1, line)) {
        if (line.size() == 0) {
            continue;
        }
        i = split(s, line, " \t");
        if (i < 4) {
            cerr << "HistObj unable to parse file : " << filename << endl;
            cerr << "line  : " << line << endl;
            exit(123);
        }
        Nbin++;
        //int i1 = string2Int(s[2]);
        xc.push_back(string2Double(s[1]));
        n.push_back(string2Double(s[2]));
        c.push_back(string2Double(s[3]));
        if (i > 4) {
            binlabels.push_back(s[4]);
        }
        if (n.back() > nmax) {
            mode = xc.back();
            nmax = n.back();
        }
        // mode after first bin (tends to collect junk in first bin)
        if ((Nbin > 1) && (n.back() > nmax1)) {
            mode1 = xc.back();
            nmax1 = n.back();
        }
        if ((c.back() >= (Ntot / 2.0)) && (!donemedian)) {
            median = xc.back();
            donemedian = true;
        }
    }
    file1.close();
};
//========================================================================
// Histogram class destructor 
//========================================================================

HistObj::~HistObj() {
    // cout << title << " boom " << endl;
}

// init with bins only

void HistObj::Initialize(const int Nb, const double Xlow, const double Xhigh) {
    // initialize
    this->Nbin = Nb;
    this->xlow = Xlow;
    this->xhigh = Xhigh;
    this->Ntot = 0;
    this->Nunder = 0;
    this->Nover = 0;
    this->Nin = 0;
    this->sumx = 0;
    this->sumxx = 0;
    this->mode=0;
    this->median=0;
    this->mean=0;
    this->std=0;
    this->normalize=false;
    this->mode1=0;
    this->collapsed=false;
	this->expanded=false;
    this->title="";
    this->xlabel="";    
    // bin width
    this->dx = (Xhigh - Xlow) / Nb;
    this->n.clear();
    this->c.clear();
    this->xc.clear();
    // allocate space in n, c, and xc vectors
    for (int i = 0; i<this->Nbin; i++) {
        this->n.push_back(0);
        this->c.push_back(0);
        this->xc.push_back(Xlow + dx / 2 + dx * i);
    }
}
//========================================================================
// Fill entries into histogram method
//========================================================================

void HistObj::Fill(const vector<double> & x, const int Nb, const double Xlow, const double Xhigh) {
    this->Initialize(Nb, Xlow, Xhigh);
    // fill n, and accumulate sumx, and sumxx
    for (int i = 0; i<this->Ntot; i++) {
        int bin = int(floor((x[i] - Xlow) / dx));
        if (bin < 0) {
            this->Nunder++;
        } else if (bin >= Nb) {
            this->Nover++;
        } else {
            this->Nin++;
            this->n[bin]++;
            this->sumx += x[i];
            this->sumxx += pow(x[i], 2);
        }
    }
    // calc cumulative sum of entries c
    this->Finalize();
};

//========================================================================
// Fill single entry into histogram method
//========================================================================

void HistObj::Fill1(double x1) {
    // initialize
    if (this->Nbin < 1) {
        cerr << "Histogram not initialized for Fill1" << endl;
    }
    this->Ntot += 1;
    // fill n, and accumulate sumx, and sumxx
    int bin = int(floor((x1 - this->xlow) / this->dx));
    if (bin < 0) {
        this->Nunder += 1;
    } else if (bin >= this->Nbin) {
        this->Nover += 1;
    } else {
        this->Nin += 1;
        this->n[bin] += 1;
        this->sumx += x1;
        this->sumxx += pow(x1, 2);
    }
};

void HistObj::Fill1(int x1) {
    this->Fill1(double(x1));
}

void HistObj::Fill1(short x1) {
    this->Fill1(double(x1));
}

//========================================================================
// Fill single entry into histogram method
//========================================================================

void HistObj::FillW(double x1, double w1) {
    // initialize
    if (this->Nbin < 1) {
        cerr << "Histogram not initialized for Fill1" << endl;
    }
    this->Ntot += w1;
    // fill n, and accumulate sumx, and sumxx
    int bin = int(floor((x1 - this->xlow) / this->dx));
    if (bin < 0) {
        this->Nunder += w1;
    } else if (bin >= this->Nbin) {
        this->Nover += w1;
    } else {
        this->Nin += w1;
        this->n[bin] += w1;
        this->sumx += x1*w1;
        this->sumxx += w1*pow(x1, 2);
    }
};

void HistObj::FillW(int x1, double w1) {
    this->FillW(double(x1),w1);
}

void HistObj::FillW(short x1, double w1) {
    this->FillW(double(x1),w1);
}

void HistObj::Finalize() {
    // calc cumulative sum of entries c
    this->c[0] = this->Nunder + this->n[0];
    double nmax = 0;
    double nmax1 = 0;
    this->median = -1e20;
    bool donemedian = false;
    int binHi = 0;
    for (int i = 1; i<this->Nbin; i++) {
        this->c[i] = this->c[i - 1] + this->n[i];
        // check for highest xc bin that contains something
        if ((binHi == 0)&(this->c[i] == this->Ntot)) {
            binHi = i;
        }
        if (this->n[i] > nmax) {
            this->mode = this->xc[i];
            nmax = this->n[i];
        }
        if ((i > 0) && (this->n[i] > nmax1)) {
            this->mode1 = this->xc[i];
            nmax1 = this->n[i];
        }
        if ((this->c[i] >= (this->Ntot / 2.0)) && (!donemedian)) {
            this->median = this->xc[i];
            donemedian = true;
        }
    }
    // trim bins beyond last non-zero
    if (binHi > 0) {
        Nbin = binHi + 1;
        this->n.resize(Nbin);
        this->c.resize(Nbin);
        this->xc.resize(Nbin);
    }
    // mean and stdev
    this->mean = this->sumx / this->Nin;
    this->std = sqrt(this->sumxx / this->Nin - pow(this->mean, 2));
    //return (this->Nin+this->Nover+this->Nunder)>0;
}
//========================================================================
// return p value for a given x value
// p is the fraction of the distribution at or below x (bin containing x)
//========================================================================

double HistObj::x2p(double x) {
    double p = 0;
    // calc bin number
    int bin = int(floor((x - this->xlow) / this->dx));
    if (bin < 0) {
        p = double(this->c[0]) / this->Ntot;
    } else if (bin >= this->Nbin) {
        p = 1;
    } else {
        p = double(this->c[bin]) / this->Ntot;
    }
    return p;
}
//========================================================================
// return x value for a given p value
// p is the fraction of the distribution at or below x (bin containing x)
//========================================================================

double HistObj::p2x(double p) {
    for (int i = 1; i<this->Nbin; i++) {
        double p1 = double(this->c[i]) / this->Ntot;
        if (p1 > p) {
            return this->xc[i - 1];
        }
    }
    return this->xc[this->Nbin - 1];
}
//========================================================================
// return p value for a given x value
// p is the fraction of the Trimmed distribution at or below x (bin containing x)
//========================================================================

double HistObj::x2pTrim(double x) {
    double p = 0;
    double tot = double(this->Ntot - this->Nunder - this->Nover);
    // calc bin number
    int bin = int(floor((x - this->xlow) / this->dx));
    if (bin < 0) {
        p = double(this->c[0] - this->Nunder) / tot;
    } else if (bin >= this->Nbin) {
        p = 1;
    } else {
        p = double(this->c[bin] - this->Nunder) / tot;
    }
    return p;
}
//========================================================================
// return x value for a given p value
// p is the fraction of the Trimmed distribution at or below x (bin containing x)
//========================================================================

double HistObj::p2xTrim(double p) {
    double tot = double(this->Ntot - this->Nunder - this->Nover);
    for (int i = 1; i<this->Nbin; i++) {
        double p1 = double(this->c[i] - this->Nunder) / tot;
        if (p1 > p) {
            return this->xc[i - 1];
        }
    }
    return this->xc[this->Nbin - 1];
}

HistObj HistObj::collapse(int N) {
	// collapse bins from Nbin to N and return new histo
	if (this->Nbin < N) {
		return *this;
	}
	double xmin = 1e10;
	double xmax = -1e10;
	for (int i = 0; i<this->Nbin; i++) {
		if (this->n[i] > 0) {
			if (this->xc[i] < xmin) {
				xmin = this->xc[i];
			}
			if (this->xc[i] > xmax) {
				xmax = this->xc[i];
			}
		}
	}
	dx = (xmax - xmin) / (N - 1);
	// don't collapse this guy any more than 1 bin
	if (dx == 0) {
		return *this;
	}
	HistObj H2;
	H2.Initialize(N, xmin - dx / 2, xmax + dx / 2); // constructor with bins only
	for (int i = 0; i<this->Nbin; i++) {
		if (this->n[i] > 0) {
			int bin = int(floor((this->xc[i] - H2.xlow) / H2.dx));
			H2.n[bin] += this->n[i];
		}
	}
	H2.title= this->title;
	H2.Ntot = this->Ntot;
	H2.Nin = this->Nin;
	H2.mean = this->mean;
	H2.std = this->std;
	H2.Nover = this->Nover;
	H2.Nunder = this->Nunder;
	H2.sumx = this->sumx;
	H2.sumxx = this->sumxx;
	H2.Finalize();
	H2.collapsed=true;
	return H2;
}

HistObj HistObj::expand() {
	// expand bins - fill zero bins
	
	double xmin = 1e10;
	double xmax = -1e10;
	double dx1 = 1e10;
	for (int i = 0; i<this->Nbin; i++) {
		if (i>0) {
			if (dx1>(xc[i]-xc[i-1])) {
				 dx1=xc[i]-xc[i-1];
			}
		}
		if (this->n[i] > 0) {
			if (this->xc[i] < xmin) {
				xmin = this->xc[i];
			}
			if (this->xc[i] > xmax) {
				xmax = this->xc[i];
			}
		}
	}
  int N1 = 1+round((xmax - xmin) / dx1);
	// don't collapse this guy any more than 1 bin
	if (dx == 0) {
		return *this;
	}
	HistObj H2;
	H2.Initialize(N1, xmin - dx1 / 2, xmax + dx1 / 2); // constructor with bins only
	for (int i = 0; i<this->Nbin; i++) {
		if (this->n[i] > 0) {
			int bin = int(floor((this->xc[i] - H2.xlow) / H2.dx));
			H2.n[bin] += this->n[i];
		}
	}
	H2.title= this->title;
	H2.Ntot = this->Ntot;
	H2.Nin = this->Nin;
	H2.mean = this->mean;
	H2.std = this->std;
	H2.Nover = this->Nover;
	H2.Nunder = this->Nunder;
	H2.sumx = this->sumx;
	H2.sumxx = this->sumxx;
	H2.Finalize();
	H2.collapsed=false;
	H2.expanded=true;
	return H2;
}

//========================================================================
// i/o stream operator <<
//========================================================================

ostream & operator<<(ostream &output, const HistObj & H1) {
    output << H1.title << endl;
    output << setw(10) << "TOT" << setw(15) << "MEAN" << setw(15) << "STD"
            << setw(15) << "IN" << setw(15) << "OVER" << setw(15) << "UNDER" << endl;
    output << setw(10) << H1.Ntot << setw(15) << setprecision(8) << H1.mean
            << setw(15) << setprecision(8) << H1.std
            << setw(15) << H1.Nin << setw(15) << H1.Nover << setw(15) << H1.Nunder << endl;
	if (!H1.collapsed) {
		output << setw(10) << "bin" << setw(10) << H1.xlabel << setw(12) << "n" << setw(12) << "cum";
	} else {
		output << setw(10) << "-" << setw(10) << H1.xlabel << setw(12) << "n" << setw(12) << "cum";

	}

    if (int(H1.binlabels.size()) == H1.Nbin) {
        output << setw(12) << "label" << endl;
    } else {
        output << endl;
    }
    for (int i = 0; i < H1.Nbin; i++) {
        if (H1.n[i] > 0) {
			if (!H1.collapsed) {
				output << setw(10) << i << setw(10) << H1.xc[i] << setw(12) << setprecision(8)
                    << H1.n[i] << setw(12) << setprecision(8) << H1.c[i];
			} else {
				output << setw(10) << "-" << setw(10) << H1.xc[i] << setw(12) << setprecision(8)
				<< H1.n[i] << setw(12) << setprecision(8) << H1.c[i];
			}

            if (int(H1.binlabels.size()) >i) {
                output << setw(12) << H1.binlabels[i] << endl;
            } else {
                output << endl;
            }
        }
    }
    return output;
};

void HistObj::setTitle(const string & s1) {
    title = s1;
}

void HistObj::setXlabel(const string & s1) {
    xlabel = s1;
}

void HistObj::setBinLabels(vector<string> & vs) {
    binlabels = vs;
}

HistObj& HistObj::operator=(const HistObj &rhs) {
    Nbin = rhs.Nbin;
    xlow = rhs.xlow;
    xhigh = rhs.xhigh;
    Ntot = rhs.Ntot;
    Nin = rhs.Nin;
    Nover = rhs.Nover;
    Nunder = rhs.Nunder;
    mode = rhs.mode;
    mode1 = rhs.mode1;
    median = rhs.median;
    dx = rhs.dx;
    n = rhs.n;
    xc = rhs.xc;
    c = rhs.c;
    mean = rhs.mean;
    std = rhs.std;
    sumx = rhs.sumx;
    sumxx = rhs.sumxx;
    title = rhs.title;
    xlabel = rhs.xlabel;
    binlabels = rhs.binlabels;
    normalize = rhs.normalize;
	collapsed = rhs.collapsed;
	
    return *this;
}


//========================================================================
// multiple histogram  group class constructor 
//========================================================================

C_HistoGroups::C_HistoGroups() {
 Groups.clear();
}
	
//========================================================================
// multiple histogram  group class constructor from stored file 
//========================================================================

C_HistoGroups::C_HistoGroups(string & filename) {
	ifstream file1;
	file1.open(filename.c_str(), ios::in);
	if (!file1) {
		cerr << "C_Histos unable to open file: " << filename << endl;
		exit(124);
	}
	bool more =true;
	while (more){ 
		C_Histos h1(file1);
		if (h1.h.size()>0) {
			this->Groups.push_back(h1);
			for (int rg=0; rg<h1.ReadGroupTag.size(); rg++) {
				size_t pos = h1.ReadGroupTag[rg].find("ID:");      // position after "ID:" in str
				string str1 = h1.ReadGroupTag[rg].substr(pos+3);   // get from "ID:" to the end
				if (str1.find("\t")!=string::npos) {
					str1.erase( str1.find_first_of( "\t" ) );
				}
				if (str1.find(" ")!=string::npos) {
					str1.erase( str1.find_last_not_of( " " ) + 1 );
					str1.erase( 0, str1.find_first_not_of( " " ) );
				}
				this->ReadGroupIndex[str1]=int(this->Groups.size()-1);
			}
		} else {
			more = false;
		}
	}
}
//========================================================================
// multiple histogram  class defaults constructor 
//========================================================================
C_Histos::C_Histos() {
    ReadGroupTag.clear();
    h.clear();
};


//========================================================================
// multiple histogram class constructor from stored file 
//========================================================================

C_Histos::C_Histos(string & filename) {
    ifstream file1;
    file1.open(filename.c_str(), ios::in);
    if (!file1) {
        cerr << "C_Histos unable to open file: " << filename << endl;
        exit(125);
    }
	C_Histos h1(file1);
	this->h = h1.h;
	this->ReadGroupTag = h1.ReadGroupTag;	
	/*
    string line;
    vector<string> s;
    int i, n;
    // title
    // output <<setw(10)<<"TOT"<< setw(10)<<"MEAN"<< setw(10) <<"STD"
    // <<setw(10)<<"IN"<< setw(10)<<"OVER"<< setw(10) <<"UNDER" << endl;
    bool more = true;
    while (more) {
        more = getline(file1, line);
        if (!more) break;
        n++;
        HistObj h1;
		
		
		// scan file until "@RG" tags are found in the line...
        while ((line.find("@RG") == string::npos) & (line.find("ID") == string::npos)) {
            getline(file1, line);
        }
		// scan file until "@RG" tags are not found in the line...
        while ((line.find("@RG") != string::npos) & (line.find("ID") != string::npos)) {
			ReadGroupTag.push_back(line);
            getline(file1, line);
        }
		break;
	}		
	n = 0;
	char L = '1';
	while (more) {
		more = getline(file1, line);
		if (!more) break;
		n++;
		HistObj h1;
			
		// scan until title line
		while (line.size()<1) {
            getline(file1, line);
        }
		
        i = split(s, line, " \t");
        if (i > 0) h1.title = s[0];
        if (h1.title.size() == 0) {
            h1.title = L;
            L++;
        }
        // scan file until "TOT" and "STD" are found in the line...
        while ((line.find("TOT") == string::npos) & (line.find("STD") == string::npos)) {
            getline(file1, line);
        }
        //numbers
        getline(file1, line);
        //output <<setw(10)<<H1.Ntot<< setw(10)<<H1.mean<< setw(10) <<H1.std
        //  <<setw(10)<<H1.Nin<< setw(10)<<H1.Nover<< setw(10) <<H1.Nunder<< endl;
        i = split(s, line, " \t");
        if (i < 6) {
            cerr << "HistObj unable to parse file header: " << filename << endl;
            exit(1);
        }
        h1.Ntot = string2Double(s[0]);
        h1.mean = string2Double(s[1]);
        h1.std = string2Double(s[2]);
        h1.Nin = string2Double(s[3]);
        h1.Nover = string2Double(s[4]);
        h1.Nunder = string2Double(s[5]);
        //
        getline(file1, line);
        //output <<setw(10)<<"bin"<< setw(10)<<"x"<< setw(10)
        //       <<"n" << setw(10) <<"cum" << endl;
        i = split(s, line, " \t");
        if (i == 4) {
            h1.xlabel = s[1];
        }
        h1.Nbin = 0;
        //vector<double> xc1;
        //vector<int> n1;
        //vector<int> cn1;
        h1.mode = 0;
        h1.mode1 = 0;
        h1.median = -1e20;
        bool donemedian = false;
        double nmax = -1;
        double nmax1 = -1;
        while (getline(file1, line)) {
            if (line.size() == 0) {
                break;
            }
            i = split(s, line, " \t");
            if (i < 4) {
                cerr << "C_Histo unable to parse file : " << filename << endl;
                cerr << "line  : " << line << endl;
                exit(1);
            }
            h1.Nbin++;
            //int i1 = string2Int(s[2]);
            h1.xc.push_back(string2Double(s[1]));
            h1.n.push_back(string2Double(s[2]));
            h1.c.push_back(string2Double(s[3]));
            if (i > 4) {
                h1.binlabels.push_back(s[4]);
            }
            if (h1.n.back() > nmax) {
                h1.mode = h1.xc.back();
                nmax = h1.n.back();
            }
            if ((h1.Nbin > 1) && (h1.n.back() > nmax1)) {
                h1.mode1 = h1.xc.back();
                nmax1 = h1.n.back();
            }
            if ((h1.c.back() >= (h1.Ntot / 2.0)) && (!donemedian)) {
                h1.median = h1.xc.back();
                donemedian = true;
            }
        }
        h[h1.title] = h1;
    }
	*/
    file1.close();
};

C_Histos::C_Histos(ifstream & file1) {
	string line;
	vector<string> s;
	int i, n;
	// title
	// output <<setw(10)<<"TOT"<< setw(10)<<"MEAN"<< setw(10) <<"STD"
	// <<setw(10)<<"IN"<< setw(10)<<"OVER"<< setw(10) <<"UNDER" << endl;
	bool more = true;
	n = 0;
	while (more) {
		more = getline(file1, line);
		if (!more) break;
		n++;
		HistObj h1;
		
		
		// scan file until "@RG" tags are found in the line...
		//while ((line.find("RG") == string::npos) & (line.find("ID") == string::npos)) {
		while (line.find("RG") == string::npos)  {
			more=getline(file1, line);
			if (!more) exit(126);
				
		}
		// scan file until "@RG" tags are not found in the line...
		//while ((line.find("RG") != string::npos) & (line.find("ID") != string::npos)) {
		while (line.find("RG") != string::npos) {
			if  (line.find("ID") == string::npos) {
				cerr << " missing ID: tag in stat file RG header" << endl;
				cerr << line<< endl;
				exit(119);
			}
			ReadGroupTag.push_back(line);
			more=getline(file1, line);
			if (!more) exit(127);
		}		
		break;
	}	
	
	// check for optional Mapping quality threshold line
	if (line.find("Mapping quality threshold") != string::npos)  {
		more=getline(file1, line);
		if (!more) exit(120);
	}		
	
	//
	
	n = 0;
	char L = '1';
	while (more) {
		
		
		streampos fp=file1.tellg();
		
		more = getline(file1, line);
		if (!more) break;
		
		if  ((line.find("@RG") != string::npos) & (line.find("ID") != string::npos)) {
			
			file1.seekg(fp);
			break;
		}
		
		n++;
		HistObj h1;
		
		// scan until title line
		while (line.size()<1) {
			getline(file1, line);
		}
		
		i = split(s, line, " \t");
		if (i > 0) h1.title = s[0];
		if (h1.title.size() == 0) {
			h1.title = L;
			L++;
		}
		// scan file until "TOT" and "STD" are found in the line...
		while ((line.find("TOT") == string::npos) & (line.find("STD") == string::npos)) {
			getline(file1, line);
		}
		//numbers
		getline(file1, line);
		//output <<setw(10)<<H1.Ntot<< setw(10)<<H1.mean<< setw(10) <<H1.std
		//  <<setw(10)<<H1.Nin<< setw(10)<<H1.Nover<< setw(10) <<H1.Nunder<< endl;
		i = split(s, line, " \t");
		if (i < 6) {
			cerr << "HistObj unable to parse file header: " << endl;
			exit(128);
		}
		h1.Ntot = string2Double(s[0]);
		h1.mean = string2Double(s[1]);
		h1.std = string2Double(s[2]);
		h1.Nin = string2Double(s[3]);
		h1.Nover = string2Double(s[4]);
		h1.Nunder = string2Double(s[5]);
		//
		getline(file1, line);
		//output <<setw(10)<<"bin"<< setw(10)<<"x"<< setw(10)
		//       <<"n" << setw(10) <<"cum" << endl;
		i = split(s, line, " \t");
		if (i == 4) {
			h1.xlabel = s[1];
		}
		h1.Nbin = 0;
		//vector<double> xc1;
		//vector<int> n1;
		//vector<int> cn1;
		h1.mode = 0;
		h1.mode1 = 0;
		h1.median = -1e20;
		bool donemedian = false;
		double nmax = -1;
		double nmax1 = -1;
		while (getline(file1, line)) {
			if (line.size() == 0) {
				break;
			}
			i = split(s, line, " \t");
			if (i < 4) {
				cerr << "C_Histo unable to parse file : " << endl;
				cerr << "line  : " << line << endl;
				exit(129);
			}
			h1.Nbin++;
			//int i1 = string2Int(s[2]);
			h1.xc.push_back(string2Double(s[1]));
			h1.n.push_back(string2Double(s[2]));
			h1.c.push_back(string2Double(s[3]));
			if (i > 4) {
				h1.binlabels.push_back(s[4]);
			}
			if (h1.n.back() > nmax) {
				h1.mode = h1.xc.back();
				nmax = h1.n.back();
			}
			if ((h1.Nbin > 1) && (h1.n.back() > nmax1)) {
				h1.mode1 = h1.xc.back();
				nmax1 = h1.n.back();
			}
			if ((h1.c.back() >= (h1.Ntot / 2.0)) && (!donemedian)) {
				h1.median = h1.xc.back();
				donemedian = true;
			}
		}
		HistObj h2 = h1.expand();
		h[h1.title] = h2;
	}
};

/*
 ostream & operator<<(ostream &output, const C_Histos & Hs) {
    map<string, HistObj,less<string> >::iterator ih; // = Hs.h.begin();   
    for ( ih=Hs.h.begin() ; ih != Hs.h.end(); ih++ ) {
           output << (*ih).second << endl;
    }
}
*/
    
//========================================================================
// 1D descriptive statistics class constructor
//========================================================================

StatObj::StatObj() {
    N = 0;
}

StatObj::StatObj(const vector<double> & x, const int Nb, const double Xlow, const double Xhigh) {
    this->Fill(x, Nb, Xlow, Xhigh);
}

//========================================================================
// initialize Stat Object
//========================================================================

void StatObj::Initialize(const int Nb, const double Xlow, const double Xhigh) {
    this->N = 0;
    this->sumx = 0;
    this->sumxx = 0;
    this->h.n.clear();
    this->h.c.clear();
    this->h.xc.clear();
    this->h.Initialize(Nb, Xlow, Xhigh);
}
//========================================================================
// Fill 1 event
//========================================================================

void StatObj::Fill1(double x1) {
    this->N++;
    this->sumx += x1;
    this->sumxx += pow(x1, 2);
    this->h.Fill1(x1);
}
// for int input

void StatObj::Fill1(int x1) {
    this->Fill1(double(x1));
}
// for short input

void StatObj::Fill1(short x1) {
    this->Fill1(double(x1));
}
//========================================================================
// Finalize
//========================================================================

void StatObj::Finalize() {
    this->mean = this->sumx / this->N;
    this->std = sqrt(this->sumxx / this->N - pow(this->mean, 2));
    this->h.Finalize();
}
//========================================================================
// calc stats 
//========================================================================

void StatObj::Fill(const vector<double> & x, const int Nb, const double Xlow, const double Xhigh) {
    this->Initialize(Nb, Xlow, Xhigh);
    this->N = int(x.size());
    for (int i = 0; i<this->N; i++) {
        this->sumx += x[i];
        this->sumxx += pow(x[i], 2);
        this->h.Fill1(x[i]);
    }
    this->Finalize();
}

void StatObj::Fill(const vector<int> & x, const int Nb, const double Xlow, const double Xhigh) {
    this->Initialize(Nb, Xlow, Xhigh);
    this->N = int(x.size());
    for (int i = 0; i<this->N; i++) {
        this->sumx += x[i];
        this->sumxx += pow(x[i], 2);
        this->h.Fill1(x[i]);
    }
    this->Finalize();
}

void StatObj::Fill(const vector<short> & x, const int Nb, const double Xlow, const double Xhigh) {
    this->Initialize(Nb, Xlow, Xhigh);
    this->N = int(x.size());
    for (int i = 0; i<this->N; i++) {
        this->sumx += x[i];
        this->sumxx += pow(x[i], 2);
        this->h.Fill1(x[i]);
    }
    this->Finalize();
}

double StatObj::threshold(const double x) {
    return this->h.p2xTrim(x);
}

//========================================================================
// i/o stream operator <<
//========================================================================

ostream & operator<<(ostream &output, const StatObj & S1) {
    if (S1.h.title.size() > 0) {
        output << S1.h.title << endl;
        //if (S1.h.xlabel.size() > 0) output << S1.h.xlabel;
        //output << endl;
    }
    output << setw(10) << "TOT" << setw(15) << "MEAN" << setw(15) << "STD" << endl;
    output << setw(10) << S1.N << setw(15) << S1.mean << setw(15) << S1.std << endl;
    return output;
}
