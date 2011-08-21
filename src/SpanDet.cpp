//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// SpanDet.cpp
// Chip Stewart
// Copyright 2008 Boston College. All rights reserved.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include "SpanDet.h"

C_SpannerSV::C_SpannerSV(C_pairedfiles & data,  RunControlParameters & pars) {
    C_contigs::const_iterator iterContig;   
    int Nset = data.set.size();
    
    for (int iset = 0; iset< Nset; iset++) {
        string setName = data.set[iset].getSetName();
        
        // x18Mar2011 int Qmin = data.set[iset].pars.getQmin();
        
        libraries = data.set[iset].libraries;
        //readFractions=libraries.readFractionSamples();
        for (iterContig = data.set[iset].contig.begin();	iterContig != data.set[iset].contig.end(); iterContig++) {
            string name = iterContig->first;
            
            string area = pars.getOutputDir();
            int k = 1+area.find_last_of("/");
            string subdir = area.substr(k);
            if (area.size()>0) { area = area+"/";}
            string prefix = pars.getPrefix();
            if (prefix.size()>0) { prefix = prefix+".";}
            string sn1 = setName;
            // dont need redundant setName if subdirectory is already the setName
            // if (sn1==subdir) sn1="";
            if (sn1.size()>0) { sn1 = sn1+".";}
            string basename = area+prefix+sn1+name;
            if (setName==name) basename = area+prefix+name;
            
            // replace evil character "|" with benign "_" 
            size_t found=basename.find("|");
            while (found!=string::npos) {
                basename.replace(found,1,"_");
                found=basename.find("|");
            }
            
            string fname;
            
            //-----------------------------------------------------------------------
            // bug out when missing pair info
            //-----------------------------------------------------------------------
            if (data.set[iset].contig[name].pairStats.N==0) { continue ; }
            
            //-----------------------------------------------------------------------
            // Clustering
            //-----------------------------------------------------------------------
            C_SpannerCluster clus(data.set[iset].contig[name], libraries, pars);
            clus.writeall();
            
			/*
             //-----------------------------------------------------------------------
             // calc read starts for CNV if not already done...
             //-----------------------------------------------------------------------
             if (data.set[iset].contig[name].read_start.n.size()==0) {
             // check for existing .read.start.span file in input area...
             fname = area+prefix+sn1+name+".read.start.span";
             ifstream fin;
             fin.open(fname.c_str());
             if( !fin ) {
             data.set[iset].contig[name].calcStarts(Qmin);
             fname = basename+".read.start.span";
             cout << "write\t " << fname << endl;
             data.set[iset].contig[name].writeDepth(fname,
             data.set[iset].contig[name].read_start);  
             } else { 
             cout << "load\t " << fname << endl;
             data.set[iset].contig[name].read_start=data.set[iset].contig[name].loadDepth(fname);
             data.set[iset].contig[name].read_start.calcStats();
             }         
             }
             
             //-----------------------------------------------------------------------
             // nominal coverage distributions for CNV p-values
             //-----------------------------------------------------------------------
             int nua[8] = {50,100,200,500,1000,2000,5000,10000};
             vector<int> nu0(nua,nua+8);
             C_NominalCov nomcov1(data.set[iset],  nu0); 
             nomcov = nomcov1;
             
             //-----------------------------------------------------------------------
             // calc fragment depth for retro insertions if not already done...
             //-----------------------------------------------------------------------
             if (data.set[iset].contig[name].frag_depth.n.size()==0) {
             fname = area+prefix+sn1+name+".frag.depth.span";
             ifstream fin;
             fin.open(fname.c_str());
             if( !fin ) {
             // change this to use library LMlow LMhigh
             int lmLow = int(pars.getFragmentLengthLo());
             int lmHigh = int(pars.getFragmentLengthHi());
             int L =  data.set[iset].contig[name].Length;
             data.set[iset].contig[name].calcFragDepth(0,L,lmLow,lmHigh);
             fname = basename+".frag.depth.span";
             cout << "write\t  " << fname << endl;
             data.set[iset].contig[name].writeDepth(fname,
             data.set[iset].contig[name].frag_depth);  
             } else { 
             cout << "load \t " << fname << endl;
             data.set[iset].contig[name].frag_depth = data.set[iset].contig[name].loadDepth(fname);
             data.set[iset].contig[name].frag_depth.calcStats();
             }       
             
             }
			 */
			
            //-----------------------------------------------------------------------
            // check for mobile element insertions
            //-----------------------------------------------------------------------
            vector<string> patternElement = pars.getMobileElements();
			if (patternElement.size()<1) {
				patternElement.push_back("moblist_");
			}
			C_retroElements retrolist(data.set[iset].contig[name].anchors, patternElement[0]);
            int NE=retrolist.e.size();
            
            for (int j=0; j<NE; j++) {
                int e = retrolist.e[j];
                string retroType = retrolist.name[j];
                if (retroType[2]=='.') { retroType = retroType.substr(0,2);}
                
                //-----------------------------------------------------------------------
                // Retro mobile element Clustering
                //-----------------------------------------------------------------------
                
                C_SpannerRetroCluster rclus(data.set[iset].contig[name],libraries,pars,e,retroType);
                
                //-----------------------------------------------------------------------
                // mobile element masking
                //-----------------------------------------------------------------------
                string maskfile1 = pars.getMobiMaskFile();
                size_t ia=maskfile1.find("$");
                if (ia!=string::npos) { 
                    string maskfile=maskfile1;
                    maskfile.replace(ia, 1, retroType);
                    if (!FileExists(maskfile)) {
                        maskfile=maskfile1;
                        string rtype=upperCase(retroType);
                        maskfile.replace(ia, 1, rtype);
                        //cout << " uppercase annotations  " << maskfile << endl;
                    }
                    if (!FileExists(maskfile)) {
                        maskfile=maskfile1;
                        string rtype=lowerCase(retroType);
                        maskfile.replace(ia, 1, rtype);
                        //cout << " lowercase annotations  " << maskfile << endl;
                    }
                    if (!FileExists(maskfile)) {
                        maskfile=maskfile1;
                        string rtype=Titlecase(retroType);
                        maskfile.replace(ia, 1, rtype);
                        //cout << " Titlecase annotations  " << maskfile << endl;
                    } 
                    
                    
                    if (FileExists(maskfile)) {
                        C_BedChr mask(maskfile,name);
                        mask.name = retroType;
                        ret.Mask=mask;
                        
                        cout << " loaded " << mask.p.size() ;
                        cout << " annotations from " << maskfile << endl;
                        
                        int LMmax = int(libraries.maxLF()/2.0);
                        // mask out clusters near annotated elements
                        rclus.Mask(mask,LMmax);
                        
                    }
                }
                
                //-----------------------------------------------------------------------
                // write Element Insertion clusters
                //-----------------------------------------------------------------------
                rclus.writeall();
                
                //-----------------------------------------------------------------------
                // ELement Insertions
                //-----------------------------------------------------------------------
                int Nret = findRet(data.set[iset].contig[name], rclus, pars);
                ret.contigName=name;
                ret.setName=setName;
                cout << " find " << Nret << " " << retroType << " element insertions" << endl;
                fname = basename+"."+retroType+".svcf";
                cout << "write " << retroType << " element insertions : " << fname << endl;
                ret.print(fname);
				C_SVR ret1;
				ret=ret1;
            }
            
            
            //-----------------------------------------------------------------------
            // Cross links
            //-----------------------------------------------------------------------
            int Ncrx = findCross(data.set[iset].contig[name], clus, pars);
            crx.typeName="crosslinks";
            crx.contigName=name;
            crx.setName=setName;
            cout << " find " << Ncrx << " cross chromosomal links" << endl;
            fname = basename+".cross.svcf";
            cout << "write cross chromosomal links : " << fname << endl;
            crx.print(fname);
            
            //-----------------------------------------------------------------------
            // Inversions
            //-----------------------------------------------------------------------
            int Ninv = findInv(data.set[iset].contig[name], clus, pars);
            inv.typeName="inversions";
            inv.contigName=name;
            inv.setName=setName;
            cout << " find " << Ninv << " inversions" << endl;
            fname = basename+".inv.svcf";
            cout << "write inversions : " << fname << endl;
            inv.print(fname);
            
            
            //-----------------------------------------------------------------------
            // Deletions
            //-----------------------------------------------------------------------
            int Ndel = findDel(data.set[iset].contig[name], clus, pars);
            del.typeName="deletions";
            del.contigName=name;
            del.setName=setName;
            cout << " find " << Ndel << " deletions" << endl;
            //int Ndup = findDel(data.set[iset].contig[name], clus, pars);
            //del.typeName="deletions";
            fname = basename+".del.svcf";
            cout << "write deletions : " << fname << endl;
            del.print(fname);      
            
            //-----------------------------------------------------------------------
            // Duplications
            //-----------------------------------------------------------------------
            int Ndup = findDup(data.set[iset].contig[name], clus, pars);
            dup.typeName="duplications";
            dup.contigName=name;
            dup.setName=setName;
            cout << " find " << Ndup << " duplications" << endl;
            fname = basename+".dup.svcf";
            cout << "write duplications : " << fname << endl;
            dup.print(fname); 
            
            
            //-----------------------------------------------------------------------
            // Done - print summary 
            //-----------------------------------------------------------------------
            fname = basename+".summary.span.txt";
            cout << "write summary : " << fname << endl;                    
            printSummary(fname, data.set[iset].contig[name]);
        }
    }
}


void C_SpannerSV::printSummary(string & outfilename, C_contig & c1) {
    // format: optimized to loadFragments.m matlab script  
    // open output binary file. bomb if unable to open
    fstream output(outfilename.c_str(), ios::out);
    if (!output) {
        cerr << "Unable to open summary file: " << outfilename << endl;
        return;
    }
    output << del.typeName << "\t " <<del.setName << "\t "<<del.contigName << endl;
    output << del.evt.size() << " events " << endl;
    output << dup.typeName << "\t " <<dup.setName << "\t "<<dup.contigName << endl;
    output << dup.evt.size() << " events " << endl;
    output << inv.typeName << "\t " <<inv.setName << "\t "<<inv.contigName << endl;
    output << inv.evt.size() << " events " << endl;
    output << "Set " << c1.setName << "\tContig " << c1.getContigName() << endl;
    output << "repeats " << endl;
    output <<   c1.repeat.Stats << endl;
    output <<   c1.repeat.Stats.h << endl;
    output << "pairs " << endl;
    output <<   c1.pairStats << endl;
    output <<   c1.pairStats.h << endl;
    output << "cross " << endl;
    output <<   c1.crossStats << endl;
    output <<   c1.crossStats.h << endl;
    output << "coverage histograms " << endl;
    output <<   nomcov;
    output.close();
}


int C_SpannerSV::findDel(C_contig  & contig, C_SpannerCluster & clus,  
                         RunControlParameters & pars) {
    del.typeName="deletions";
    //int Nc = clus.longpairc.cluster.size();
    C_cluster2d_elements::iterator it;
    C_cluster2d_element1 c1;
    vector<C_localpair> cpair;
    // vector of fraglength shift
	vector<int> DiffLF;
    // vector of mapping qualities 5' end
	vector<int> mapQ5;
    // vector of mapping qualities 3' end
	vector<int> mapQ3;
 	
    //----------------------------------------------------------------------------
    // samples name vector sorted & unique in ret
    //----------------------------------------------------------------------------
    C_librarymap::iterator il;
    vector<string>::iterator is;
    del.samples.clear();
    for ( il=libraries.libmap.begin() ; il != libraries.libmap.end(); il++ )
    {
        C_libraryinfo lib1 = (*il).second; 
        del.samples.push_back(lib1.Info.SampleName);
    }  
    sort(del.samples.begin(),del.samples.end());
    is = unique (del.samples.begin(), del.samples.end()); 
    del.samples.resize( is - del.samples.begin() );  
    
    //----------------------------------------------------------------------------
    // SVCF Info (list here and in C_SVR << )
    //----------------------------------------------------------------------------
    C_SVCF_TAG tag1;
	tag1.INFO=true;
	tag1.Id="SVLEN";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Difference in length between REF and ALT alleles";
    del.SVCF.INFO.push_back(tag1);
    tag1.Id="CIPOS";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Confidence interval around POS";
    del.SVCF.INFO.push_back(tag1);
    tag1.Id="END";
	tag1.Number=1;
	tag1.Descr="Deletion end coordinate";
    del.SVCF.INFO.push_back(tag1);
    tag1.Id="CIEND";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Confidence interval around END";
    del.SVCF.INFO.push_back(tag1);
    tag1.Id="CISVLEN";
	tag1.Descr="Confidence interval around SVLEN";
    del.SVCF.INFO.push_back(tag1);
    tag1.Id="NALTF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of ALT supporting fragments";
    del.SVCF.INFO.push_back(tag1);
    tag1.Id="QLF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Aberrant fragment length metric";
    del.SVCF.INFO.push_back(tag1);
    tag1.Id="QOUT";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Outlier fragment metric";
    del.SVCF.INFO.push_back(tag1);
	tag1.Id="PRF";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Range of F cluster reads";
    del.SVCF.INFO.push_back(tag1);
	tag1.Id="PRR";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Range of R cluster reads";
    del.SVCF.INFO.push_back(tag1);
    tag1.Id="MQ5";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Median mapping quality at 5' end";
    del.SVCF.INFO.push_back(tag1);
    tag1.Id="MQ3";
	tag1.Descr="Median mapping quality at 3' end";
    del.SVCF.INFO.push_back(tag1);
    tag1.INFO=false;
	tag1.FORMAT=true;
	tag1.Id="NF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of ALT supporting fragments/sample";
    del.SVCF.FMT.push_back(tag1);
    tag1.FORMAT=false;
	tag1.ALT=true;
	tag1.Id="DEL";
	tag1.Descr="Deletion relative to REF";
    del.SVCF.ALT.push_back(tag1);
	
	//string s1[] = {"L","NF","NP","UP","UL","PA","PB","PC","PD","AL","NR","ER","MR"};
    //vector<string> Info1(s1, s1 + 13);
    //del.SVCF.Info=Info1;
    
    // SVCF Format  (list here and in C_SVR << operator method)
    // string s2[] = {"NF","NP","N3","N5","NR","ER"}; //"CN"};
    //string s2[] = {"NF","NP","NR","ER"}; //"CN"};
    //string s2[] = {"NF"}; //"CN"};
    //vector<string> Format1(s2, s2 + 1);
    //del.SVCF.Format=Format1;
    
    // SVCF samples
    del.SVCF.Samples=del.samples;    
    
    // SVCF EventType
    del.SVCF.EventType=del.typeName;
    
    
    
    //----------------------------------------------------------------------------
    // use library fragment length properties
    //----------------------------------------------------------------------------
    // use library LF, LFmax, and LFsig for selection rather than pars
    //double LF = pars.getFragmentLength();
    //double LFmax = pars.getFragmentLengthHi();
    //HistObj LFhist = pars.getFragHist();
    // double LFsig= LFhist.std;
    //----------------------------------------------------------------------------
    // expected number of pairs that span a given deletion
    //----------------------------------------------------------------------------
    //double Nexp = contig.localpairs.size()*double(LF-2*contig.aLR)/double(contig.Length);
    //double Nexp = contig.frag_depth.Stats.h.mean;
    // contig alignability estimate
    //float totSites = contig.Length-int(contig.totalNoCovBases);
    //float acontig  = float (totSites-contig.totalRepeatBases) / totSites;
    //----------------------------------------------------------------------------
    // loop over long pair clusters to identify candidate deletions
    //----------------------------------------------------------------------------
    for ( it=clus.longpairc.cluster.begin() ; it != clus.longpairc.cluster.end(); it++ ) {
        int i = (*it).first;
        c1 = clus.longpairc.cluster[i];   
        int Np = c1.inp.size();
        //-------------------------------------------------------------------------
        // initialize cluster limits:    p5a--p5b      p3a-p3b
        //-------------------------------------------------------------------------
        int p5b = 0;              // trailing edge of 5' cluster
        int p5a = contig.Length;  // leading edge of 5' cluster
        int p3b = 0;              // trailing edge of 3' cluster
        int p3a = contig.Length;  // leading edge of 3' cluster
        int p5aMax=0;             // leading edge of trailing read in 5' cluster 
        int p3aMax=0;             // leading edge of trailing read in 3' cluster
        int LFmax=0;              // longest library in this cluster
        unsigned int RGCmax=0;    // ReadGroupCode for longest library
        int LFrange=0;            // LF width of longest library
        double qAberr=0;          // product of fragment consistency with LF
        C_RGmap RGmap;            // declare empty Readgroup counter
        
        // loop over fragments in cluster
        cpair.resize(Np);
        DiffLF.resize(Np);
        mapQ5.resize(Np);
        mapQ3.resize(Np);
        for (int j=0; j<Np; j++) {
            int k = c1.inp[j];
            C_localpair lp1 = clus.longpair[k];
            cpair[j]=lp1;
            int p = lp1.pos;                // start of 5' cluster
            if (p<p5a) p5a=p;
            if (p>p5aMax) p5aMax=p;
            p = lp1.pos+lp1.len1;           // end of 5' cluster
            if (p>p5b) p5b=p;
            p = lp1.pos+lp1.lm-lp1.len2;    // start of 3' cluster
            if (p<p3a) p3a=p;
            if (p>p3aMax) p3aMax=p;
            p = lp1.pos+lp1.lm;             // end of 3' cluster
            if (p>p3b) p3b=p;
            unsigned int RGC1=lp1.ReadGroupCode;
            RGmap[RGC1]++;
            int LMhigh1=libraries.libmap[RGC1].LMhigh;
            if (LMhigh1>LFmax) {
                RGCmax=RGC1;
                LFmax=LMhigh1;
                LFrange=LMhigh1-libraries.libmap[RGC1].LMlow;
            }
            double pAb1=libraries.libmap[RGC1].fragHist.x2pTrim(double(lp1.lm));
            qAberr+=double(p2q(1-pAb1))/Np;
            // frag length shift
            DiffLF[j]=lp1.lm-int(libraries.libmap[RGC1].fragHist.median);
            // mapQ5
            mapQ5[j]=lp1.q1;
            // mapQ3
            mapQ3[j]=lp1.q2;
        }
        
		vector <int> DiffLF1 = DiffLF;
		sort(DiffLF1.begin(),DiffLF1.end());
		int DLF1=DiffLF1[Np/2];
		
		//if ( 2*(Np/2)==Np)  {
        if ( Np % 2 == 0)  {                
			DLF1= int((DiffLF1[Np/2]+DiffLF1[(Np/2)-1])/2.0);
		}
        
		sort(mapQ5.begin(),mapQ5.end());
		int mQ5=mapQ5[Np/2];		
		sort(mapQ3.begin(),mapQ3.end());
		int mQ3=mapQ3[Np/2];		
		if ( 2*(Np/2)==Np)  {
			mQ5= int((mapQ5[Np/2]+mapQ5[-1+(Np/2)])/2.0);
			mQ3= int((mapQ3[Np/2]+mapQ3[-1+(Np/2)])/2.0);
		}
        
        // check for very strange read lengths or indexing screwups
        //if ((p5b+1)!=p0) {
        //if ( abs((p5b+1)-p0)>(int(20*contig.sLR)) ) {
        //   cerr << "problem with cluster index edge calc " << p0 << " ~ " << p5b << endl;
        //}
        //-------------------------------------------------------------------------
        // mapping length ~  extent of deletion 
        // for some reason the average deletion is ~2 reads too big 
        //-------------------------------------------------------------------------
        // int len1 = int(c1.mean[1]-LF-2*contig.aLR);
        //-------------------------------------------------------------------------
        // length <~ distance from the start of 3' cluster and end of 5' cluster 
        //------------------------------------------------------------------------- 
        // int len = int(p3a-p5b-1);
        //-------------------------------------------------------------------------
        // estimate gap size where breakpoint should be 
        //-------------------------------------------------------------------------
        // int avegap = ((p5aMax-p5a)+(p3aMax-p3a))/(2*Np);
        // average gap over this contig
        //int avegapC = double(contig.Length)/int(2*contig.localpairs.size());
        //int avegapC = double(contig.Length)/double(contig.totalUniqueReads);
        
        //-------------------------------------------------------------------------
        // correct p0 and len for expected gaps between reads 
        //-------------------------------------------------------------------------
        // int p0 = p5b+1+avegap/2;
        // int len = int(p3a-p5b-1)-avegap/2;
        int p0 = p5b+1;
		int p1 = p3b;
        int len = int(DLF1);
        /*
         don't extrapolate without threshold info
         //-------------------------------------------------------------------------
         // extrapolate cluster effective lm below threshold
         //-------------------------------------------------------------------------
         int Low = int(c1.low[1]);      // lowest lm in cluster
         int dL = int(c1.high[1]-Low);  // range of lm in cluster
         double EsigLm = LFhist.std;   // expected lm stdev
         double sigLm = c1.std[1]; // observed lm stdev
         int dLx = int(dL*((EsigLm/sigLm) - 1.0));
         bool extend = ((Low-dLx)< int(LFmax));
         if (extend) {
         len = len -dLx;
         }
         */
        //-------------------------------------------------------------------------
        // bail on candidate with negative length...
        //-------------------------------------------------------------------------
        if (len<1) continue; 
        //-------------------------------------------------------------------------
        // check for overlap with previous evt
        //-------------------------------------------------------------------------
        bool overlap = false;
        C_SV1 e0;
        if (del.evt.size()>0) { 
            e0 = del.evt.back();
            int pp = e0.pend;
            overlap=(pp>p0);
            overlap = overlap&&(abs(int(e0.pos)-int(p0))<LFmax);
            
            /*
             if (overlap) {
             cout << " overlap: prev pos " << e0.pos << "\t prev end " 
             << e0.pos+e0.length << "\t this pos " << p0 << endl;
             }
             */
        }
        //cout << p0 << " del len1 " << len1 << " \t len " << len << "\t avegap " << avegap << "\t Np " << Np << endl;
        //len = len1;
        //
        // coverage within deletion region
        //int p1 = p0+len-1;
        //C_SVcoverage1 cov(contig, p0, p1,nomcov) ;
        // coverage within 5' cluster region
        //C_SVcoverage1 cov5(contig, p5a, p5b,nomcov) ;
        // coverage within 3' cluster region
        //C_SVcoverage1 cov3(contig, p3a, p3b,nomcov) ;
        
        C_SV1 e1;
        e1.pos = p0;
        e1.pend = p1;
        e1.anchor = contig.getAnchorIndex();
        e1.length = len;
        //e1.q = p2q(cov.p);
		// adhoc satuating functions
		e1.q=int(((Np*100.0)/(Np+10.0)));
		e1.q5=mQ5;
		e1.q3=mQ3;
        //e1.cov=cov;
        //e1.cov5=cov5;
        //e1.cov3=cov3;
        e1.cls = c1;
        e1.pair=cpair;
        e1.p5[0]=p5a;
        e1.p5[1]=p5b;
        e1.p3[0]=p3a;
        e1.p3[1]=p3b;
        
        
        // brkpoint uncertainty depends on average gap between read starts in cluster
        //e1.posU=(avegap>avegapC? avegap: avegapC);
        e1.CIpos[0]=-(p5aMax-p5a)/Np;
		e1.CIpos[1]=(p3aMax-p3a)/Np;
        e1.CIend[0]=-(p5aMax-p5a)/Np;
		e1.CIend[1]=(p3aMax-p3a)/Np;
		e1.CIlen[0]=-(DLF1-DiffLF1[0])/Np;
		e1.CIlen[1]=(DiffLF1[DiffLF1.size()-1]-DLF1)/Np;
        
        // correction to pos,len for deletion bracketed by repeat region
        //int DL=int(len-int(floor(0.5+e1.cls.mean[1])));
        int DL=len-droundi(e1.cls.mean[1]);
        
        if (DL>(4*int(e1.CIpos[1]))) {
            e1.CIlen[0]=droundi(-e1.cls.std[1]/sqrt(e1.cls.N));
            e1.CIlen[1]=droundi(e1.cls.std[1]/sqrt(e1.cls.N));
            e1.length=droundi(e1.cls.mean[1]);
   			e1.CIpos[0]=int(-DL/2);
	  		e1.CIpos[1]=int(DL/2);
            e1.pos=e1.pos+int(DL/2);
        }
        
        
        // prob(outlier) in cluster is x2p of fragment dist with range of cluster
        double aRange=((p3b-p3a)+(p5b-p5a))/2.0;
        double pOut=libraries.libmap[RGCmax].fragHist.x2pTrim(aRange);
        e1.qOutlier=p2q(1.0-pOut);
        
        // prob(Aberrant LM) for independent fragments
        e1.qAberrantLM=int(qAberr);
        
        // ReadGroups in event
        e1.ReadGroupMap=RGmap;
        
        // end alignability
        //e1.a5 = cov5.Nsite/(float(int(cov5.p1)-int(cov5.p0))*acontig);
        //e1.a3 = cov3.Nsite/(float(int(cov3.p1)-int(cov3.p0))*acontig);
        
        //-------------------------------------------------------------------------
        // if overlapped merge this cluster with prev event
        //-------------------------------------------------------------------------
        if (overlap) {
            C_SV1 e2 = merge(e0,e1,contig, pars, 0);
            // remove existing event at end of list
            del.evt.pop_back();
            e1=e2;
        } 
        
        int Npairs = e1.pair.size();
        // bkp piles up at the corner of Np-Lm space close to 0,LF
        bool significant = double(Npairs) > 0; //(0.01*Nexp - (e1.cls.mean[1])/LFrange);
        // simple cut on number of pairs in event
        significant = significant & (Npairs>=pars.getMinClustered());    
        // simple cut on event length
        significant = significant & (int(e1.length)>=pars.getMinLength());    
        // quantized copy number
        //e1.copynumber = short(2*e1.cov.N/e1.cov.eN +0.5);
        //-------------------------------------------------------------------------
        // new deletion event
        //-------------------------------------------------------------------------
        e1.type="DEL";
        if (significant) { //&(e1.copynumber<(2+pars.getCNslosh()))) { 
            del.evt.push_back(e1);
        }
    }
    
    // sort list
    del.evt.sort();
    
    // calc sample supporting and spanning fragments
    del.finalize(contig,libraries,pars);
    
    // SVCF event count
    del.SVCF.NEvent=int(del.evt.size());  
    
    return int(del.evt.size());
}       


int C_SpannerSV::findDup(C_contig  & contig, C_SpannerCluster & clus,  
                         RunControlParameters & pars) {
    dup.typeName="duplications";
    C_cluster2d_elements::iterator it;
    C_cluster2d_element1 c1;
    vector<C_localpair> cpair;
	// vector of fraglength shift
	vector<int> DiffLF;
	// vector of mapping qualities 5' end
	vector<int> mapQ5;
    // vector of mapping qualities 3' end
	vector<int> mapQ3;
 	
    //----------------------------------------------------------------------------
    // samples name vector sorted & unique in ret
    //----------------------------------------------------------------------------
    C_librarymap::iterator il;
    vector<string>::iterator is;
    dup.samples.clear();
    for ( il=libraries.libmap.begin() ; il != libraries.libmap.end(); il++ )
    {
        C_libraryinfo lib1 = (*il).second; 
        dup.samples.push_back(lib1.Info.SampleName);
    }  
    sort(dup.samples.begin(),dup.samples.end());
    is = unique (dup.samples.begin(), dup.samples.end()); 
    dup.samples.resize( is - dup.samples.begin() );  
    
    //----------------------------------------------------------------------------
    // SVCF Info (list here and in C_SVR << )
    //----------------------------------------------------------------------------
    /*
     string s1[] = {"L","NF","NP","UP","UL","PA","PB","PC","PD","AL","NR","ER","MR"};
     vector<string> Info1(s1, s1 + 13);
     dup.SVCF.Info=Info1;
     
     // SVCF Format  (list here and in C_SVR << operator method)
     string s2[] = {"NF"}; //,"NP","NR","ER"}; //"CN"};
     vector<string> Format1(s2, s2 + 1);
     dup.SVCF.Format=Format1;
     */
	C_SVCF_TAG tag1;
	tag1.INFO=true;
	tag1.Id="SVLEN";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Difference in length between REF and ALT alleles";
    dup.SVCF.INFO.push_back(tag1);
    tag1.Id="CIPOS";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Confidence interval around POS";
    dup.SVCF.INFO.push_back(tag1);
	tag1.Id="END";
	tag1.Number=1;
	tag1.Descr="Duplication end coordinate";
    dup.SVCF.INFO.push_back(tag1);
    tag1.Id="CIEND";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Confidence interval around END";
    dup.SVCF.INFO.push_back(tag1);
	tag1.Id="CISVLEN";
	tag1.Descr="Confidence interval around SVLEN";
    dup.SVCF.INFO.push_back(tag1);
	
    tag1.Id="NALTF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of ALT supporting fragments";
    dup.SVCF.INFO.push_back(tag1);
    tag1.Id="QLF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Aberrant fragment length metric";
    dup.SVCF.INFO.push_back(tag1);
    tag1.Id="QOUT";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Outlier fragment metric";
    dup.SVCF.INFO.push_back(tag1);
	tag1.Id="PRF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Range of F cluster reads";
    dup.SVCF.INFO.push_back(tag1);
	tag1.Id="PRR";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Range of R cluster reads";
    dup.SVCF.INFO.push_back(tag1);
	
	tag1.Id="MQ5";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Median mapping quality at 5' end";
    del.SVCF.INFO.push_back(tag1);
    tag1.Id="MQ3";
	tag1.Descr="Median mapping quality at 3' end";
    del.SVCF.INFO.push_back(tag1);
	
    tag1.INFO=false;
	tag1.FORMAT=true;
	tag1.Id="NF";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Number of ALT supporting fragments/sample";
    dup.SVCF.FMT.push_back(tag1);
    tag1.FORMAT=false;
	tag1.ALT=true;
	tag1.Id="DUP:TANDEM";
	tag1.Descr="Tandem Duplication relative to REF";
    dup.SVCF.ALT.push_back(tag1);
	tag1.Id="INS";
	tag1.Descr="Insertion relative to REF";
    dup.SVCF.ALT.push_back(tag1);
	
    // SVCF samples
    dup.SVCF.Samples=dup.samples;    
    
    // SVCF EventType
    dup.SVCF.EventType=dup.typeName;
    
    
    
    //----------------------------------------------------------------------------
    // fragment length properties
    //----------------------------------------------------------------------------
    // double LF = pars.getFragmentLength();
    // double LFmin = pars.getFragmentLengthLo();
    // HistObj LFhist = pars.getFragHist();
    // double LFsig= LFhist.std;
    //----------------------------------------------------------------------------
    // expected number of pairs that span a given deletion
    //----------------------------------------------------------------------------
    // double Nexp = contig.localpairs.size()*double(LF-2*contig.aLR)/double(contig.Length);
    // global fragments of all libraries
    // double Nexp = contig.frag_depth.Stats.h.mean;
    // float totSites = contig.Length-int(contig.totalNoCovBases);
    // float acontig  = float (totSites-contig.totalRepeatBases) / totSites;
    //----------------------------------------------------------------------------
    // loop over long pair clusters to identify candidate duplications
    //----------------------------------------------------------------------------
    for ( it=clus.shortpairc.cluster.begin() ; it != clus.shortpairc.cluster.end(); it++ ) {
        int i = (*it).first;
        c1 = clus.shortpairc.cluster[i];   
        //-------------------------------------------------------------------------
        // highest start position mark start of duplication 
        //-------------------------------------------------------------------------
        //int p0 = int(c1.high[0]+contig.aLR+1);
        int Np = c1.inp.size();
        //-------------------------------------------------------------------------
        // initialize cluster limits:    p5a--p5b      p3a-p3b
        //-------------------------------------------------------------------------
        int p5b = 0;              // trailing edge of 5' cluster
        int p5a = contig.Length;  // leading edge of 5' cluster
        int p3b = 0;              // trailing edge of 3' cluster
        int p3a = contig.Length;  // leading edge of 3' cluster
        //
        int p5aMax=0;             // leading edge of trailing read in 5' cluster
        int p3aMax=0;             // leading edge of trailing read in 3' cluster
        int LFmax=0;              // longest library in this cluster
        unsigned int RGCmax=0;             // ReadGroupCode for longest library
        int LFrange=0;            // LF width of longest library
        double qAberr=0;          // product of fragment consistency with LF
        C_RGmap RGmap;            // declare empty Readgroup counter
        
        //
        cpair.resize(Np);
		DiffLF.resize(Np);
		mapQ5.resize(Np);
        mapQ3.resize(Np);
        for (int j=0; j<Np; j++) {
            int k = c1.inp[j];
            C_localpair lp1 = clus.shortpair[k];
            cpair[j]=lp1;
            int p = lp1.pos;                // start of 5' cluster
            if (p<p5a) p5a=p;
            if (p>p5aMax) p5aMax=p;
            p = lp1.pos+lp1.len1;           // end of 5' cluster
            if (p>p5b) p5b=p;
            p = lp1.pos+lp1.lm-lp1.len2;    // start of 3' cluster
            if (p<p3a) p3a=p;
            if (p>p3aMax) p3aMax=p;
            p = lp1.pos+lp1.lm;             // end of 3' cluster
            if (p>p3b) p3b=p;
            
            unsigned int RGC1=lp1.ReadGroupCode;
            RGmap[RGC1]++;
            int LMhigh1=libraries.libmap[RGC1].LMhigh;
            if (LMhigh1>LFmax) {
                RGCmax=RGC1;
                LFmax=LMhigh1;
                LFrange=LMhigh1-libraries.libmap[RGC1].LMlow;
            }
            double pAb1=libraries.libmap[RGC1].fragHist.x2pTrim(double(lp1.lm));
            qAberr+=double(p2q(pAb1))/Np;
            
            DiffLF[j]=droundi(libraries.libmap[RGC1].fragHist.median)-lp1.lm;
			
			// mapQ5
			mapQ5[j]=lp1.q1;
			// mapQ3
			mapQ3[j]=lp1.q2;
			
        }
        
		vector <int> DiffLF1 = DiffLF;
		sort(DiffLF1.begin(),DiffLF1.end());
		int DLF1=DiffLF1[int(Np/2)];
		
		if ( 2*(Np/2)==Np)  {
			DLF1= droundi((DiffLF1[Np/2]+DiffLF1[(Np/2)-1])/2.0);
		}
        
        
		sort(mapQ5.begin(),mapQ5.end());
		int mQ5=mapQ5[Np/2];		
		sort(mapQ3.begin(),mapQ3.end());
		int mQ3=mapQ3[Np/2];		
		if ( 2*(Np/2)==Np)  {
			mQ5= droundi((mapQ5[Np/2]+mapQ5[-1+(Np/2)])/2.0);
			mQ3= droundi((mapQ3[Np/2]+mapQ3[-1+(Np/2)])/2.0);
		}
		
		
		//int avegap = ((p5aMax-p5a)+(p3aMax-p3a))/(2*Np);
        // average gap over this contig
        //int avegapC = double(contig.Length)/int(2*contig.localpairs.size());
        //int avegapC = double(contig.Length)/double(contig.totalUniqueReads);
        //-------------------------------------------------------------------------
        // correct p0 and len for expected gaps between reads
        //-------------------------------------------------------------------------
        // int p0 = p3a+1-avegap/2;
        // int len = int(p5b-p3a-1)+avegap/2;
        int p0 = p3a+1;
		int p1 = p5b;
        int len = DLF1;
        
        /*
         int p0 = p3a+1;
         //-------------------------------------------------------------------------
         // length <~ distance from the start of 3' cluster and end of 5' cluster 
         //------------------------------------------------------------------------- 
         int len = int(p5b-p3a-1);
         //-------------------------------------------------------------------------
         // estimate gap size where breakpoint should be 
         //-------------------------------------------------------------------------
         int avegap = ((p5b-p5a)+(p3b-p3a))/(2*Np);
         //-------------------------------------------------------------------------
         // correct p0 and len for expected gaps between reads 
         //-------------------------------------------------------------------------
         p0 = p3a+1-avegap/2;
         len = len+avegap/2;   
         
         //-------------------------------------------------------------------------
         // extrapolate cluster effective lm below threshold
         //-------------------------------------------------------------------------
         int LmHi = int(c1.high[1]);      // highest lm in cluster
         int dL = LmHi-int(c1.low[1]);    // range of lm in cluster
         double EsigLm = LFhist.std;      // expected lm stdev
         double sigLm = c1.std[1]; // observed lm stdev
         int Lmx = LmHi + int(dL*( (EsigLm/sigLm) - 1.0));
         bool extend = (Lmx> int(LFmin));
         if (extend) {
         len = len - (LmHi-Lmx);
         }    
         */
        //-------------------------------------------------------------------------
        // don't bail out for negative length events
        //-------------------------------------------------------------------------
        // if (len<1) continue; 
        //-------------------------------------------------------------------------
        // check for overlap with previous evt
        //-------------------------------------------------------------------------
        bool overlap = false;
        C_SV1 e0;
        if (dup.evt.size()>0) { 
            e0 = dup.evt.back();
            int pp = e0.pend;
            overlap=pp>p0;
            overlap = overlap&&(abs(int(e0.pos)-int(p0))<LFmax);
            /*
             if (overlap) {
             cout << " overlap: prev pos " << e0.pos << "\t prev end " 
             << e0.pos+e0.length << "\t this pos " << p0 << endl;
             }
             */
        }
        //-------------------------------------------------------------------------
        // estimate end of CNV region 
        //-------------------------------------------------------------------------  
        //int p1 = p0+len-1;
        //-------------------------------------------------------------------------
        // calc read coverage between p0 and p1
        //-------------------------------------------------------------------------  
        // C_SVcoverage1 cov(contig, p0, p1,nomcov) ;
        // coverage within 5' cluster region
        // C_SVcoverage1 cov5(contig, p5a, p5b,nomcov) ;
        // coverage within 3' cluster region
        // C_SVcoverage1 cov3(contig, p3a, p3b,nomcov) ;
        //-------------------------------------------------------------------------
        // candidate SV duplication event 
        //-------------------------------------------------------------------------  
        C_SV1 e1;
        e1.pos = p0;
		e1.pend= p1;
        e1.anchor = contig.getAnchorIndex();
        e1.length = len;
        // prob that null CNV exceeds this coverage is 1-cov.p
        //cov.p=fabs(1.0-cov.p);
        e1.q = droundi((100.*Np)/(Np+5)); //p2q(cov.p);
		e1.q5=mQ5;
		e1.q3=mQ3;
        //    
        //e1.cov=cov;  
        e1.cls = c1;
        e1.pair=cpair;
        e1.p5[0]=p5a;
        e1.p5[1]=p5b;
        e1.p3[0]=p3a;
        e1.p3[1]=p3b;
        
        
        // brkpoint uncertainty depends on average gap between read starts in cluster
        //e1.posU=(avegap>avegapC? avegap: avegapC);
		e1.CIpos[0]=-(p5aMax-p5a)/Np;
		e1.CIpos[1]=(p3aMax-p3a)/Np;
        e1.CIend[0]=-(p5aMax-p5a)/Np;
		e1.CIend[1]=(p3aMax-p3a)/Np;
		e1.CIlen[0]=-(DLF1-DiffLF1[0])/Np;
		e1.CIlen[1]=(DiffLF1[DiffLF1.size()-1]-DLF1)/Np;
        
        // correction to pos,len for duplication bracketed by repeat region
        
        int DL=(len+droundi(e1.cls.mean[1]));
        if (DL>int(4*e1.CIpos[1])) {
            e1.CIlen[0]=droundi(e1.cls.std[1]/sqrt(e1.cls.N));
            e1.CIlen[1]=droundi(e1.cls.std[1]/sqrt(e1.cls.N));
            e1.length=droundi(-e1.cls.mean[1]);
  			e1.CIpos[0]=-DL/2;
	  		e1.CIpos[1]=DL/2;
            e1.pos=e1.pos+DL/2;
        }
        
        
        // prob(outlier) in cluster is x2p of fragment dist with range of cluster
        double aRange=((p3b-p3a)+(p5b-p5a))/2.0;
        double pOut=libraries.libmap[RGCmax].fragHist.x2pTrim(aRange);
        e1.qOutlier=p2q(1.0-pOut);
        
        // prob(Aberrant LM) for independent fragments
        e1.qAberrantLM=int(qAberr);
        
        // end alignability
        //e1.a5 = cov5.Nsite/(float(int(cov5.p1)-int(cov5.p0))*acontig);
        //e1.a3 = cov3.Nsite/(float(int(cov3.p1)-int(cov3.p0))*acontig);
        
        // ReadGroups in event
        e1.ReadGroupMap=RGmap;
        
        //-------------------------------------------------------------------------
        // if overlapped merge this cluster with prev event
        //-------------------------------------------------------------------------
        if (overlap) {
            C_SV1 e2 = merge(e0,e1,contig, pars, 1);
            // remove existing event at end of list
            dup.evt.pop_back();
            e1=e2;
        } 
        int Npairs = e1.pair.size();
        
        // bkp piles up at the corner of Np-Lm space close to 0,LF
        //bool significant = double(Npairs) > (0.01*Nexp + (e1.cls.mean[1])/LFrange);
        bool significant = double(Npairs) > 0; //(0.01*Nexp + (e1.cls.mean[1])/LFrange);
        // simple cut on number of pairs in event
        significant = significant & (Npairs>=pars.getMinClustered());
        // simple cut on event length
        significant = significant & (abs(int(e1.length))>=pars.getMinLength());
        // quantized copy number
        // e1.copynumber = short(2*e1.cov.N/e1.cov.eN +0.5);
        //-------------------------------------------------------------------------
        // new duplication event
        //-------------------------------------------------------------------------
        e1.type="TDUP";
        if (e1.length>0) { 
			e1.type="INS";
        }
		//if (significant&(e1.copynumber>(2-pars.getCNslosh()))) {
        if (significant) { 
            dup.evt.push_back(e1);
        }
    }
    
    // sort list
    dup.evt.sort();
    
    // calc sample supporting and spanning fragments
    dup.finalize(contig,libraries,pars);
    
    // SVCF event count
    dup.SVCF.NEvent=int(dup.evt.size());  
    
    
    return int(dup.evt.size());
}       

int C_SpannerSV::findInv(C_contig  & contig, C_SpannerCluster & clus,  RunControlParameters & pars) {
    C_SVV inv3 = findInvDir(contig, clus, pars,3);
    C_SVV inv5 = findInvDir(contig, clus, pars,5);
    inv.typeName="inversions";
    
    //----------------------------------------------------------------------------
    // samples name vector sorted & unique in ret
    //----------------------------------------------------------------------------
    C_librarymap::iterator it;
    vector<string>::iterator is;
    inv.samples.clear();
    for ( it=libraries.libmap.begin() ; it != libraries.libmap.end(); it++ )
    {
        C_libraryinfo lib1 = (*it).second; 
        inv.samples.push_back(lib1.Info.SampleName);
    }  
    sort(inv.samples.begin(),inv.samples.end());
    is = unique (inv.samples.begin(), inv.samples.end()); 
    inv.samples.resize( is - inv.samples.begin() );  
    
    //----------------------------------------------------------------------------
    // SVCF Info (list here and in C_SVR << )
    //----------------------------------------------------------------------------
	/*
     string s1[] = {"L","NF","NP","UP","UL","PA","PB","PC","PD","AL","NR","ER","MR"};
     vector<string> Info1(s1, s1 + 13);
     inv.SVCF.Info=Info1;
     
     // SVCF Format  (list here and in C_SVR << operator method)
     string s2[] = {"NF","NP","N3","N5","NR","ER"}; //"CN"};
     vector<string> Format1(s2, s2 + 6);
     inv.SVCF.Format=Format1;
     */
	C_SVCF_TAG tag1;
	tag1.ALT=true;
	tag1.Id="INV";
	tag1.Descr="Inversion relative to REF";
    inv.SVCF.ALT.push_back(tag1);	
	tag1.ALT=false;
	tag1.INFO=true;
	tag1.Id="SVLEN";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Difference in length between REF and ALT alleles";
    inv.SVCF.INFO.push_back(tag1);
    tag1.Id="CIPOS";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Confidence interval around POS";
    inv.SVCF.INFO.push_back(tag1);
    tag1.Id="NFF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of F ALT supporting fragments";
    inv.SVCF.INFO.push_back(tag1);
    tag1.Id="NFR";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of R ALT supporting fragments";
    inv.SVCF.INFO.push_back(tag1);
    tag1.Id="QC";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Outlier fragment metric";
    inv.SVCF.INFO.push_back(tag1);
	tag1.Id="PC1";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Range of F cluster reads";
    inv.SVCF.INFO.push_back(tag1);
	tag1.Id="PC2";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Range of R cluster reads";
    inv.SVCF.INFO.push_back(tag1);
    tag1.INFO=false;
	tag1.FORMAT=true;
	tag1.Id="NFF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of F ALT supporting fragments";
    inv.SVCF.FMT.push_back(tag1);
	tag1.Id="NFR";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of R ALT supporting fragments";
    inv.SVCF.FMT.push_back(tag1);
  	
    // SVCF samples
    inv.SVCF.Samples=inv.samples;    
    
    // SVCF EventType
    inv.SVCF.EventType=inv.typeName;
    
    // x18Mar2011  double Nexp = contig.frag_depth.Stats.h.mean/2.;
    
	//----------------------------------------------------------------------------
    // expected number of pairs that span a given inversion
    //----------------------------------------------------------------------------
    //double LF = pars.getFragmentLength();
    //double Nexp = contig.localpairs.size()*double(LF-2*contig.aLR)/double(contig.Length);
    //int N3 = inv3.evt.size();
    //int N5 = inv5.evt.size();  
    int Nmerge = 0;
    list<C_SVV1>::iterator i3,i5,i51;
    C_SVV1 e3,e5,e1;
    
    /*
     int N3=inv3.evt.size();
     int N5=inv5.evt.size();
     cout << "inv5\t" << N5 << "\t inv3\t" << N3 << endl;
     */
    
    //----------------------------------------------------------------------------
    // loop over both views (5' and 3') to merge matched inversion events
    //----------------------------------------------------------------------------
    for ( i3=inv3.evt.begin() ; i3 != inv3.evt.end(); i3++ ) {
        e3 = *i3;
        double f1=0;
        for ( i5=inv5.evt.begin() ; i5 != inv5.evt.end(); i5++ ) {
            e5 = *i5;
            if (e5.p3[1]<e3.p5[0]) {
                continue;
            }
            if (e5.p5[0]>e3.p3[1]) break;
            //-------------------------------------------------------------------------
            // check for overlap 
            //-------------------------------------------------------------------------
            float p3low=float(e3.pos)-float(e3.posU);
            float p3high=float(e3.pos)+float(e3.length)+float(e3.posU);
            float p5low=float(e5.pos)-float(e5.posU);
            float p5high=float(e5.pos)+float(e5.length)+float(e5.posU);
            
            bool overlap=! ( (p5high<p3low) || (p5low>p3high) ); // not (before or after)
            //-------------------------------------------------------------------------
            // check for matching lengths - mutual fractional overlap (1 merge, 0 not)
            //-------------------------------------------------------------------------
            if (overlap) {
                int p0=(e5.pos<e3.pos? e5.pos: e3.pos);
                int p1=((e5.pos+e5.length)>(e3.pos+e3.length)? e5.pos+e5.length: e3.pos+e5.length);
                int L = p1-p0;
                //double f=double(e3.length)*double(e5.length)/(L*(double(e3.length)+double(e5.length)));
                double f=double(2*L)/(double(e3.length)+double(e5.length));
                if (f>f1) {
                    f1=f;
                    i51=i5;
                }
            }
        }
        e1=e3;
        if (f1>0.05) {    // some mutual overlap threshold
            e5 = *i51;
            e1 = merge(e3,e5,contig);     
            // remove existing event in inv5 list
            inv5.evt.erase(i51); 
            Nmerge++;
        } 
        int Npairs = e1.pair3.size()+e1.pair5.size();
        // bkp piles up at the corner of Np-Lm space close to 0,LF
        // bool significant = double(Npairs) > (0.01*Nexp); //  + (e1.cls.mean[1]-LF)/LFsig);
        bool significant = double(Npairs) > 0;
        // simple cut on number of pairs in event
        significant = significant & (Npairs>=pars.getMinClustered());    
        // simple cut on event length
        significant = significant & (int(e1.length)>=pars.getMinLength());    
        // quantized copy number
        //e1.copynumber = short(2*e1.cov.N/e1.cov.eN +0.5);
        
        //  length here or over there 
        float L3 = double(e3.length)+0.1;
        float L5 = double(e5.length)+0.1;
        // reicprocal difference with length over there (allow inverted)
        float rL = (fabs(L5)>fabs(L3)? fabs(L5-L3)/fabs(L5): fabs(L5-L3)/fabs(L3));    
        // q value is set by supporting fragments (max 40)
        double qLr  = (e1.length>=0 ? rL/(0.2+rL): 0.1);    
        // q value is set by supporting fragments (max 40)
        double qNF = double(Npairs)/(5.0+double(Npairs));    
        // q value from length here
        double qL  = (e1.length>=0 ? double(e1.length)/(100.0+fabs(double(e1.length))): 0.1);   
        // q value combined
        double q  = sqrt(qNF*qL*qLr);
        
        // q value is set by supporting fragments (max 40)
        e1.q = int(40.0*q/(0.1+q)); 
        
        //-------------------------------------------------------------------------
        // new inversion event
        //-------------------------------------------------------------------------
        //if ((double(e1.cov.N)>(0.1*e1.cov.eN))&significant) { 
        if (significant) { 
            inv.evt.push_back(e1);
        }
    }
    //-------------------------------------------------------------------------
    // remaining inv5 events
    //-------------------------------------------------------------------------
    for ( i5=inv5.evt.begin() ; i5 != inv5.evt.end(); i5++ ) {
        e1 = *i5;
        int Npairs = e1.pair5.size();
        // bkp piles up at the corner of Np-Lm space close to 0,LF
        // bool significant = double(Npairs) > (0.01*Nexp); //  + (e1.cls5.mean[1]-LF)/LFsig); 
        bool significant = double(Npairs) > 0;
        // simple cut on number of pairs in event
        significant = significant & (Npairs>=pars.getMinClustered());    
        // simple cut on event length
        significant = significant & (int(e1.length)>=pars.getMinLength());    
        // quantized copy number
        //e1.copynumber = short(2*e1.cov.N/e1.cov.eN +0.5);
        
        //-------------------------------------------------------------------------
        // new inversion event
        //-------------------------------------------------------------------------
        //if ((double(e1.cov.N)>(0.1*e1.cov.eN))&significant) { 
        if (significant) { 
            inv.evt.push_back(e1);
        }
    }
    
    // sort events
    inv.evt.sort();
    
    // calc sample supporting and spanning fragments
    inv.finalize(contig,libraries,pars);
    
    // SVCF event count
    inv.SVCF.NEvent=int(inv.evt.size());  
    
    return int(inv.evt.size());
}       

C_SVV C_SpannerSV::findInvDir(C_contig  & contig, C_SpannerCluster & clus,  RunControlParameters & pars, int dir) {
    C_SVV inv1;
    if (!( (dir==3)||(dir==5) )) return inv1;
    C_cluster2d_elements::iterator it;
    C_cluster2d_elements::iterator it0;
    C_cluster2d_elements::iterator it1;
    C_cluster2d_element1 c1;
    vector<C_localpair> cpair;
    //----------------------------------------------------------------------------
    // fragment length properties
    //----------------------------------------------------------------------------
    //double Nexp = contig.frag_depth.Stats.h.mean;
    // contig alignability estimate
    //float totSites = contig.Length-int(contig.totalNoCovBases);
    //float acontig  = float (totSites-contig.totalRepeatBases) / totSites;
    //double LF = pars.getFragmentLength();
    //double Nexp = contig.localpairs.size()*double(LF-2*contig.aLR)/double(contig.Length);
    //----------------------------------------------------------------------------
    // loop over long pair clusters to identify candidate deletions
    //----------------------------------------------------------------------------
    it0=clus.invert3c.cluster.begin();
    it1=clus.invert3c.cluster.end();
    inv1.typeName ="inv3";  
    if (dir==5)  {
        inv1.typeName ="inv5";  
        it0=clus.invert5c.cluster.begin();
        it1=clus.invert5c.cluster.end();
    }
    for ( it=it0 ; it != it1; it++ ) {
        int i = (*it).first;
        if (dir==5) {
            c1 = clus.invert5c.cluster[i];   
        } else {
            c1 = clus.invert3c.cluster[i];   
        }
        //-------------------------------------------------------------------------
        // number of pairs in cluster
        //-------------------------------------------------------------------------
        int Np = c1.inp.size();
        //-------------------------------------------------------------------------
        // initialize cluster limits:    p5a--p5b      p3a-p3b
        //-------------------------------------------------------------------------
        int p5b = 0;              // trailing edge of 5' cluster
        int p5a = contig.Length;  // leading edge of 5' cluster
        int p3b = 0;              // trailing edge of 3' cluster
        int p3a = contig.Length;  // leading edge of 3' cluster
        int p5aMax=0;             // leading edge of trailing read in 5' cluster
        int p3aMax=0;             // leading edge of trailing read in 3' cluster
        int LFmax=0;              // longest library in this cluster
        unsigned int RGCmax=0;    // ReadGroupCode for longest library
        int LFrange=0;            // LF width of longest library
        double qAberr=0;          // product of fragment consistency with LF
        C_RGmap RGmap;            // declare empty Readgroup counter
        //
        cpair.resize(Np);
        for (int j=0; j<Np; j++) {
            int k = c1.inp[j];
            C_localpair lp1;
            if(dir==5) {
                lp1 = clus.invert5[k];
            } else {
                lp1 = clus.invert3[k];
            }
            cpair[j]=lp1;
            int p = lp1.pos;                // start of 5' cluster
            if (p<p5a) p5a=p;
            if (p>p5aMax) p5aMax=p;
            p = lp1.pos+lp1.len1;           // end of 5' cluster
            if (p>p5b) p5b=p;
            p = lp1.pos+lp1.lm-lp1.len2;    // start of 3' cluster
            if (p<p3a) p3a=p;
            if (p>p3aMax) p3aMax=p;
            p = lp1.pos+lp1.lm;             // end of 3' cluster
            if (p>p3b) p3b=p;
            //
            unsigned int RGC1=lp1.ReadGroupCode;
            RGmap[RGC1]++;
            int LMhigh1=libraries.libmap[RGC1].LMhigh;
            if (LMhigh1>LFmax) {
                RGCmax=RGC1;
                LFmax=LMhigh1;
                LFrange=LMhigh1-libraries.libmap[RGC1].LMlow;
            }
            double pAb1=libraries.libmap[RGC1].fragHist.x2pTrim(double(lp1.lm));
            qAberr+=double(p2q(1.0-pAb1))/Np;
        }
        
        
        
        //-------------------------------------------------------------------------
        // estimate gap size where breakpoint should be (Derek's exponential dist?)
        //-------------------------------------------------------------------------
        //int avegap = ((p5b-p5a)+(p3b-p3a))/(2*Np);
        int avegap = ((p5aMax-p5a)+(p3aMax-p3a))/(2*Np);
        // average gap over this contig
        //int avegapC = double(contig.Length)/int(2*contig.localpairs.size());
        //int avegapC = double(contig.Length)/double(contig.totalUniqueReads);
        //-------------------------------------------------------------------------
        // start of inversion is begining of leading cluster 5'
        // length <~ distance from the start of 3' cluster and end of 5' cluster
        //-------------------------------------------------------------------------
        int p0 = p5a+1;
        int len = int(p3a-p5a-1);
        if (dir==5) {
            p0= p5b+1;
            len = int(p3b-p5b-1);
        }
        
        
        //-------------------------------------------------------------------------
        // correct p0 and len for expected gaps between reads 
        //-------------------------------------------------------------------------
        //p0 = p0-avegap/2;
        //len = len+avegap/2;
        //-------------------------------------------------------------------------
        // bail out for negative length events
        //-------------------------------------------------------------------------
        if (len<1) continue; 
        //-------------------------------------------------------------------------
        // check for overlap with previous evt
        //-------------------------------------------------------------------------
        bool overlap = false;
        C_SVV1 e0;
        if (inv1.evt.size()>0) { 
            e0 = inv1.evt.back();
            int pp = e0.pos+e0.length;
            overlap=pp>p0;
            overlap = overlap&&(abs(int(e0.pos)-int(p0))<LFmax);
            /*
             if (overlap) {
             cout << " overlap: prev pos " << e0.pos << "\t prev end " 
             << e0.pos+e0.length << "\t this pos " << p0 << endl;
             }
             */
        }
        //-------------------------------------------------------------------------
        // estimate end of INV region 
        //-------------------------------------------------------------------------  
        // x18Mar2011 int p1 = p0+len-1;
        //-------------------------------------------------------------------------
        // calc read coverage between p0 and p1
        //-------------------------------------------------------------------------  
        //C_SVcoverage1 cov(contig, p0, p1,nomcov) ;
        // coverage within 5' cluster region
        //C_SVcoverage1 cov5(contig, p5a, p5b,nomcov) ;
        // coverage within 3' cluster region
        //C_SVcoverage1 cov3(contig, p3a, p3b,nomcov) ;
        //-------------------------------------------------------------------------
        // candidate SV duplication event 
        //-------------------------------------------------------------------------  
        C_SVV1 e1;
        e1.pos = p0;
        e1.anchor = contig.getAnchorIndex();
        e1.length = len;
        // prob that null CNV exceeds this coverage is 1-cov.p
        //cov.p=fabs(1.0-cov.p);
        e1.q = droundi((100.*Np)/(5+Np)); //p2q(cov.p);
        //    
        //e1.cov=cov;  
        if (dir==5)  {
            e1.cls5 = c1;
            e1.pair5=cpair;
            //e1.cov5=cov5;
            // ReadGroups in event
            e1.ReadGroupMap5=RGmap;
            // end alignability
            //e1.a5[0] = cov5.Nsite/(float(int(cov5.p1)-int(cov5.p0))*acontig);
            //e1.a3[0] = cov3.Nsite/(float(int(cov3.p1)-int(cov3.p0))*acontig);
            
            e1.NfragCluster[0]=Np;
            
        } else {
            e1.cls3 = c1;
            e1.pair3=cpair;
            //e1.cov3=cov3;
            // ReadGroups in event
            e1.ReadGroupMap3=RGmap;
            // end alignability
            //e1.a5[1] = cov5.Nsite/(float(int(cov5.p1)-int(cov5.p0))*acontig);
            //e1.a3[1] = cov3.Nsite/(float(int(cov3.p1)-int(cov3.p0))*acontig);
            
            e1.NfragCluster[1]=Np;
            
        }
        e1.p5[0]=p5a;
        e1.p5[1]=p5b;
        e1.p3[0]=p3a;
        e1.p3[1]=p3b;
        
        // brkpoint uncertainty depends on average gap between read starts in cluster
        //e1.posU=(avegap>avegapC? avegap: avegapC);
        e1.posU=avegap;
        e1.lenU=e1.posU;
        
        // correction to pos,len for deletion bracketed by repeat region
        int DL=(len-droundi(c1.high[1]));
        if (DL>(2*int(e1.posU))) {
            e1.lenU=droundi(c1.std[1]/sqrt(c1.N));
            e1.length=droundi(c1.high[1]);
            e1.posU=DL/2;
            e1.pos=e1.pos+DL/2;
        }
        
        // prob(outlier) in cluster is x2p of fragment dist with range of cluster
        double aRange=((p3b-p3a)+(p5b-p5a))/2.0;
        double pOut=libraries.libmap[RGCmax].fragHist.x2pTrim(aRange);
        e1.qOutlier=p2q(1.0-pOut);
        
        // prob(Aberrant LM) for independent fragments
        e1.qAberrantLM=int(qAberr);
        
        // fragments spanning break region 
        /*
		 e1.NfragCov=int(contig.frag_depth.n[p0]);
         if (p0>LFmax) {
         e1.NfragCovOut[0]=int(contig.frag_depth.n[p0-int(LFmax)]);
         } else {
         e1.NfragCovOut[0]=0;
         }
         if ((int(e1.pos)+len+LFmax)<int(contig.Length)) {
         e1.NfragCovOut[1]=int(contig.frag_depth.n[e1.pos+len+int(LFmax)]);
         } else {
         e1.NfragCovOut[1]=0;
         }
         e1.NfragCovExp=int(Nexp);
         */
		
        //-------------------------------------------------------------------------
        // if overlapped merge this cluster with prev event
        //-------------------------------------------------------------------------
        if (overlap) { 
            C_SVV1 e2 = merge(e0,e1,contig);
            // remove existing event at end of list
            inv1.evt.pop_back();
            e1 = e2;
        } 
        int Npairs = e1.pair5.size()+e1.pair3.size();
        // bkp piles up at the corner of Np-Lm space close to 0,LF
        bool significant = double(Npairs) > 0; //(0.01*Nexp); // + (e1.cls.mean[1]-LF)/LFsig);
        // simple cut on number of pairs in event
        significant = significant & (Npairs>=pars.getMinClustered());
        // simple cut on event length
        significant = significant & (int(e1.length)>=pars.getMinLength());    
        // quantized copy number
        //e1.copynumber = short(2*e1.cov.N/e1.cov.eN +0.5);
        
        // q value is set by supporting fragments (max 40)
        double qNF = 40*double(Npairs)/(5.0+double(Npairs));    
        // q value from length here
        //double qL  = (e1.length>0 ? double(e1.length)/(100.0+fabs(double(e1.length))): 0.1);           
        double qL  = double(e1.length)/(100.0+fabs(double(e1.length)));   
        // q value combined
        double q  = sqrt(qNF*qL);
        
        // q value is set by supporting fragments (max 40)
        e1.q = int(20.0*q/(0.1+q)); 
        
        //-------------------------------------------------------------------------
        // new inversion event
        //-------------------------------------------------------------------------
        //if ((double(e1.cov.N)>(0.1*e1.cov.eN))&significant) {
		if (significant) {
            inv1.evt.push_back(e1);
        }
    }
    
    // sort list
    inv1.evt.sort();
    
    return inv1;
}       

int C_SpannerSV::findCross(C_contig  & contig, C_SpannerCluster & clus,  RunControlParameters & pars) {
    
    C_SVX cross3 = findCrossDir(contig, clus, pars,3);
    C_SVX cross5 = findCrossDir(contig, clus, pars,5);
    
    crx.typeName="crossChromosome";
    
    //----------------------------------------------------------------------------
    // samples name vector sorted & unique in ret
    //----------------------------------------------------------------------------
    C_librarymap::iterator it;
    vector<string>::iterator is;
    crx.samples.clear();
    for ( it=libraries.libmap.begin() ; it != libraries.libmap.end(); it++ )
    {
        C_libraryinfo lib1 = (*it).second; 
        crx.samples.push_back(lib1.Info.SampleName);
    }  
    sort(crx.samples.begin(),crx.samples.end());
    is = unique (crx.samples.begin(), crx.samples.end()); 
    crx.samples.resize( is - crx.samples.begin() );  
    
    //----------------------------------------------------------------------------
    // SVCF Info (list here and in C_SVR << )
    //----------------------------------------------------------------------------
	/*
     string s1[] = {"L","NF","NP","UP","UL","PA","PB","PC","PD","AL","NR","ER","MR"};
     vector<string> Info1(s1, s1 + 13);
     crx.SVCF.Info=Info1;
     
     // SVCF Format  (list here and in C_SVR << operator method)
     string s2[] = {"NF","NP","N3","N5","NR","ER"}; //"CN"};
     vector<string> Format1(s2, s2 + 6);
     crx.SVCF.Format=Format1;
     */
	C_SVCF_TAG tag1;
	tag1.ALT=true;
	tag1.Id="X";
	tag1.Descr="Cross chromosomal adjacency relative to REF";
    crx.SVCF.ALT.push_back(tag1);	
	tag1.ALT=false;
	tag1.INFO=true;
	tag1.Id="MATECHR";
	tag1.Number=1;
	tag1.Type="String";
	tag1.Descr="Chromosome of mate end";
    crx.SVCF.INFO.push_back(tag1);
	tag1.Id="MATEPOS";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Position of mate end";
    crx.SVCF.INFO.push_back(tag1);
    tag1.Id="CIPOS";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Confidence interval around POS";
    crx.SVCF.INFO.push_back(tag1);
    tag1.Id="NFF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of F ALT supporting fragments";
    crx.SVCF.INFO.push_back(tag1);
    tag1.Id="NFR";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of R ALT supporting fragments";
    crx.SVCF.INFO.push_back(tag1);
    tag1.Id="QC";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Outlier fragment metric";
    crx.SVCF.INFO.push_back(tag1);
	tag1.Id="PC1";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Range of F cluster reads";
    crx.SVCF.INFO.push_back(tag1);
	tag1.Id="PC2";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Range of R cluster reads";
    crx.SVCF.INFO.push_back(tag1);
    tag1.INFO=false;
	tag1.FORMAT=true;
	tag1.Id="NFF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of F ALT supporting fragments";
    crx.SVCF.FMT.push_back(tag1);
	tag1.Id="NFR";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of R ALT supporting fragments";
    crx.SVCF.FMT.push_back(tag1);
    
	
    // SVCF samples
    crx.SVCF.Samples=crx.samples;    
    
    // SVCF EventType
    crx.SVCF.EventType=crx.typeName;
    
    
	//double Nexp = contig.frag_depth.Stats.h.mean/2.;
    //----------------------------------------------------------------------------
    // expected number of pairs that span a given deletion
    //----------------------------------------------------------------------------
    int Nmerge = 0;
    list<C_SVX1>::iterator i3,i5,i51;
    C_SVX1 e3,e5,e1;
    
    //----------------------------------------------------------------------------
    // loop over both views (5' and 3') to merge matched inversion events
    //----------------------------------------------------------------------------
    for ( i3=cross3.evt.begin() ; i3 != cross3.evt.end(); i3++ ) {
        e3 = *i3;
        double f1=0;
        for ( i5=cross5.evt.begin() ; i5 != cross5.evt.end(); i5++ ) {
            e5 = *i5;
            if (e5.anchor2!=e3.anchor2) continue;
            if (e5.p5[0]>e3.p3[1]) break;
            //-------------------------------------------------------------------------
            // check for matching lengths 
            //-------------------------------------------------------------------------
            // length over here
            float LH = e3.p3[0]-e5.p5[1];
            // length over there (allow inverted)
            float LT = fabs(double(e5.p3a2[0]-e3.p5a2[1]));
            /// diff
            float DL= fabs(LH-LT);
            // relative diff
            float RL= (LT<LH? DL/LT: DL/LH);
            
            bool matched =( (RL<0.25) & (fabs(DL)<1000) ); // < 50% different in 1kb
            //-------------------------------------------------------------------------
            // check for matching lengths - mutual fractional overlap (1 merge, 0 not)
            //-------------------------------------------------------------------------
            if (matched) {
                f1=1;
                i51=i5;
            }
        }
        e1=e3;
        if (f1>0.05) {    // some mutual overlap threshold
            e5 = *i51;
            e1 = merge(e3,e5,contig);     
            // remove existing event in inv5 list
            cross5.evt.erase(i51); 
            Nmerge++;
        } 
        int Npairs = e1.cross3.size()+e1.cross5.size();
        // bkp piles up at the corner of Np-Lm space close to 0,LF
        //bool significant = double(Npairs) > (0.01*Nexp); //  + (e1.cls.mean[1]-LF)/LFsig);
        // simple cut on number of pairs in event
        // significant = significant & (Npairs>=pars.getMinClustered());    
        bool  significant = (Npairs>=pars.getMinClustered());    
        // no simple cut on event length
        // significant = significant & (int(e1.length)>pars.getMinLength());    
        // quantized copy number
        //e1.copynumber = short(2*e1.cov.N/e1.cov.eN +0.5);
        
        //  length here or over there 
        float Lh = fabs(double(e1.length))+0.1;
        float Lt = double(e1.p3a2[0]-e1.p5a2[1])+0.1;
        // reicprocal difference with length over there (allow inverted)
        float rL = (fabs(Lh)>fabs(Lt)? fabs(Lh-Lt)/fabs(Lh): fabs(Lh-Lt)/fabs(Lt));    
        // q value is set by supporting fragments (max 40)
        double qLr  = (e1.length>=0 ? rL/(0.2+rL): 0.1);    
        // q value is set by supporting fragments (max 40)
        double qNF = double(Npairs)/(5.0+double(Npairs));    
        // q value from length here
        double qL  = (e1.length>=0 ? double(e1.length)/(100.0+fabs(double(e1.length))): 0.1);   
        // q value combined
        double q  = sqrt(qNF*qL*qLr);
        
        // q value is set by supporting fragments (max 40)
        e1.q = int(40.0*q/(0.1+q)); 
        
        //-------------------------------------------------------------------------
        // new cross event
        //-------------------------------------------------------------------------
        //if ((double(e1.cov.N)>(0.1*e1.cov.eN))&significant) { 
        if (significant) { 
            crx.evt.push_back(e1);
        }
    }
    //-------------------------------------------------------------------------
    // remaining cross5 events
    //-------------------------------------------------------------------------
    for ( i5=cross5.evt.begin() ; i5 != cross5.evt.end(); i5++ ) {
        e1 = *i5;
        int Npairs = e1.cross5.size();
        // bkp piles up at the corner of Np-Lm space close to 0,LF
        //bool significant = double(Npairs) > (0.01*Nexp); //  + (e1.cls5.mean[1]-LF)/LFsig); 
        // simple cut on number of pairs in event
        //significant = significant & (Npairs>=pars.getMinClustered());    
        bool significant = (Npairs>=pars.getMinClustered());    
        // quantized copy number
        // e1.copynumber = short(2*e1.cov.N/e1.cov.eN +0.5);
        //-------------------------------------------------------------------------
        // new cross event
        //-------------------------------------------------------------------------
        //if ((double(e1.cov.N)>(0.1*e1.cov.eN))&significant) { 
        if (significant) { 
   			crx.evt.push_back(e1);
        }
    }
    
    // sort events
    crx.evt.sort();
    
    // calc sample supporting and spanning fragments
    crx.finalize(contig,libraries,pars);
    
    // SVCF event count
    crx.SVCF.NEvent=int(crx.evt.size());
    
    return crx.SVCF.NEvent;
    
}       

C_SVX C_SpannerSV::findCrossDir(C_contig  & contig, C_SpannerCluster & clus,  RunControlParameters & pars, int dir) {
    
    C_SVX crx1;
    if (!( (dir==3)||(dir==5) )) return crx1;
    C_cluster2d_elements::iterator it;
    C_cluster2d_elements::iterator it0;
    C_cluster2d_elements::iterator it1;
    C_cluster2d_element1 c1;
    vector<C_crosspair> xPair;
    
    //----------------------------------------------------------------------------
    // fragment length properties
    //----------------------------------------------------------------------------
    // double Nexp = contig.frag_depth.Stats.h.mean;
    // contig alignability estimate
    //float totSites = contig.Length-int(contig.totalNoCovBases);
    //float acontig  = float (totSites-contig.totalRepeatBases) / totSites;
    
    //----------------------------------------------------------------------------
    // loop over long pair clusters to identify candidate deletions
    //----------------------------------------------------------------------------
    it0=clus.cross3c.cluster.begin();
    it1=clus.cross3c.cluster.end();
    crx1.typeName ="cross3";  
    if (dir==5)  {
        crx1.typeName ="cross5";  
        it0=clus.cross5c.cluster.begin();
        it1=clus.cross5c.cluster.end();
    }
    
    //----------------------------------------------------------------------------
    // loop over clusters
    //----------------------------------------------------------------------------
    for ( it=it0 ; it != it1; it++ ) {
        int i = (*it).first;
        if (dir==5) {
            c1 = clus.cross5c.cluster[i];   
        } else {
            c1 = clus.cross3c.cluster[i];   
        }
        //-------------------------------------------------------------------------
        // number of pairs in cluster
        //-------------------------------------------------------------------------
        int Np = c1.inp.size();
        //-------------------------------------------------------------------------
        // initialize cluster limits:    p5a--p5b      p3a-p3b
        //-------------------------------------------------------------------------
        int pHb = 0;              // trailing edge of here cluster
        int pHa = contig.Length;  // leading edge of there cluster
        int pTb = 0;              // trailing edge of here cluster
        int pTa = 500000000L;     // contig.Length;  // leading edge of there cluster
        int pHaMax=0;             // leading edge of trailing read in here cluster
        int pTaMax=0;             // leading edge of trailing read in there cluster
        short aT=-1;              // there chromosome 
        int LFmax=0;              // longest library in this cluster
        unsigned int RGCmax=0;    // ReadGroupCode for longest library
        int LFrange=0;            // LF width of longest library
        C_RGmap RGmap;            // declare empty Readgroup counter
        //
        xPair.resize(Np);
        
        //--------------------------------------------------------------------------
        // 
        //--------------------------------------------------------------------------  
        for (int j=0; j<Np; j++) {
            int k = c1.inp[j];
            C_crosspair xp1;
            if(dir==5) {
                xp1 = clus.cross5[k];
            } else {
                xp1 = clus.cross3[k];
            }
            xPair[j]=xp1;
            int p = xp1.read[0].pos;                // here start of cluster
            if (p<pHa) pHa=p;
            if (p>pHaMax) pHaMax=p;
            p = xp1.read[0].pos+xp1.read[0].len;    // end of cluster
            if (p>pHb) pHb=p;
            
            p = xp1.read[1].pos;                    // there start of cluster 2
            if (p<pTa) pTa=p;
            if (p>pTaMax) pTaMax=p;
            p = xp1.read[1].pos+xp1.read[1].len;    // end of cluster 2
            if (p>pTb) pTb=p;
            
            if ((aT>=0)&(aT!=xp1.read[1].anchor)) {
                cerr << "inconsistent cross cluster" << endl;
            } 
            aT=xp1.read[1].anchor;
            //
            unsigned int RGC1=xp1.ReadGroupCode;
            RGmap[RGC1]++;
            int LMhigh1=libraries.libmap[RGC1].LMhigh;
            if (LMhigh1>LFmax) {
                RGCmax=RGC1;
                LFmax=LMhigh1;
                LFrange=LMhigh1-libraries.libmap[RGC1].LMlow;
            }
        }
        
        //-------------------------------------------------------------------------
        // estimate gap size where breakpoint should be (Derek's exponential dist?)
        //-------------------------------------------------------------------------
        //int avegap = ((p5b-p5a)+(p3b-p3a))/(2*Np);
        int avegap = ((pHaMax-pHa)+(pTaMax-pTa))/(2*Np);
        // average gap over this contig
        //int avegapC = double(contig.Length)/int(2*contig.localpairs.size());
        // int avegapC = double(contig.Length)/double(contig.totalUniqueReads);
        
        //-------------------------------------------------------------------------
        // breakpoint 
        //-------------------------------------------------------------------------
        int pB = pHb+1;
        if (dir==3) {
            pB= pHa-1;
        }
        
        //-------------------------------------------------------------------------
        // estimate end of CRX region 
        //-------------------------------------------------------------------------  
        // x18Mar2011  int p0 = pB-100;
        // x18Mar2011  int p1 = pB+100;
        
		//-------------------------------------------------------------------------
        // calc read coverage between p0 and p1
        //-------------------------------------------------------------------------  
        //C_SVcoverage1 cov(contig, p0, p1,nomcov) ;
        
		// coverage within 5' cluster region
        //C_SVcoverage1 covH(contig, pHa, pHb,nomcov) ;
        //-------------------------------------------------------------------------
        // candidate SV duplication event 
        //-------------------------------------------------------------------------  
        C_SVX1 e1;
        e1.pos = pB;
        e1.anchor = contig.getAnchorIndex();
        e1.length = 0;
        e1.pos2 = (pTa<pTb? pTa: pTb);
        e1.anchor2 = aT+1;
        e1.length2 = 0;
        // prob that null CNV exceeds this coverage is 1-cov.p
        //cov.p=fabs(1.0-cov.p);
        
        // q value is set by supporting fragments (max 20)
        e1.q = int(20.0*Np/(5.0+double(Np)));    
        
        //e1.cov=cov;  
        
        if (dir==5)  {
            e1.cls5 = c1;
            e1.cross5=xPair;
            //e1.cov5=covH;
            // ReadGroups in event
            e1.ReadGroupMap5=RGmap;
            // end alignability
            //e1.a5[0] = covH.Nsite/(float(int(covH.p1)-int(covH.p0))*acontig);      
            e1.NfragCluster[0]=Np;
            e1.p5[0]=pHa;
            e1.p5[1]=pHb;
            e1.p3a2[0]=pTa;
            e1.p3a2[1]=pTb;
        } else {
            e1.cls3 = c1;
            e1.cross3=xPair;
            //e1.cov3=covH;
            // ReadGroups in event
            e1.ReadGroupMap3=RGmap;
            // end alignability
            //e1.a3[1] = covH.Nsite/(float(int(covH.p1)-int(covH.p0))*acontig);
            e1.NfragCluster[1]=Np;
            e1.p3[0]=pHa;
            e1.p3[1]=pHb;
            e1.p5a2[0]=pTa;
            e1.p5a2[1]=pTb;    
        }
        
        // brkpoint uncertainty depends on average gap between read starts in cluster
		
        //e1.posU=(avegap>avegapC? avegap: avegapC);
        e1.posU=avegap;
		
        // prob(outlier) in cluster is x2p of fragment dist with range of cluster
        double aRange=((pTb-pTa)+(pHb-pHa))/2.0;
        double pOut=libraries.libmap[RGCmax].fragHist.x2pTrim(aRange);
        e1.qOutlier=p2q(1.0-pOut);
        
        // fragments spanning break region 
        /*
         e1.NfragCov=int(contig.frag_depth.n[pB]);
         if (pB>LFmax) {
         e1.NfragCovOut[0]=int(contig.frag_depth.n[pB-int(LFmax)]);
         } else {
         e1.NfragCovOut[0]=0;
         }
         if ((int(e1.pos)+LFmax)<int(contig.Length)) {
         e1.NfragCovOut[1]=int(contig.frag_depth.n[e1.pos+int(LFmax)]);
         } else {
         e1.NfragCovOut[1]=0;
         }
         e1.NfragCovExp=int(Nexp);
         */
		
        int Npairs = e1.cross5.size()+e1.cross3.size();
        // bkp piles up at the corner of Np-Lm space close to 0,LF
        //bool significant = double(Npairs) > (0.01*Nexp); // + (e1.cls.mean[1]-LF)/LFsig);
        // simple cut on number of pairs in event
        bool significant =  (Npairs>=pars.getMinClustered());
        // quantized copy number
        //e1.copynumber = short(2*e1.cov.N/e1.cov.eN +0.5);
        //-------------------------------------------------------------------------
        // new inversion event
        //-------------------------------------------------------------------------
        if (significant) { //(double(e1.cov.N)>(0.01*e1.cov.eN))&significant) {
            crx1.evt.push_back(e1);
        }
    }
    
    // sort list
    crx1.evt.sort();
    
    return crx1;
}       




//-------------------------------------------------------------------------
// merge two SV1 events
//------------------------------------------------------------------------- 
C_SV1 C_SpannerSV::merge(C_SV1 & e0, C_SV1 & e1, C_contig  & contig, RunControlParameters & pars, int type) {
    
    // init empty event
    C_SV1 e2;
    e2.merge=true;
    e2.anchor=e0.anchor;
	// vector of fraglength shift
	vector<int> DiffLF;
	
    // ranges
    e2.p5[0] = (e0.p5[0]<e1.p5[0]? e0.p5[0]: e1.p5[0]);
    e2.p5[1] = (e0.p5[1]>e1.p5[1]? e0.p5[1]: e1.p5[1]);
    e2.p3[0] = (e0.p3[0]<e1.p3[0]? e0.p3[0]: e1.p3[0]);
    e2.p3[1] = (e0.p3[1]>e1.p3[0]? e0.p3[1]: e1.p3[1]);
    
    // for deletions (type=0) & duplications (type=1)
    //int p0 = (e0.pos>e1.pos? e0.pos: e1.pos);
    //int len = int(e.p3[0]-e.p5[1]-1);
    //if (type==1) { len = int(e2.p5[1]-e.p3[0]-1);  }
    int N0 = e0.pair.size();
    int N1 = e1.pair.size();
    int Np = e0.pair.size()+e1.pair.size();
    
    // merge read group maps
    // e2.ReadGroupMap=e0.ReadGroupMap;
    e2.ReadGroupMap.clear();
    std::map<unsigned int,unsigned int, less<unsigned int> >::iterator it;
    for ( it=e0.ReadGroupMap.begin() ; it != e0.ReadGroupMap.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        unsigned int NF1  = e0.ReadGroupMap[RGC1];
        e2.ReadGroupMap[RGC1]+=NF1;
    }
    for ( it=e1.ReadGroupMap.begin() ; it != e1.ReadGroupMap.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        unsigned int NF1  = e1.ReadGroupMap[RGC1];
        e2.ReadGroupMap[RGC1]+=NF1;
    }
    
    int LFmax=0;              // longest library in this cluster
    unsigned int RGCmax=0;    // ReadGroupCode for longest library
    int LFrange=0;
    for ( it=e2.ReadGroupMap.begin() ; it != e2.ReadGroupMap.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        //unsigned int NF1  = e2.ReadGroupMap[RGC1];
        int LMhigh1=libraries.libmap[RGC1].LMhigh;
        if (LMhigh1>LFmax) {
            RGCmax=RGC1;
            LFmax=LMhigh1;
            LFrange=LMhigh1-libraries.libmap[RGC1].LMlow;
        }
    }
    
    // merge clusters
    e2.pair = e0.pair;
    e2.pair.insert(e2.pair.end(), e1.pair.begin(), e1.pair.end() );
    C_cluster2d_element1 c1;
    c1.N = e2.pair.size();
    
    // loop over fragments
    c1.low[0]=contig.Length;
    c1.low[1]=contig.Length;
    int p5aMax=0;             // leading edge of trailing read in 5' cluster
    int p3aMax=0;             // leading edge of trailing read in 3' cluster
    
    for (int i=0; i<c1.N; i++)  {
        int p = e2.pair[i].pos;
        if (p>p5aMax) p5aMax=p;
        p = p+e2.pair[i].lm-e2.pair[i].len2;    // start of 3' cluster
        if (p>p3aMax) p3aMax=p;
        for (int j=0; j<2; j++) {
            double x = e2.pair[i].pos;
            if (j==1) x = e2.pair[i].lm;
            c1.mean[j]+=x;
            c1.std[j]+=x*x;
            c1.low[j]=(x>c1.low[j] ? c1.low[j]: x);
            c1.high[j]=(x<c1.high[j] ? c1.high[j]: x);
            c1.inp.push_back(i);
        }
    }
    for (int j=0; j<2; j++) {
        c1.mean[j]=c1.mean[j]/c1.N;
        c1.std[j]=sqrt(c1.std[j]/c1.N-(c1.mean[j]*c1.mean[j]));
    }
    
    e2.cls=c1;
    
    //-------------------------------------------------------------------------
    // estimate gap size where breakpoint should be 
    //-------------------------------------------------------------------------
    //int avegap = ((e.p5[1]-e.p5[0])+(e.p3[1]-e.p3[0]))/(2*Np);
    int avegap = ((p5aMax-e2.p5[0])+(p3aMax-e2.p3[0]))/(2*Np);
    // average gap over this contig
    //int avegapC = double(contig.Length)/int(2*contig.localpairs.size());
    //int avegapC = double(contig.Length)/double(contig.totalUniqueReads);
    //-------------------------------------------------------------------------
    // correct p0 and len for expected gaps between reads 
    //-------------------------------------------------------------------------
    // e2.pos = e2.p5[1]+1+avegap/2;
    // e2.length = int(e2.p3[0]-e2.p5[1]-1)-avegap/2;
    e2.pos = e2.p5[1]+1;
    e2.length = droundi(c1.mean[1]); //int(e2.p3[0]-e2.p5[1]-1);
    if (type==1) { 
        e2.pos = e2.p3[0]+1-avegap/2;
        e2.length = int(e2.p3[1]-e2.p3[0]-1)+avegap/2;
    }
    
    
    //-------------------------------------------------------------------------
    // coverage
    //-------------------------------------------------------------------------
    // x18Mar2011  int p1 = e2.pos+e2.length-1;
    //C_SVcoverage1 cov(contig, e2.pos+1, p1,nomcov) ;
    
    /*
     if (type==1) { cov.p = fabs(1.0-cov.p);  }
     e2.q = p2q(cov.p);
     e2.cov=cov;
     e2.copynumber = short(2*e2.cov.N/e2.cov.eN +0.5);
     
     // coverage within 5' cluster region
     C_SVcoverage1 cov5(contig, e2.p5[0], e2.p5[1],nomcov) ;
     // coverage within 3' cluster region
     C_SVcoverage1 cov3(contig, e2.p3[0], e2.p3[1],nomcov) ;
     e2.cov5=cov5;
     e2.cov3=cov3;
     */
	
	cerr << "merge SV1" << endl;
    cerr << e0.pos << endl;
    cerr << e1.pos << endl;
	cerr << endl;
	
    // brkpoint uncertainty depends on average gap between read starts in cluster
    //e1.posU=(avegap>avegapC? avegap: avegapC);
    e2.CIpos[0]=(p5aMax-e2.p5[0])/Np;
	e2.CIpos[1]=(p3aMax-e2.p3[0])/Np;
    e2.CIend[0]=(p5aMax-e2.p5[0])/Np;
	e2.CIend[1]=(p3aMax-e2.p3[0])/Np;
    e2.CIlen[0]=0;
	e2.CIlen[1]=0;
	// correction to pos,len for deletion bracketed by repeat region
    int DL=(e2.length-droundi(e2.cls.mean[1]));
    if (DL>(4*int(e2.CIpos[0]))) {
        e2.CIlen[0]=droundi(e2.cls.std[1]/sqrt(e2.cls.N));
		e2.CIlen[1]=droundi(e2.cls.std[1]/sqrt(e2.cls.N));
		e2.length=droundi(e2.cls.mean[1]);
		e2.CIpos[0]=-DL/2;
		e2.CIpos[1]=DL/2;
		e2.pos=e2.pos+DL/2;
	}
	
    
    
    double aRange=((e2.p3[1]-e2.p3[0])+(e2.p5[1]-e2.p5[0]))/2.0;
    double pOut=libraries.libmap[RGCmax].fragHist.x2pTrim(aRange);
    
    e2.qOutlier=p2q(1.0-pOut);
    
    //prob(Aberrant LM) for independent fragments
    e2.qAberrantLM=(N0*e0.qAberrantLM+N1*e1.qAberrantLM)/Np;
    
    // average alignability over contig
    //float totSites = contig.Length-int(contig.totalNoCovBases);
    //float acontig  = float (totSites-contig.totalRepeatBases) / totSites;
    // end alignability
    //e2.a5 = cov5.Nsite/(float(int(cov5.p1)-int(cov5.p0))*acontig);
    //e2.a3 = cov3.Nsite/(float(int(cov3.p1)-int(cov3.p0))*acontig);
    
    return e2;
}

//------------------------------------------------------------------------------
// merge adjacent inv clusters
//------------------------------------------------------------------------------
C_SVV1 C_SpannerSV::merge(C_SVV1 & e0, C_SVV1 & e1, C_contig  & contig) {
    
    //----------------------------------------------------------------------------
    // flip e0 to 5 and e1 to R
    //----------------------------------------------------------------------------
    if ( (e1.pair5.size()>0) && (e0.pair5.size()==0) ) {
        // put F cluster in e0, R in e1
        C_SVV1 etemp=e0;
        e0=e1;
        e1=etemp;
    } 
    
    // fragment counts from each event
    int N0 = e0.pair5.size()+e0.pair3.size();
    int N1 = e1.pair5.size()+e1.pair3.size();
    int Np = N0+N1;
    
    //----------------------------------------------------------------------------
    // declare merged event, start with e0
    //----------------------------------------------------------------------------
    C_SVV1 ee;
    ee=e0;
    
    ee.p3[0] = (e0.p3[0]<e1.p3[0]? e0.p3[0]: e1.p3[0]);
    ee.p3[1] = (e0.p3[1]>e1.p3[0]? e0.p3[1]: e1.p3[1]);
    ee.p5[0] = (e0.p5[0]<e1.p5[0]? e0.p5[0]: e1.p5[0]);
    ee.p5[1] = (e0.p5[1]>e1.p5[1]? e0.p5[1]: e1.p5[1]);
    
    // alignability estimates for 3' end of F (0) cluster (ee.a3[0]) or R cluster (1) 
    //ee.a3[0]=e0.a3[0];
    ee.a3[1]=e1.a3[1];
    //ee.a5[0]=e0.a5[0];
    ee.a5[1]=e1.a5[1];
    
    //----------------------------------------------------------------------------
    // check for clustered F's or R's
    //----------------------------------------------------------------------------
    bool bothR= (e0.pair5.size()==0) && (e1.pair5.size()==0) ;
    bool bothF= (e0.pair3.size()==0) && (e1.pair3.size()==0) ;
    if (bothF) {     
        //ee.pair5 = e0.pair5;
        ee.pair5.insert(ee.pair5.end(), e1.pair5.begin(), e1.pair5.end() );
        ee.merge5=true;
        
    } else if (bothR) {
        // ee.pair3 = e0.pair3;
        ee.pair3.insert(ee.pair3.end(), e1.pair3.begin(), e1.pair3.end() );
        ee.merge3=true;
        
    } else { // merge F and R
        
        // recalc cluster limits for both F and R
        // start of p5 end
        ee.p5[0] = e0.p5[0];
        // estimate break at inv start
        ee.p5[1] = (e0.p5[1]*N0+e1.p5[0]*N1)/Np;
        // estimate break at inv end
        ee.p3[0] = (e0.p3[1]*N0+e1.p3[0]*N1)/Np;
        ee.p3[1] = e1.p3[1];
        //ee.cls5  = e0.cls5;
        
        ee.cls3  = e1.cls3;
        /* / debug vvv
         int N05=e0.pair5.size();
         int N15=e1.pair5.size();
         int N03=e0.pair3.size();
         int N13=e1.pair3.size();
         cout<< N05 <<"\t" << N15 <<"\t" << N03 <<"\t" << N13 <<endl;
         // debug ^^^
         */
        ee.pair5 = e0.pair5;
        ee.pair5.insert(ee.pair5.end(), e1.pair5.begin(), e1.pair5.end() );
        ee.pair3 = e0.pair3;
        ee.pair3.insert(ee.pair3.end(), e1.pair3.begin(), e1.pair3.end() );
        /* / debug vvv
         int N3=ee.pair3.size();
         int N5=ee.pair5.size();
         cout<< N5 <<"\t" << N3 <<endl;
         // debug ^^^
         */
        ee.merge3=e1.merge3;
        
    }
    
    
    //----------------------------------------------------------------------------
    // position of event
    //----------------------------------------------------------------------------
    ee.pos = ee.p5[1]+1;
    ee.length = int(ee.p3[0])-int(ee.pos);
    
    
    //ee.pos = (e0.pos*N0+e1.pos*N5)/Np;  // weighted average
    //ee.length = int(e0.length*N0+e1.length*N1)/Np;
    ee.NfragCluster[0] = e0.pair5.size()+e1.pair5.size();
    ee.NfragCluster[1] = e0.pair3.size()+e1.pair3.size();
    int p0 = (int(ee.pos)>100? int(ee.pos)-100: 0);  
    // x18Mar2011  int p1 = (int(ee.pos)<(contig.Length-101)? ee.pos+100: contig.Length-1);
    
	/*
     //-------------------------------------------------------------------------
     // calc read coverage between p0 and p1
     //-------------------------------------------------------------------------
     C_SVcoverage1 cov(contig, p0, p1,nomcov) ;
     // coverage within 5' cluster region
     C_SVcoverage1 cov5(contig, ee.p5[0], ee.p5[1],nomcov) ;
     // coverage within 3' cluster region
     C_SVcoverage1 cov3(contig, ee.p3[0], ee.p3[1],nomcov) ;
     
     ee.cov=cov;
     ee.cov5=cov5;
     ee.cov3=cov3; 
     */
	
    // merge read group info
    std::map<unsigned int,unsigned int, less<unsigned int> >::iterator it;
    /* already did e0 
     ee.ReadGroupMap3.clear();
     for ( it=e0.ReadGroupMap3.begin() ; it != e0.ReadGroupMap3.end(); it++ ) {
     unsigned int RGC1 = (*it).first;
     unsigned int NF1  = e0.ReadGroupMap3[RGC1];
     ee.ReadGroupMap3[RGC1]+=NF1;
     }
     */
    for ( it=e1.ReadGroupMap3.begin() ; it != e1.ReadGroupMap3.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        unsigned int NF1  = e1.ReadGroupMap3[RGC1];
        ee.ReadGroupMap3[RGC1]+=NF1;
    }
    // merge read group info
    /*
     ee.ReadGroupMap5.clear();
     for ( it=e0.ReadGroupMap5.begin() ; it != e0.ReadGroupMap5.end(); it++ ) {
     unsigned int RGC1 = (*it).first;
     unsigned int NF1  = e0.ReadGroupMap5[RGC1];
     ee.ReadGroupMap5[RGC1]+=NF1;
     }
     */
    for ( it=e1.ReadGroupMap5.begin() ; it != e1.ReadGroupMap5.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        unsigned int NF1  = e1.ReadGroupMap5[RGC1];
        ee.ReadGroupMap5[RGC1]+=NF1;
    }
    
    // Max lib LF in this event
    int LFmax=0;              // longest library in this cluster
    unsigned int RGCmax=0;    // ReadGroupCode for longest library
    int LFrange=0;
    for ( it=ee.ReadGroupMap3.begin() ; it != ee.ReadGroupMap3.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        //unsigned int NF1  = ee.ReadGroupMap3[RGC1];
        int LMhigh1=libraries.libmap[RGC1].LMhigh;
        if (LMhigh1>LFmax) {
            RGCmax=RGC1;
            LFmax=LMhigh1;
            LFrange=LMhigh1-libraries.libmap[RGC1].LMlow;
        }
    }
    for ( it=ee.ReadGroupMap5.begin() ; it != ee.ReadGroupMap5.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        //unsigned int NF1  = ee.ReadGroupMap5[RGC1];
        int LMhigh1=libraries.libmap[RGC1].LMhigh;
        if (LMhigh1>LFmax) {
            RGCmax=RGC1;
            LFmax=LMhigh1;
            LFrange=LMhigh1-libraries.libmap[RGC1].LMlow;
        }
    }
    
    //----------------------------------------------------------------------------
    // clustering  
    //----------------------------------------------------------------------------
    
    vector<C_localpair> pair;
    if (bothF||bothR) {
        int p5aMax=0;           // 5' cluster
        int p3aMax=0;           // 3' cluster
        C_cluster2d_element1 c1;
        if (bothF) {
            pair=ee.pair5;
        }
        if (bothR) {
            pair=ee.pair3;
        }
        c1.N=pair.size();
        for (int j=0; j<2; j++) {
            c1.low[j]=contig.Length;
            for (int i=0; i<c1.N; i++)  {
                double x = pair[i].pos+pair[i].lm/2;
                //if (bothF) {x+= pair[i].len1;}
                if (j==1) {
                    int LF1 = libraries.libmap[pair[i].ReadGroupCode].LM;
                    x = pair[i].lm-LF1;
                } else {
                    if (int(pair[i].pos)>p5aMax) p5aMax=pair[1].pos;
                    int p3x=pair[i].pos+pair[i].lm-pair[i].len2;
                    if (p3x>p3aMax) p3aMax=p3x;
                }
                c1.mean[j]+=x;
                c1.std[j]+=x*x;
                c1.low[j]=(x>c1.low[j] ? c1.low[j]: x);
                c1.high[j]=(x<c1.high[j] ? c1.high[j]: x);
                c1.inp.push_back(i);
            }
            c1.mean[j]=c1.mean[j]/c1.N;
            c1.std[j]=sqrt(c1.std[j]/c1.N-(c1.mean[j]*c1.mean[j]));
        }
        if (bothF) {
            ee.cls5=c1;
        }  else {
            ee.cls3=c1;
        }   
    }
    
    
    
    //----------------------------------------------------------------------------
    // local copy number
    //----------------------------------------------------------------------------
    //ee.copynumber = short(2.0*ee.cov.N/ee.cov.eN +0.5);  
    
    
    //----------------------------------------------------------------------------
    // fragment coverage at and around insertion
    //----------------------------------------------------------------------------
    p0=ee.pos;
    
	// x18Mar2011  int len=ee.length;
    /*
	 ee.NfragCov=int(contig.frag_depth.n[p0]);
     if (p0>LFmax) {
     ee.NfragCovOut[0]=int(contig.frag_depth.n[p0-int(LFmax)]);
     } else {
     ee.NfragCovOut[0]=0;
     }
     if ((p0+len+LFmax)<contig.Length) {
     ee.NfragCovOut[1]=int(contig.frag_depth.n[p0+len+int(LFmax)]);
     } else {
     ee.NfragCovOut[0]=0;
     }
     double Nexp = contig.frag_depth.Stats.h.mean;
     ee.NfragCovExp=int(Nexp);
     */
	
    //
    // brkpoint uncertainty
    ee.posU=droundi(float(Np)/(float(N0)/e0.posU+float(N1)/e1.posU));
    ee.lenU=ee.posU;
    
    double aRange=((ee.p3[1]-ee.p3[0])+(ee.p5[1]-ee.p5[0]))/2.0;
    double pOut=libraries.libmap[RGCmax].fragHist.x2pTrim(aRange);
    
    ee.qOutlier=p2q(1.0-pOut);
    
    // prob(Aberrant LM) for independent fragments
    ee.qAberrantLM=(N0*e0.qAberrantLM+N1*e1.qAberrantLM)/Np;
    
    ee.q = 0;
    
    C_SVV1 eout = ee;
    return eout;
}


//------------------------------------------------------------------------------
// merge adjacent cross events
//------------------------------------------------------------------------------
C_SVX1 C_SpannerSV::merge(C_SVX1 & e0, C_SVX1 & e1, C_contig  & contig) {
    
    //----------------------------------------------------------------------------
    // flip e0 to 5 and e1 to R
    //----------------------------------------------------------------------------
    if ( (e1.cross5.size()>0) && (e0.cross5.size()==0) ) {
        // put F cluster in e0, R in e1
        C_SVX1 etemp=e0;
        e0=e1;
        e1=etemp;
    } 
    
    // fragment counts from each event
    int N0 = e0.cross5.size()+e0.cross3.size();
    int N1 = e1.cross5.size()+e1.cross3.size();
    int Np = N0+N1;
    
    //----------------------------------------------------------------------------
    // declare merged event, start with e0
    //----------------------------------------------------------------------------
    C_SVX1 ee;
    ee=e0;
    
    
    //----------------------------------------------------------------------------
    // check for clustered F's or R's
    //----------------------------------------------------------------------------
    bool bothR= (e0.cross5.size()==0) && (e1.cross5.size()==0) ;
    bool bothF= (e0.cross3.size()==0) && (e1.cross3.size()==0) ;
    if (bothF) {     
        //ee.pair5 = e0.pair5;
        ee.cross5.insert(ee.cross5.end(), e1.cross5.begin(), e1.cross5.end() );
        ee.merge5=true;
        ee.p5[0] = (e0.p5[0]<e1.p5[0]? e0.p5[0]: e1.p5[0]);
        ee.p5[1] = (e0.p5[1]>e1.p5[1]? e0.p5[1]: e1.p5[1]);
        ee.p3a2[0] = (e0.p3a2[0]<e1.p3a2[0]? e0.p3a2[0]: e1.p3a2[0]);
        ee.p3a2[1] = (e0.p3a2[1]>e1.p3a2[1]? e0.p3a2[1]: e1.p3a2[1]);
        
    } else if (bothR) {
        // ee.pair3 = e0.pair3;
        ee.cross3.insert(ee.cross3.end(), e1.cross3.begin(), e1.cross3.end() );
        ee.merge3=true;
        ee.p3[0] = (e0.p3[0]<e1.p3[0]? e0.p3[0]: e1.p3[0]);
        ee.p3[1] = (e0.p3[1]>e1.p3[0]? e0.p3[1]: e1.p3[1]);
        ee.p5a2[0] = (e0.p5a2[0]<e1.p5a2[0]? e0.p5a2[0]: e1.p5a2[0]);
        ee.p5a2[1] = (e0.p5a2[1]>e1.p5a2[1]? e0.p5a2[1]: e1.p5a2[1]);
        
    } else { // merge F and R
        
        ee.p3[0] = e1.p3[0];
        ee.p3[1] = e1.p3[1];
        ee.p5[0] = e0.p5[0];
        ee.p5[1] = e0.p5[1];
        ee.p3a2[0] = e0.p3a2[0];
        ee.p3a2[1] = e0.p3a2[1];
        ee.p5a2[0] = e1.p5a2[0];
        ee.p5a2[1] = e1.p5a2[1];
        
        // alignability estimates for 3' end of F (0) cluster (ee.a3[0]) or R cluster (1) 
        //ee.a3[0]=e0.a3[0];
        ee.a3[0]=e1.a3[0];
        ee.a3[1]=e1.a3[1];
        //ee.a5[0]=e0.a5[0];
        ee.a5[0]=e0.a5[0];
        ee.a5[1]=e0.a5[1];
        ee.cls3  = e1.cls3;
        ee.cross3 = e1.cross3;
        ee.length = int(ee.p3[0])-int(ee.p5[1]);
        ee.length2 = int(ee.p3a2[0])-int(ee.p5a2[1]);
    }
    
    
    //----------------------------------------------------------------------------
    // position of event
    //----------------------------------------------------------------------------
    ee.pos = ee.p5[1]+1;
    
    
    ee.NfragCluster[0] = e0.cross5.size()+e1.cross5.size();
    ee.NfragCluster[1] = e0.cross3.size()+e1.cross3.size();
    int p0 = (int(ee.pos)>100? int(ee.pos)-100: 0);  
    int p1 = (int(ee.pos)<(contig.Length-101)? ee.pos+100: contig.Length-1);
    if (ee.length>0) {
        p0=ee.pos+1;
        p1=p0+ee.length-1;
    }
    /*
     //-------------------------------------------------------------------------
     // calc read coverage between p0 and p1
     //-------------------------------------------------------------------------
     C_SVcoverage1 cov(contig, p0, p1,nomcov) ;
     // coverage within 5' cluster region
     C_SVcoverage1 cov5(contig, ee.p5[0], ee.p5[1],nomcov) ;
     // coverage within 3' cluster region
     C_SVcoverage1 cov3(contig, ee.p3[0], ee.p3[1],nomcov) ;
     
     ee.cov=cov;
     ee.cov5=cov5;
     ee.cov3=cov3;
     */
	
    // merge read group info
    std::map<unsigned int,unsigned int, less<unsigned int> >::iterator it;
    for ( it=e1.ReadGroupMap3.begin() ; it != e1.ReadGroupMap3.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        unsigned int NF1  = e1.ReadGroupMap3[RGC1];
        ee.ReadGroupMap3[RGC1]+=NF1;
    }
    // merge read group info
    for ( it=e1.ReadGroupMap5.begin() ; it != e1.ReadGroupMap5.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        unsigned int NF1  = e1.ReadGroupMap5[RGC1];
        ee.ReadGroupMap5[RGC1]+=NF1;
    }
    
    // Max lib LF in this event
    int LFmax=0;              // longest library in this cluster
    unsigned int RGCmax=0;    // ReadGroupCode for longest library
    int LFrange=0;
    for ( it=ee.ReadGroupMap3.begin() ; it != ee.ReadGroupMap3.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        //unsigned int NF1  = ee.ReadGroupMap3[RGC1];
        int LMhigh1=libraries.libmap[RGC1].LMhigh;
        if (LMhigh1>LFmax) {
            RGCmax=RGC1;
            LFmax=LMhigh1;
            LFrange=LMhigh1-libraries.libmap[RGC1].LMlow;
        }
    }
    for ( it=ee.ReadGroupMap5.begin() ; it != ee.ReadGroupMap5.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        //unsigned int NF1  = ee.ReadGroupMap5[RGC1];
        int LMhigh1=libraries.libmap[RGC1].LMhigh;
        if (LMhigh1>LFmax) {
            RGCmax=RGC1;
            LFmax=LMhigh1;
            LFrange=LMhigh1-libraries.libmap[RGC1].LMlow;
        }
    }
    
    //----------------------------------------------------------------------------
    // clustering  
    //----------------------------------------------------------------------------
    
    vector<C_crosspair> cross;
    if (bothF||bothR) {
        int p5aMax=0;           // 5' cluster
        int p3aMax=0;           // 3' cluster
        C_cluster2d_element1 c1;
        if (bothF) {
            cross=ee.cross5;
        }
        if (bothR) {
            cross=ee.cross3;
        }
        c1.N=cross.size();
        for (int j=0; j<2; j++) {
            c1.low[j]=contig.Length;
            for (int i=0; i<c1.N; i++)  {
                double x = (cross[i].read[0].pos);
                if (j==1) {
                    x = (cross[i].read[1].pos);
                } else {
                    if (int(cross[i].read[0].pos)>p5aMax) p5aMax=cross[1].read[0].pos;
                    int p3x=cross[i].read[1].pos;
                    if (p3x>p3aMax) p3aMax=p3x;
                }
                c1.mean[j]+=x;
                c1.std[j]+=x*x;
                c1.low[j]=(x>c1.low[j] ? c1.low[j]: x);
                c1.high[j]=(x<c1.high[j] ? c1.high[j]: x);
                c1.inp.push_back(i);
            }
            c1.mean[j]=c1.mean[j]/c1.N;
            c1.std[j]=sqrt(c1.std[j]/c1.N-(c1.mean[j]*c1.mean[j]));
        }
        if (bothF) {
            ee.cls5=c1;
        }  else {
            ee.cls3=c1;
        }   
    }
    
    
    
    //----------------------------------------------------------------------------
    // local copy number
    //----------------------------------------------------------------------------
    //ee.copynumber = short(2.0*ee.cov.N/ee.cov.eN +0.5);  
    
    
    //----------------------------------------------------------------------------
    // fragment coverage at and around insertion
    //----------------------------------------------------------------------------
    p0=ee.pos;
    
	// x18Mar2011  int len=ee.length;
	
	/*
	 ee.NfragCov=int(contig.frag_depth.n[p0]);
     if (p0>LFmax) {
     ee.NfragCovOut[0]=int(contig.frag_depth.n[p0-int(LFmax)]);
     } else {
     ee.NfragCovOut[0]=0;
     }
     
     if ((p0+len+LFmax)<contig.Length) {
     ee.NfragCovOut[1]=int(contig.frag_depth.n[p0+len+int(LFmax)]);
     } else {
     ee.NfragCovOut[0]=0;
     }
     double Nexp = contig.frag_depth.Stats.h.mean;
     ee.NfragCovExp=int(Nexp);
	 */
    
    //
    // brkpoint uncertainty
    ee.posU=droundi(float(Np)/(float(N0)/e0.posU+float(N1)/e1.posU));
    ee.lenU=ee.posU;
    
    double aRange=((ee.p3[1]-ee.p3[0])+(ee.p5[1]-ee.p5[0]))/2.0;
    double pOut=libraries.libmap[RGCmax].fragHist.x2pTrim(aRange);
    
    ee.qOutlier=p2q(1.0-pOut);
    
    ee.q = 0;
    
    C_SVX1 eout = ee;
    return eout;
}




//-------------------------------------------------------------------------
// element insertion event
//------------------------------------------------------------------------- 
int C_SpannerSV::findRet(C_contig  & contig, C_SpannerRetroCluster & rclus,  
                         RunControlParameters & pars) {
    ret.evt.clear();
    ret.typeName=rclus.typeName;
    C_cluster2d_elements::iterator i1;
    C_cluster2d_element1 c1;
    
    //----------------------------------------------------------------------------
    // samples name vector sorted & unique in ret
    //----------------------------------------------------------------------------
    C_librarymap::iterator it;
    vector<string>::iterator is;
    ret.samples.clear();
    for ( it=libraries.libmap.begin() ; it != libraries.libmap.end(); it++ )
    {
        C_libraryinfo lib1 = (*it).second; 
        ret.samples.push_back(lib1.Info.SampleName);
    }  
    sort(ret.samples.begin(),ret.samples.end());
    is = unique (ret.samples.begin(), ret.samples.end()); 
    ret.samples.resize( is - ret.samples.begin() );  
    
    //----------------------------------------------------------------------------
    // SVCF Info (list here and in C_SVR << )
    //----------------------------------------------------------------------------
	/*
     string s1[] = {"L","NF","NP","UP","UL","PA","PB","PC","PD","AL","NR","ER","MR"};
     vector<string> Info1(s1, s1 + 13);
     ret.SVCF.Info=Info1;
     
     // SVCF Format  (list here and in C_SVR << operator method)
     string s2[] = {"NF","NP","N3","N5","NR","ER"}; //"CN"};
     vector<string> Format1(s2, s2 + 6);
     ret.SVCF.Format=Format1;
     */
	C_SVCF_TAG tag1;
	
	tag1.ALT=true;
	tag1.Id="INS:"+ret.typeName;
	tag1.Descr="Insertion relative to REF: "+ret.typeName;
    ret.SVCF.ALT.push_back(tag1);
	//
	tag1.ALT=false;
	tag1.INFO=true;
	tag1.Id="SVLEN";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Difference in length between REF and ALT alleles";
    ret.SVCF.INFO.push_back(tag1);
    tag1.Id="CIPOS";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Confidence interval around POS";
    ret.SVCF.INFO.push_back(tag1);
    tag1.Id="NFF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of F ALT supporting fragments";
    ret.SVCF.INFO.push_back(tag1);
    tag1.Id="QF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Aberrant fragment length metric";
    ret.SVCF.INFO.push_back(tag1);
    tag1.Id="QC";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Outlier fragment metric";
    ret.SVCF.INFO.push_back(tag1);
	tag1.Id="PC1";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Range of F cluster reads";
    ret.SVCF.INFO.push_back(tag1);
	tag1.Id="PC2";
	tag1.Number=2;
	tag1.Type="Integer";
	tag1.Descr="Range of R cluster reads";
    ret.SVCF.INFO.push_back(tag1);
    tag1.INFO=false;
	tag1.FORMAT=true;
	tag1.Id="NFF";
	tag1.Number=1;
	tag1.Type="Integer";
	tag1.Descr="Number of F ALT supporting fragments";
    ret.SVCF.FMT.push_back(tag1);
	tag1.Id="NFR";
	tag1.Descr="Number of R ALT supporting fragments";
    ret.SVCF.FMT.push_back(tag1);
    tag1.FILTER=true;
	tag1.FORMAT=false;
	tag1.Id="1_SIDE";
	tag1.Descr="Insertion supporting fragments from only one side";
    ret.SVCF.FILT.push_back(tag1);
	
	
    // SVCF samples
    ret.SVCF.Samples=ret.samples;    
    
    // SVCF EventType
    ret.SVCF.EventType=ret.typeName;
    
    //----------------------------------------------------------------------------
    // fragment length properties
    //----------------------------------------------------------------------------
    int LFmax=libraries.maxLF();
    
    //----------------------------------------------------------------------------
    // expected number of pairs that span a given deletion
    //----------------------------------------------------------------------------
    double Nexp = contig.frag_depth.Stats.h.mean;
    
    //----------------------------------------------------------------------------
    // remove constrained clusters
    //----------------------------------------------------------------------------
    C_NNcluster2d c5c = cullRet(contig, rclus.e5c, rclus.e5,  pars);
    C_NNcluster2d c3c = cullRet(contig, rclus.e3c, rclus.e3,  pars);
    
    //----------------------------------------------------------------------------
    // remove clusters consistent with annotated element
    //----------------------------------------------------------------------------
    //C_NNcluster2d c5c = cullRetA(contig, rclus.e5c, rclus.e5,  pars);
    //C_NNcluster2d c3c = cullRetA(contig, rclus.e3c, rclus.e3,  pars);
    
    //----------------------------------------------------------------------------
    // loop over retro clusters to identify candidate insertions
    //----------------------------------------------------------------------------
    // fiveprime first
    for ( i1=c5c.cluster.begin() ; i1 != c5c.cluster.end(); i1++ ) {
        int j = (*i1).first;
        C_SVR1 e1=makeRetEvent(contig, c5c.cluster[j], rclus.e5, Nexp  );
        ret.evt.push_back(e1);
    }
    // threeprime next
    for ( i1=c3c.cluster.begin() ; i1 != c3c.cluster.end(); i1++ ) {
        int j = (*i1).first;
        C_SVR1 e1=makeRetEvent(contig, c3c.cluster[j], rclus.e3, Nexp  );
        ret.evt.push_back(e1);
    }
    
    //----------------------------------------------------------------------------
    // sort to mix 5' and 3' in position order
    //----------------------------------------------------------------------------  
    ret.evt.sort();
    
    //----------------------------------------------------------------------------
    // merge adjacent events
    //----------------------------------------------------------------------------
    list<C_SVR1>::iterator ie,ie0,ie1,ie2;
    C_SVR1 e0,e1,ee;
    if (ret.evt.size()==0) { return 0;}
    // start with first event nearest to position 0
    e0=*ret.evt.begin();
    // bump ie0 pointer to second event
    ie0=ret.evt.begin();
    ie0++;
    // know when to stop
    ie1=ret.evt.end();
    
    
    //----------------------------------------------------------------------------
    // empty event list
    //----------------------------------------------------------------------------
    list<C_SVR1> evt;
    e0.type=ret.typeName;
    evt.push_back(e0);  // add first event to prime the event pump
    for ( ie=ie0 ; ie != ie1; ++ie ) {
        e1 = *ie;
        ie2=ie;
        ie2--;
        e0 = *ie2;
        if (abs( int(e1.pos)-int(e0.pos) )<LFmax) { 
            ee = merge(e0,e1,contig);     
            ee.type=ret.typeName;
            
            // remove existing event at end of list
            evt.pop_back(); 
            evt.push_back(ee);
        } else {
            e1.type=ret.typeName;
            evt.push_back(e1);
        }   
    }
    ret.evt=evt;    
    // sort list
    ret.evt.sort();
    
    // calc sample supporting and spanning fragments
    ret.finalize(contig,libraries,pars);
    
    // SVCF event count
    ret.SVCF.NEvent=int(ret.evt.size());
    
    return ret.SVCF.NEvent;
    
}       

//-------------------------------------------------------------------------
// make element insertion event from a cluster
//------------------------------------------------------------------------- 
C_SVR1  C_SpannerSV::makeRetEvent(C_contig  & contig,  
                                  C_cluster2d_element1 & c0, vector<C_umpair> & r0, double Nexp  ) {
    
    vector<C_umpair> retro1;  
    
    //-------------------------------------------------------------------------
    // Count of fragments
    //-------------------------------------------------------------------------
    int Np1 = c0.inp.size();
    // pmed1 for Deniz's median 
    list <unsigned int> pmed1;
    // constraint counts for Gabor...
    int Ncon=0;
    bool fiveprime=false;
    
    // ReadGroup info
    unsigned int RGCmax=0;    // ReadGroupCode for longest library
    int LFmax=0;              // LF width of longest library
    int LFrange=0;            // LF width of longest library
    C_RGmap RGmap;            // declare empty Readgroup counter
    unsigned int p5a=contig.Length;
    unsigned int p5b=0;   
    unsigned int p3a=contig.Length;  
    unsigned int p3b=0;  
    
    for (int j=0; j<Np1; j++) {
        int k = c0.inp[j];
        C_umpair r1 = r0[k];
        fiveprime=(r1.read[0].sense=='F');
        int pU;  // edge of U end closest to insertion
        unsigned int RGC1=r1.ReadGroupCode;
        RGmap[RGC1]++;
        int LM1=libraries.libmap[RGC1].LM;
        int LMhigh1=libraries.libmap[RGC1].LMhigh;
        int LMlow1=libraries.libmap[RGC1].LMlow;
        if (LMhigh1>LFmax) {
            RGCmax=RGC1;
            LFmax=LMhigh1;
            LFrange=LMhigh1-LMlow1;
        }
        
        unsigned int pUa=r1.read[0].pos;                   // low edge of U read
        unsigned int pUb=r1.read[0].pos+r1.read[0].len;    // high edge of U read
        
        if (fiveprime) {                                   // 5' cluster
            unsigned int pMa=r1.read[0].pos+LM1-r1.read[1].len;// low estimate M read
            unsigned int pMb=r1.read[0].pos+LM1;               // high estimate M read
            if (p5a>pUa) p5a = pUa;
            if (p5b<pUb) p5b = pUb;
            if (p3a>pMa) p3a = pMa;
            if (p3b<pMb) p3b = pMb;
            pU=pUb;
        } else {
            unsigned int pMa=r1.read[0].pos-LM1;                   // low estimate M read
            unsigned int pMb=r1.read[0].pos-LM1+r1.read[1].len;    // high estimate M read
            if (p5a>pMa) p5a = pMa;
            if (p5b<pMb) p5b = pMb;
            if (p3a>pUa) p3a = pUa;
            if (p3b<pUb) p3b = pUb;
            pU=pUa;
        }
        pmed1.push_back(pU);
        retro1.push_back(r1);
        
        if (r1.constrain(LMlow1,LMhigh1)) Ncon++;
    }
    
    //-------------------------------------------------------------------------
    // insertion breakpoint 
    //-------------------------------------------------------------------------
    int p0 = int(fiveprime? p5b: p3a );
    
    //-------------------------------------------------------------------------
    // new element insertion event
    //-------------------------------------------------------------------------
    C_SVR1 e1;
    e1.pos = p0;
    e1.anchor = contig.getAnchorIndex();
    // make length zero for insertion 
    e1.length = 0;
    //C_SVcoverage1 cov(contig, p0-100, p0+100,nomcov) ;
    //e1.cov=cov;  
    e1.p5[0]=p5a;
    e1.p5[1]=p5b;
    e1.p3[0]=p3a;
    e1.p3[1]=p3b;
    
    // average gap over this event
    int avegap = droundi(float((p5b-p5a)+(p3b-p3a))/float(2*Np1));
    // average gap over this contig
    //int avegapC = double(contig.Length)/double(contig.totalUniqueReads);
    // contig alignability estimate
    //float totSites = float(contig.Length-int(contig.totalNoCovBases));
    //float acontig  = float (totSites-contig.totalRepeatBases) / totSites;
    
    if (fiveprime) {
        e1.cls5=c0;
        e1.retro5=retro1;
        e1.NfragCluster[0]=Np1;
        // constraint counts
        e1.Nconstrain[0]=Ncon;
        e1.ReadGroupMap5=RGmap;    
        // coverage within U cluster region
        //C_SVcoverage1 covU(contig, p5a, p5b,nomcov) ;
        //e1.a5 = covU.Nsite/(float(int(covU.p1)-int(covU.p0))*acontig);  
    } else {
        e1.cls3=c0;
        e1.retro3=retro1;
        e1.NfragCluster[1]=Np1;
        // constraint counts
        e1.Nconstrain[1]=Ncon;
        e1.ReadGroupMap3=RGmap;    
        // coverage within U cluster region
        //C_SVcoverage1 covU(contig, p3a, p3b,nomcov) ;
        //e1.a3 = covU.Nsite/(float(int(covU.p1)-int(covU.p0))*acontig);  
    }
    // quantized copy number
    //e1.copynumber = short(2*e1.cov.N/e1.cov.eN +0.5);
    
    // brkpoint uncertainty depends on average gap between read starts in cluster
    //e1.posU=(avegap>avegapC? avegap: avegapC);
    e1.posU=avegap;
    //e1.lenU=e1.posU; // no more info
    // retro event info
    /*
     e1.NfragCov=int(contig.frag_depth.n[p0]);
     if (p0>LFmax) {
     e1.NfragCovOut[0]=int(contig.frag_depth.n[int(p0-LFmax)]);
     } else {
     e1.NfragCovOut[0]=0;
     }
     if ((p0+LFmax)<contig.Length) {
     e1.NfragCovOut[1]=int(contig.frag_depth.n[int(p0+LFmax)]);
     } else {
     e1.NfragCovOut[0]=0;
     }
     e1.NfragCovExp=int(Nexp);
     */
    //    
    // double nf = e1.NfragCov;
    // double pf = contig.frag_depth.Stats.h.x2pTrim(nf);
    // e1.q = p2q(pf);
    //
    // q value is set by supporting fragments (max 25)
    e1.q = int(25.0*Np1/(5.0+double(Np1)));    
    
    //-------------------------------------------------------------------------
    // median pos
    //-------------------------------------------------------------------------
    pmed1.sort();
    list <unsigned int>::iterator imed;
    int Nmed=pmed1.size()/2;
    int nmed=0;
    e1.pmedian =0;
    for (imed=pmed1.begin(); imed!=pmed1.end(); ++imed) {
        nmed++;
        if (nmed>Nmed) {
            e1.pmedian = *imed;
            break;
        }
    }
    return e1;
}       


//-------------------------------------------------------------------------
// cull element cluster 
//------------------------------------------------------------------------- 
C_NNcluster2d C_SpannerSV::cullRet(C_contig  & contig, C_NNcluster2d & clu1,  
                                   vector<C_umpair> & r1, RunControlParameters & pars) {
    
    C_NNcluster2d clu=clu1;
    C_cluster2d_elements::iterator it;
    C_cluster2d_element1 c1;
    //----------------------------------------------------------------------------
    // fragment length properties
    //----------------------------------------------------------------------------
    //double LMlow = pars.getFragmentLengthLo();
    //double LMhigh = pars.getFragmentLengthHi();
    //HistObj LFhist = pars.getFragHist();  
    //double LF= LFhist.median;
    int Nmin = pars.getMinClusteredRetro();
    
    //----------------------------------------------------------------------------
    // max allowed constrained fragments in cluster parameter (0=none, 10000=all)
    //----------------------------------------------------------------------------
    //int NmaxConstrain=1000000;
    int NmaxConstrain=1;
    
    //----------------------------------------------------------------------------
    // expected number of pairs that span a given deletion
    //----------------------------------------------------------------------------
    // double Nexp = contig.localpairs.size()*double(LF-2*contig.aLR)/double(contig.Length);
    //----------------------------------------------------------------------------
    // loop over retro clusters to identify candidate insertions
    //----------------------------------------------------------------------------
    for ( it=clu1.cluster.begin() ; it != clu1.cluster.end(); it++ ) {
        int i = (*it).first;
        c1 = clu.cluster[i];   
        //-------------------------------------------------------------------------
        // highest start position mark start of deletion 
        //-------------------------------------------------------------------------
        // int p0 = int(c1.high+round(contig.aLR)+1);
        int Np = c1.inp.size();
        //-------------------------------------------------------------------------
        // initialize cluster limits:    p5a--p5b      p3a-p3b
        //-------------------------------------------------------------------------
        int p5b = 0;              // trailing edge of 5' cluster
        int p5a = contig.Length;  // leading edge of 5' cluster
        int Nconstrain=0;
        if (Np<Nmin) { Np=0; }    // bug out if Np too small
        //
        // retro.resize(Np);
        for (int j=0; j<Np; j++) {
            int k = c1.inp[j];
            C_umpair r = r1[k];
            
            int p = r.read[0].pos;                // start of 5' cluster
            if (p<p5a) p5a=p;
            p = r.read[0].pos+r.read[0].len;      // end of 5' cluster
            if (p>p5b) p5b=p;
            
            bool constrain=false;
            double lm=1e20;
            if (r.read[0].anchor==r.read[1].anchor) {
                lm = r.read[1].pos+r.read[1].len-r.read[0].pos;
                if (r.read[0].sense=='R') {
                    lm = r.read[0].pos+r.read[0].len-r.read[1].pos;
                }
                unsigned int RGC1=r.ReadGroupCode;
                int LMhigh1=libraries.libmap[RGC1].LMhigh;
                int LMlow1=libraries.libmap[RGC1].LMlow;
                constrain=(lm<LMhigh1)&&(lm>LMlow1);
                if (constrain) { Nconstrain++; }
            }       
        }
        // remove any cluster with a constrained fragment
        if ( (Nconstrain>=NmaxConstrain)||(Np<Nmin) ){
            clu.cluster.erase(i);
        }
    }
    return clu;  
}

//------------------------------------------------------------------------------
// merge adjacent retro clusters
//------------------------------------------------------------------------------
C_SVR1 C_SpannerSV::merge(C_SVR1 & e0, C_SVR1 & e1, C_contig  & contig) {
    
    
    //----------------------------------------------------------------------------
    // expected number of spanning frags
    //----------------------------------------------------------------------------
    //double Nexp = contig.frag_depth.Stats.h.mean;
    
    //----------------------------------------------------------------------------
    // flip e0 to F ,  and e1 to R
    //----------------------------------------------------------------------------
    if (e0.NfragCluster[0]==0)  {
        // put F cluster in e0, R in e1
        C_SVR1 etemp=e0;
        e0=e1;
        e1=etemp;
    } 
    
    //----------------------------------------------------------------------------
    // declare merged event, start with e0
    //----------------------------------------------------------------------------
    C_SVR1 em;
    em=e0;   
    
    //----------------------------------------------------------------------------
    // check for clustered F's or R's
    //----------------------------------------------------------------------------
    bool bothF= (e0.NfragCluster[0]>0) && (e1.NfragCluster[0]>0) ;
    bool bothR= (e0.NfragCluster[1]>0) && (e1.NfragCluster[1]>0) ;  
    if (bothF) {     
        em.p5[0] = (e0.p5[0]<e1.p5[0]? e0.p5[0]: e1.p5[0]);
        em.p5[1] = (e0.p5[1]>e1.p5[1]? e0.p5[1]: e1.p5[1]);
        em.retro5 = e0.retro5;
        em.retro5.insert(em.retro5.end(), e1.retro5.begin(), e1.retro5.end() );
        em.Nconstrain[0]=e0.Nconstrain[0]+e1.Nconstrain[0];
        em.merge5=true;
    } else if (bothR) {
        em.p3[0] = (e0.p3[0]<e1.p3[0]? e0.p3[0]: e1.p3[0]);
        em.p3[1] = (e0.p3[1]>e1.p3[0]? e0.p3[1]: e1.p3[1]);
        em.retro3 = e0.retro3;
        em.retro3.insert(em.retro3.end(), e1.retro3.begin(), e1.retro3.end() );
        em.Nconstrain[1]=e0.Nconstrain[1]+e1.Nconstrain[1];
        em.merge3=true;
    } else {
        em.p3[0] = e1.p3[0];
        em.p3[1] = e1.p3[1];
        // ee.p5 already set above 
        em.cls3=e1.cls3;
        em.retro5.insert(em.retro5.end(), e1.retro5.begin(), e1.retro5.end() );
        em.retro3.insert(em.retro3.end(), e1.retro3.begin(), e1.retro3.end() );
        em.Nconstrain[1]=e1.Nconstrain[1];
        em.merge3=e1.merge3;
    }
    
    //----------------------------------------------------------------------------
    // position of event
    //----------------------------------------------------------------------------
    int N0=e0.NfragCluster[0]+e0.NfragCluster[1];
    int N1=e1.NfragCluster[0]+e1.NfragCluster[1];
    int Np=N0+N1;
    // pos estimate is p5b if F cluster, p3a if R
    double P0=(e0.NfragCluster[0]>0? double(e0.p5[1]) : double(e0.p3[0]));
    double P1=(e1.NfragCluster[1]>0? double(e1.p3[0]) : double(e1.p5[1]));
    
    //em.pos = int((P0*N0+P1*N1)/double(Np));  
	// leftmost insertion pos (before tsd)	
    
	em.pos = int(P0<P1? P0: P1);
    //em.length = int(P0)-int(P1); // should be ~0 unless insert is shorter than LF or repeat region
    //em.gap = int(P0)-int(P1); // should be ~0 unless insert is shorter than LF or repeat region
    em.gap = droundi(P1-P0); //em.p3[0]-em.p5[1];
    
    // average gap over this event
	em.posU= droundi(float(em.p5[1]-em.p5[0])/float(N0));
	em.lenU= droundi(float(em.p3[1]-em.p3[0])/float(N1));
    // average gap over this contig
    // int avegapC = double(contig.Length)/int(2*contig.localpairs.size());
    //int avegapC = double(contig.Length)/double(contig.totalUniqueReads);
    // brkpoint uncertainty depends on average gap between read starts in cluster
    //em.posU=avegap; //(avegap>avegapC? 2*avegap: 2*avegapC);
    //em.lenU=avegap; // meaningless 
    
    em.NfragCluster[0] = e0.retro5.size()+e1.retro5.size();
    em.NfragCluster[1] = e0.retro3.size()+e1.retro3.size();
    int p0 = (int(em.pos)>100? int(em.pos)-100: 0);  
    //int p1 = (int(em.pos)<(contig.Length-101)? em.pos+100: contig.Length-1);
    //C_SVcoverage1 cov(contig, p0, p1,nomcov);
    //em.cov=cov; 
    
    //----------------------------------------------------------------------------
    // alignability estimates for 3' end of F (a5) or R cluster (a3) 
    //----------------------------------------------------------------------------
    em.a5=e0.a5;
    em.a3=e1.a3;
    
    //----------------------------------------------------------------------------
    // ReadGroups from e1
    //----------------------------------------------------------------------------
    std::map<unsigned int,unsigned int, less<unsigned int> >::iterator it;
    for ( it=e1.ReadGroupMap3.begin() ; it != e1.ReadGroupMap3.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        unsigned int NF1  = e1.ReadGroupMap3[RGC1];
        em.ReadGroupMap3[RGC1]+=NF1;
    }
    for ( it=e1.ReadGroupMap5.begin() ; it != e1.ReadGroupMap5.end(); it++ ) {
        unsigned int RGC1 = (*it).first;
        unsigned int NF1  = e1.ReadGroupMap5[RGC1];
        em.ReadGroupMap5[RGC1]+=NF1;
    }
    
    int LFmax=0;
    //----------------------------------------------------------------------------
    // clustering (+ silly pmedian per Deniz) 
    //----------------------------------------------------------------------------
    list <unsigned int> pmed1;
    
    if (bothF||bothR) {
        C_cluster2d_element1 c1;
        c1.N=em.retro5.size();
        if (bothR) {c1.N=em.retro3.size();}
        c1.low[0]=contig.Length;
        for (int i=0; i<c1.N; i++)  {
            double x=0;
            C_umpair p1;
            if (bothR) {
                p1=em.retro3[i];
            } else {
                p1=em.retro5[i];
            }  
            int LF1 = libraries.libmap[p1.ReadGroupCode].LM;
            if (LF1>LFmax) LFmax=LF1;
            if (bothR) {
                x= p1.read[0].pos+LF1;
            } else {
                x= p1.read[0].pos+p1.read[0].len-LF1;
            }
            c1.mean[0]+=x;
            c1.std[0]+=x*x;
            c1.low[0]=(x>c1.low[0] ? c1.low[0]: x);
            c1.high[0]=(x<c1.high[0] ? c1.high[0]: x);
            c1.inp.push_back(i);        
            pmed1.push_back(int(x));
        }
        c1.mean[0]=c1.mean[0]/c1.N;
        c1.std[0]=sqrt(c1.std[0]/c1.N-(c1.mean[0]*c1.mean[0]));
        if (bothF) {
            em.cls5=c1;
        }  else {
            em.cls3=c1;
        }   
    } else {
        int Nfrag=em.retro5.size();
        for (int i=0; i<Nfrag; i++)  {
            int LF1 = libraries.libmap[em.retro5[i].ReadGroupCode].LM;
            if (LF1>LFmax) LFmax=LF1;
            double x = em.retro5[i].read[0].pos+em.retro5[i].read[0].len-LF1;
            pmed1.push_back(int(x));
        }
        Nfrag=em.retro3.size();
        for (int i=0; i<Nfrag; i++)  {
            int LF1 = libraries.libmap[em.retro3[i].ReadGroupCode].LM;
            if (LF1>LFmax) LFmax=LF1;
            double x = em.retro3[i].read[0].pos+LF1;
            pmed1.push_back(int(x));
        }
    }
    
    //----------------------------------------------------------------------------
    // local copy number
    //----------------------------------------------------------------------------
    //em.copynumber = short(2*em.cov.N/em.cov.eN +0.5);  
    
    //----------------------------------------------------------------------------
    // fragment coverage at and around insertion
    //----------------------------------------------------------------------------
    p0=em.pos;
	/*
     em.NfragCov=int(contig.frag_depth.n[p0]);
     if (p0>LFmax) {
     em.NfragCovOut[0]=int(contig.frag_depth.n[p0-int(LFmax)]);
     } else {
     em.NfragCovOut[0]=0;
     }
     if ((p0+LFmax)<contig.Length) {
     em.NfragCovOut[1]=int(contig.frag_depth.n[p0+int(LFmax)]);
     } else {
     em.NfragCovOut[0]=0;
     }
     em.NfragCovExp=int(Nexp);
     */
	
    // q value is set by NfragCov
    //double nf = em.NfragCov;
    //double pf = contig.frag_depth.Stats.h.x2pTrim(nf);
    //em.q = p2q(pf);
    
    // q value is set by supporting fragments 
    em.q = int(50.0*double(Np)/(10.0+double(Np)));    
    
    //-------------------------------------------------------------------------
    // median pos
    //-------------------------------------------------------------------------
    pmed1.sort();
    list <unsigned int>::iterator imed;
    int Nmed=pmed1.size()/2;
    int nmed=0;
    em.pmedian =0;
    for (imed=pmed1.begin(); imed!=pmed1.end(); ++imed) {
        nmed++;
        if (nmed>Nmed) {
            em.pmedian = *imed;
            break;
        }
    }
    C_SVR1 eout;
    eout=em;
    //eout.lenU=1;  // marker
    return eout;
}


char C_SpannerSV::p2q(double p) {
    double Q = -10.0*log10(p); 
    char q = char(Q);
    // overflow?  make 40 higest q - 
    double Qmax=40;
    if ((Q>Qmax)||(p<1e-4)) q=char(Qmax);
    if ((q==0)&&(p<1e-4)) q=char(Qmax);
    return q;
}

C_SVspanfrags1::C_SVspanfrags1() {
    N=0;
    NN=0;
    N5=0;
    N3=0;
    NR=0;
    ER=0.0;
}

C_SVspanfrags1::C_SVspanfrags1(unsigned int nn, unsigned int n5, unsigned int n3,unsigned int nr, float er) {
    N=n3+n5;
    NN=nn;
    N5=n5;
    N3=n3;
    NR=nr;
    ER=er;
}

C_SVspanfrags1& C_SVspanfrags1::operator=(const C_SVspanfrags1 &rhs)
{
    this->N = rhs.N;
    this->NN = rhs.NN;
    this->N5 = rhs.N5;
    this->N3 = rhs.N3;
    this->NR = rhs.NR;
    this->ER = rhs.ER;
    return  *this;
}

int C_SVspanfrags1::operator==(const C_SVspanfrags1 &rhs) const
{
    if( this->N != rhs.N) return 0;
    if( this->NN != rhs.NN) return 0;
    if( this->N5 != rhs.N5) return 0;
    if( this->N3 != rhs.N3) return 0;
    if( this->NR != rhs.NR) return 0;
    if( this->ER != rhs.ER) return 0;
    return 1;
}

C_SVCF_TAG::C_SVCF_TAG() {
	Id="";
	Number=0;
	Type="";  // Integer, Float, String, Flag
	Descr="";
	INFO=false;
	ALT=false;
	FORMAT=false;
	FILTER=false;
}

C_SVCF_TAG::C_SVCF_TAG(char I1, string & Id1, int Number1, string & Type1,string & Descr1 ) {
	
	switch (I1) {
		case 0:
			INFO=true;
			break;
		case 1:
			ALT=true;
			break;
		case 2:
			FORMAT=true;
			break;
		case 3:
			FILTER=true;
			break;
		default:
			break;
	}	
	Id=Id1;
	Number=Number1;
	Type=Type1;  // Integer, Float, String, Flag
	Descr=Descr1;
}


ostream &operator<<(ostream &output, C_SVCF_TAG & tag)
{
    //----------------------------------------------------------------------------
    //  print SVCF format header 
    //----------------------------------------------------------------------------  
    char b [200];
    if (tag.ALT) { 		
		sprintf(b,"##ALT=<ID=%s,Description=\"%s\">", tag.Id.c_str(),tag.Descr.c_str());
	} else if (tag.INFO) { 
		sprintf(b,"##INFO=<ID=%s,Number=%d,Type=%s,Description=\"%s\">", 
				tag.Id.c_str(),tag.Number,tag.Type.c_str(),tag.Descr.c_str());
	} else if (tag.FORMAT) { 
		sprintf(b,"##FORMAT=<ID=%s,Number=%d,Type=%s,Description=\"%s\">", 
				tag.Id.c_str(),tag.Number,tag.Type.c_str(),tag.Descr.c_str());
	} else if (tag.FILTER) { 
		sprintf(b,"##FILTER=<ID=%s,Description=\"%s\">", tag.Id.c_str(),tag.Descr.c_str());
	}
	string s = b;
	output << s ;
	return output;	
}	


C_SVCF::C_SVCF() {
    Version="VCFv4.0";
    // date
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (buffer,80,"%Y%m%d",timeinfo);
    Date=buffer;
    // source
    Source="Spanner:V";
    Source=Source+SpannerSvnVersion;
    // blanks
    Reference="";
    EventType="";
    NEvent=0;
    ALT.clear();
	INFO.clear();
	FMT.clear();
	FILT.clear();
	//Info.clear();
    //Format.clear();
    Samples.clear();  
}


C_SVCF::C_SVCF(RunControlParameters & pars, vector<C_SVCF_TAG> & Info1,vector<C_SVCF_TAG> & Format1, vector<C_SVCF_TAG> & Filt1, 
               vector<C_SVCF_TAG> & Alt1,vector<string> & Samples1) {
    // fill-in-the-blanks
    Reference=pars.getAnchorfile();
    //Info=Info1;
    //Format=Format1;
    Samples=Samples1;
}


ostream &operator<<(ostream &output, C_SVCF & svcf1)
{
    //----------------------------------------------------------------------------
    //  print SVCF format header 
    //----------------------------------------------------------------------------  
    char b [200];
    sprintf(b,"##format=%s", svcf1.Version.c_str());
    string s = b;
    output << s << endl;
	
    
    sprintf(b,"##fileDate=%s", svcf1.Date.c_str());
    s = b;
    output << s << endl;
    sprintf(b,"##source=%s", svcf1.Source.c_str());
    s = b;
    output << s << endl;
    if (svcf1.Reference.size()>1) {
        sprintf(b,"##reference=%s", svcf1.Reference.c_str());
        s = b;
        output << s << endl;
    }
    
    if (svcf1.EventType.size()>1) {
        sprintf(b,"##event=%s", svcf1.EventType.c_str());
        s = b;
        output << s << endl;
        sprintf(b,"##Nevt=%d", svcf1.NEvent);
        s = b;
        output << s << endl;
    }
    
	for (size_t i=0; i<svcf1.ALT.size(); i++) {
		output << svcf1.ALT[i] << endl;
	}
	for (size_t i=0; i<svcf1.INFO.size(); i++) {
		output << svcf1.INFO[i] << endl;
	}
	for (size_t i=0; i<svcf1.FMT.size(); i++) {
		output << svcf1.FMT[i] << endl;
	}
	for (size_t i=0; i<svcf1.FILT.size(); i++) {
		output << svcf1.FILT[i] << endl;
	}
    
    output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    
    size_t NS=svcf1.Samples.size();
    if (NS>0) {
        for (int i=0; i<int(NS); i++) {
            sprintf(b,"\t%s", svcf1.Samples[i].c_str());
            s = b;
            output << s;
        }
    }
    output << endl;
    return output;
    
}




C_SV1::C_SV1() {
    type='d';
    pos=0;
    anchor=0;
    length=0;
  	pend=0;
    q=0;
    p5[0]=0;
    p5[1]=0;
    p3[0]=0;
    p3[1]=0;
    copynumber=0;
  	CIpos[0]=0;
    CIpos[1]=0;
    CIend[0]=0;
    CIend[1]=0;
    CIlen[0]=0;
    CIlen[1]=0;
    qOutlier=0;
    qAberrantLM=0;
    merge=false;
    a5=0;
    a3=0;   
    q5=0;
    q3=0;   
    SampleMap.clear();
    id=0;
}

ostream &operator<<(ostream &output, C_SV1 & e1)
{
    //----------------------------------------------------------------------------
    //  text file position coordinates are 1 based (not c convention)
    //----------------------------------------------------------------------------  
    char b [200];
    string s=e1.type;
	if (s.compare("TDUP")==0) {
		s="DUP:TANDEM";
	}
    sprintf(b,"%d\t%d\t%d\t.\t<%s>\t%d\tPASS\t", e1.anchor, e1.pos+1,e1.id, s.c_str(),e1.q);
    s = b;
    output << s;
    //int u=e1.posU; 
    sprintf(b,"SVLEN=%d;CIPOS=%d,%d;NALTF=%d;", e1.length,e1.CIpos[0],e1.CIpos[1],int(e1.cls.inp.size()));
    s = b;
    output << s;
    
    sprintf(b,"END=%d;CIEND=%d,%d;", e1.pend,e1.CIend[0],e1.CIend[1]);
    s = b;
    output << s;
	
	sprintf(b,"CISVLEN=%d,%d;", e1.CIlen[0],e1.CIlen[1]);
    s = b;
    output << s;
    
	sprintf(b,"QF=%d;QC=%d;", e1.qAberrantLM, e1.qOutlier);
    s = b;
    output << s;
	
    int PA=e1.p5[0]-e1.pos;
    int PB=e1.p5[1]-e1.pos;
    int PC=e1.p3[0]-e1.pos;
    int PD=e1.p3[1]-e1.pos;
    
    sprintf(b,"PC1=%d,%d;PC2=%d,%d;",PA,PB,PC,PD);
    s = b;
    output << s;
    
	sprintf(b,"MQ5=%d;MQ3=%d",e1.q5,e1.q3);
    s = b;
    output << s;
	
	/*
     sprintf(b,"NR=%d;ER=%.1f;MR=%d;", e1.cov.N,e1.cov.eN,(e1.merge)?1:0);
     s = b;
     output << s;  
     
     sprintf(b,"A5=%.2f;A3=%.2f;", e1.a5,e1.a3);
     s = b;
     output << s;  
     */
	
    return output;
}




C_SV1& C_SV1::operator=(const C_SV1 &rhs)
{
    this->pos = rhs.pos;
    this->anchor = rhs.anchor;
    this->length = rhs.length;
    this->q = rhs.q;
    this->p5[0] = rhs.p5[0];
    this->p5[1] = rhs.p5[1];
    this->p3[0] = rhs.p3[0];
    this->p3[1] = rhs.p3[1];
    this->copynumber = rhs.copynumber;
    this->cov   = rhs.cov;
    this->cov5  = rhs.cov5;
    this->cov3  = rhs.cov3;
    this->cls   = rhs.cls;  
    this->pair  = rhs.pair;  
    this->pend = rhs.pend;
    this->CIpos[0]  = rhs.CIpos[0];
    this->CIpos[1]  = rhs.CIpos[1];
    this->CIend[0]  = rhs.CIend[0];
    this->CIend[1]  = rhs.CIend[1];
    this->CIlen[0]  = rhs.CIlen[0];
    this->CIlen[1]  = rhs.CIlen[1];
    this->qOutlier   = rhs.qOutlier;
    this->qAberrantLM= rhs.qAberrantLM;
    this->ReadGroupMap=rhs.ReadGroupMap;
    this->merge     = rhs.merge;
    this->a5        = rhs.a5;
    this->a3        = rhs.a3;
    this->q5        = rhs.q5;
    this->q3        = rhs.q3;
    this->type      = rhs.type;
    this->id        = rhs.id;
    this->SampleMap = rhs.SampleMap;
    
    return *this;
}

int C_SV1::operator==(const C_SV1 &rhs) const
{
	if( this->pos != rhs.pos) return 0;
	if( this->anchor != rhs.anchor) return 0;
	if( this->length != rhs.length) return 0;
	if( this->p5[0] != rhs.p5[0]) return 0;
	if( this->p5[1] != rhs.p5[1]) return 0;
	if( this->p3[0] != rhs.p3[0]) return 0;
	if( this->p3[1] != rhs.p3[1]) return 0;
	if( this->q != rhs.q) return 0;
	if( this->pend != rhs.pend) return 0;
	if( this->CIpos[0] != rhs.CIpos[0]) return 0;
	if( this->CIpos[1] != rhs.CIpos[1]) return 0;
	if( this->CIend[0] != rhs.CIend[0]) return 0;
	if( this->CIend[1] != rhs.CIend[1]) return 0;
	if( this->CIlen[0] != rhs.CIlen[0]) return 0;
	if( this->CIlen[1] != rhs.CIlen[1]) return 0;
	if( this->qOutlier != rhs.qOutlier) return 0;
	if( this->qAberrantLM != rhs.qAberrantLM) return 0;
	if (this->copynumber != rhs.copynumber) return 0;
	if (this->merge != rhs.merge) return 0;
	if (this->ReadGroupMap.size() != rhs.ReadGroupMap.size()) return 0;
	if( this->a5 != rhs.a5) return 0;
	if( this->a3 != rhs.a3) return 0;
	if( this->q5 != rhs.q5) return 0;
	if( this->q3 != rhs.q3) return 0;
	return 1;
}

int C_SV1::operator<(const C_SV1 &rhs) const
{
    if( this->anchor < rhs.anchor ) return 1;
    if( this->anchor == rhs.anchor && this->pos < rhs.pos ) return 1;
    if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->length < rhs.length ) return 1;
    if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->length == rhs.length  && this->q < rhs.q) return 1;
    return 0;
}


C_SV::C_SV() {
    typeName="x";
    contigName="";
    setName="";
    evt.clear(); 
}

C_SV::C_SV(C_contig  & c1,  RunControlParameters & par) {
    evt.clear(); 
    typeName="x";
    contigName =c1.getContigName();
    setName =c1.setName;
}

void C_SV::finalize(C_contig  & contig, C_libraries & libraries, RunControlParameters & pars) {
    
    //----------------------------------------------------------------------------
    // sample by sample genotype info
    //----------------------------------------------------------------------------
    genotype(contig,libraries,pars);
    
    //----------------------------------------------------------------------------
    // loop over events / clusters within events to calculate:
    //----------------------------------------------------------------------------
    
    // check for events, libraries
    if (evt.size()==0) { return;}
    if (libraries.libmap.size()==0) { return;}
    
    list<C_SV1>::iterator ie,ie1,ie2;
    // first event
    ie1=evt.begin();
    // last event
    ie2=evt.end();  
    // loop over events
    int ne=0;
    
    // max a5 or a3 alignability metrics
    double aMax=0;
    
    for ( ie=ie1 ; ie != ie2; ++ie ) {
        ne++;
        
        // id
        (*ie).id=ne;
        
        // find max alignability scale
        if ((*ie).a5>aMax) { aMax=(*ie).a5;}
        if ((*ie).a3>aMax) { aMax=(*ie).a3;}
        
    }
    
    // normalize alignability to <=1.0
    for ( ie=evt.begin(); ie != evt.end(); ++ie ) {
        (*ie).a5=(*ie).a5/aMax;
        (*ie).a3=(*ie).a3/aMax;
    }
    
} 




void C_SV::genotype(C_contig  & contig, C_libraries & libraries, RunControlParameters & pars)
{
    // check for events, libraries
    if (evt.size()==0) { return;}
    if (libraries.libmap.size()==0) { return;}
    
    map<string, double, less<string> >  rX = libraries.readFractionSamples();
	
    list<C_SV1>::iterator ie,ie1,ie2;
    //C_SVR1 e1;
    unsigned int ReadGroupCode1;
    string SAM;
    int NF1;
    
    ie1=evt.begin();
	
    // know when to stop
    ie2=evt.end();  
    
    list<C_localpair>::iterator i;
	
	/*
     // get Minimum Q mapping value for unique read
     int Qmin  = pars.getQmin();
     
     // loop over pairs
     int np=0;
     for(i=contig.localpairs.begin(); i != contig.localpairs.end(); ++i) {
     np++;
     // skip low mapping quality fragments
     if (((*i).q1<Qmin)&((*i).q2<Qmin)) {
     continue;
     }
     // fragment length limits for this library
     ReadGroupCode1 = (*i).ReadGroupCode;  
     //int LM1=libraries.libmap[ReadGroupCode1].LM;
     int LMhigh1=libraries.libmap[ReadGroupCode1].LMhigh;
     int LMlow1=libraries.libmap[ReadGroupCode1].LMlow;
     
     // fragment length
     int lm = (*i).lm;
     // demand usual fragment
     if ((lm<LMlow1)|(lm>LMhigh1)) {continue; }
     
     SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
     
     // add only the non-read part of the fragment
     int p0 = (*i).pos;
     // end of 5' read
     int p1 = (*i).pos+(*i).len1;
     //int p2 = p1+(*i).lm-(*i).len2;
     // begining of 3' read at p0+lm-len2
     int p2 = p0+(*i).lm-(*i).len2;
     // end of 3' read
     int p3 = p0+(*i).lm;
     
     // loop over events
     int ne=0;
     for ( ie=ie1 ; ie != ie2; ++ie ) {
     
     ne++;
     
     // two spanning windows eP1a-eP1b and eP2b-eP2b  
     int eP1a=(*ie).pos-(*ie).posU/2;     
     int eP1b=int((*ie).pos)+int((*ie).posU/2);
     int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
     int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
     
     // check for reads starts inside event window eP1b-eP2a
     if ((p0>=eP1b)&(p0<=eP2a)&((*i).q1>=Qmin)) {       
     (*ie).SampleMap[SAM].NR++;
     }
     // other end of fragment
     if ((p2>=eP1b)&(p2<=eP2a)&((*i).q1>=Qmin)) {       
     (*ie).SampleMap[SAM].NR++;
     }
     
     // beginning of event well beyond this fragment
     if ((eP1a-1000)>p3) { break;}        
     
     // end of event not yet reaching this fragment
     if (eP2b<p1) { continue;}       
     
     // select spanning fragments - both above Qmin 
     if (((*i).q1<Qmin)|((*i).q2<Qmin)) {
     continue;
     } 
     
     // spanning  breakpoint
     if ((p1<eP1a)&(p2>eP1b)) {
     (*ie).SampleMap[SAM].NN++;
     }
     
     }
     
     // start next loop over events where this one left off - 3
     if (ie!=ie1)  ie1=ie--;
     if (ie1!=evt.begin())  ie1--;
     if (ie1!=evt.begin())  ie1--;
     
     }
     
     ie1=evt.begin();
     */
	
    std::map<unsigned int,unsigned int, less<unsigned int> >::iterator it;  
    // iterate over read groups in 5' end;
    for ( ie=ie1 ; ie != ie2; ++ie ) {
        //e1 = *ie;
        for ( it=(*ie).ReadGroupMap.begin() ; it != (*ie).ReadGroupMap.end(); it++ ) {
			ReadGroupCode1 = (*it).first;
			NF1  = (*ie).ReadGroupMap[ReadGroupCode1];
			SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
			//(*ie).SampleMap[SAM].N5+=NF1;
			(*ie).SampleMap[SAM].N+=NF1;
        }     
        
		/*
         // loop over libraries for estimated number of reads 
         //map<string, double, less<string> > rX = libraries.readFractionSamples();
         map<string, double, less<string> >::iterator irX;
         for ( irX=rX.begin() ; irX != rX.end(); irX++ ) {
         SAM=(*irX).first;
         (*ie).SampleMap[SAM].ER=double((*ie).cov.eN)*rX[SAM];
         }
         */
    }
	
    /*
     // fill depth from cross pairs (not so many...)
     int N = contig.crosspairs.size();
     list<C_crosspair>::iterator ix;
     ie1=evt.begin();
     for(ix=contig.crosspairs.begin(); ix != contig.crosspairs.end(); ++ix) {
     // skip low mapping quality fragments
     if ((*ix).read[0].q<Qmin) {continue;}
     // fragment length limits for this library
     ReadGroupCode1 = (*ix).ReadGroupCode;  
     SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
     
     // read start
     int p1=int((*ix).read[0].pos);
     
     // loop over events
     int nx=0;
     for ( ie=ie1 ; ie != ie2; ++ie ) {
     nx++;
     
     int eP1a=(*ie).pos-(*ie).posU/2;     
     int eP1b=int((*ie).pos)+int((*ie).posU/2);
     int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
     int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
     
     // beginning of event well beyond this fragment
     if ((eP1a-1000)>p1) { break;}        
     
     // end of event not yet reaching this fragment
     if (eP2b<p1) { continue;}       
     
     // check for reads starts inside event window eP1b-eP2a
     if ((p1>=eP1b)&(p1<=eP2a)&((*ix).read[0].q>=Qmin)) {       
     (*ie).SampleMap[SAM].NR++;
     }
     
     }
     
     // start next loop over events where this one left off - 3
     if (ie!=ie1)  ie1=ie--;
     if (ie1!=evt.begin())  ie1--;
     if (ie1!=evt.begin())  ie1--;
     }
     
     // fill depth from dangle starts 
     N = contig.dangle.size();
     list<C_singleEnd>::iterator is; 
     for(is=contig.dangle.begin(); is != contig.dangle.end(); ++is) {
     // skip low mapping quality fragments
     if ((*is).q<Qmin) {continue;}
     ReadGroupCode1 = (*is).ReadGroupCode;  
     SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
     
     // read start
     int p1=int((*is).pos);
     
     // loop over events
     int nd=0;
     for ( ie=ie1 ; ie != ie2; ++ie ) {
     nd++;
     
     int eP1a=(*ie).pos-(*ie).posU/2;     
     int eP1b=int((*ie).pos)+int((*ie).posU/2);
     int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
     int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
     
     // beginning of event well beyond this fragment
     if ((eP1a-1000)>p1) { break;}        
     
     // end of event not yet reaching this fragment
     if (eP2b<p1) { continue;}       
     
     // check for reads starts inside event window eP1b-eP2a
     if ((p1>=eP1b)&(p1<=eP2a)&((*is).q>=Qmin)) {       
     (*ie).SampleMap[SAM].NR++;
     }
     
     }
     
     // start next loop over events where this one left off - 3
     if (ie!=ie1)  ie1=ie--;
     if (ie1!=evt.begin())  ie1--;
     if (ie1!=evt.begin())  ie1--;
     }
     // fill depth from umpairs 
     N = contig.umpairs.size();
     list<C_umpair>::iterator iu;
     
     for(iu=contig.umpairs.begin(); iu != contig.umpairs.end(); ++iu) {
     if ((*iu).read[0].q<Qmin) {continue;}
     ReadGroupCode1 = (*iu).ReadGroupCode;  
     SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
     
     // read start
     int p1=int((*iu).read[0].pos);
     
     // loop over events
     int nd=0;
     for ( ie=ie1 ; ie != ie2; ++ie ) {
     nd++;
     
     int eP1a=(*ie).pos-(*ie).posU/2;     
     int eP1b=int((*ie).pos)+int((*ie).posU/2);
     int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
     int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
     
     // beginning of event well beyond this fragment
     if ((eP1a-1000)>p1) { break;}        
     
     // end of event not yet reaching this fragment
     if (eP2b<p1) { continue;}       
     
     // check for reads starts inside event window eP1b-eP2a
     if ((p1>=eP1b)&(p1<=eP2a)&((*iu).read[0].q>=Qmin)) {       
     (*ie).SampleMap[SAM].NR++;
     }
     
     }
     
     // start next loop over events where this one left off - 3
     if (ie!=ie1)  ie1=ie--;
     if (ie1!=evt.begin())  ie1--;
     if (ie1!=evt.begin())  ie1--;
     }
     
	 
     // fill depth from unique ends 
     N = contig.singleton.size();
     for(is=contig.singleton.begin(); is != contig.singleton.end(); ++is) {
     if ((*is).q<Qmin) {continue;}
     ReadGroupCode1 = (*is).ReadGroupCode;  
     SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
     
     // read start
     int p1=int((*is).pos);
     
     // loop over events
     int nd=0;
     for ( ie=ie1 ; ie != ie2; ++ie ) {
     nd++;
     int eP1a=(*ie).pos-(*ie).posU/2;     
     int eP1b=int((*ie).pos)+int((*ie).posU/2);
     int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
     int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
     
     // beginning of event well beyond this fragment
     if ((eP1a-1000)>p1) { break;}        
     
     // end of event not yet reaching this fragment
     if (eP2b<p1) { continue;}       
     
     // check for reads starts inside event window eP1b-eP2a
     if ((p1>=eP1b)&(p1<=eP2a)&((*is).q>=Qmin)) {       
     (*ie).SampleMap[SAM].NR++;
     }
     }
     // start next loop over events where this one left off - 3
     if (ie!=ie1)  ie1=ie--;
     if (ie1!=evt.begin())  ie1--;
     if (ie1!=evt.begin())  ie1--;
     }
     */
    
}

ostream &operator<<(ostream &output, C_SV & sv)
{
	
    // SVCF header
    output  << sv.SVCF;
	
    size_t NS = sv.samples.size();
	
    char b[200];
    string s;
    
    // format 
    size_t NFMT=sv.SVCF.FMT.size();
    string FMT=sv.SVCF.FMT[0].Id;
    for (int f=1; f<int(NFMT); f++) {
		FMT=FMT+":"+sv.SVCF.FMT[f].Id;
    }
    // events
    list<C_SV1>::iterator i;
    for(i=sv.evt.begin(); i != sv.evt.end(); ++i) {
		output << (*i);
		
		output << "\t" << FMT ;
		
		for (int ns=0; ns<int(NS); ns++) {
			s=sv.samples[ns];
			// double cn = 2.0*double((*i).SampleMap[s].NR)/((*i).SampleMap[s].ER+0.01);
			// sprintf(b,"\t%d:%d:%d:%d:%d:%.1f",(*i).SampleMap[s].N,(*i).SampleMap[s].NN,(*i).SampleMap[s].N5,(*i).SampleMap[s].N3,(*i).SampleMap[s].NR,(*i).SampleMap[s].ER); //cn);
			//sprintf(b,"\t%d:%d:%d:%.1f",(*i).SampleMap[s].N,(*i).SampleMap[s].NN,(*i).SampleMap[s].NR,(*i).SampleMap[s].ER); //cn);
			sprintf(b,"\t%d",(*i).SampleMap[s].N); //cn);
			s = b;
			output << s;
		}
		output  << endl;  
    }  
    return output;
}


/*
 ostream &operator<<(ostream &output, C_SV & sv)
 {
 output << sv.typeName << "\t " << sv.setName << "\t " << sv.contigName << endl;
 int Nevt = sv.evt.size();
 output << Nevt << endl;
 char b[200];
 sprintf(b,"%10s\t%10s\t%4s\t%4s", "pos","length","anc","q");
 string s = b;
 output << s;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s", "Nread","p0","p1","Nsite","Nrepeat","eN","p");
 s = b;
 output << s;
 sprintf(b,"\t%10s", "Nf");
 s = b;
 output << s;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s", "<p>","Pstd","Plow","Phigh");
 s = b;
 output << s;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s", "<lm>","Lstd","Llow","Lhigh");
 s = b;
 output << s;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s", "p5a","p5b","p3a","p3b");
 s = b;
 output << s; 
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s", "CN","posU","lenU","Qoutlier","Qaberr","Source","Merge");
 s = b;
 output << s ;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s", "a5","a3","Nr5","Nr3","N5","N3");
 s = b;
 output << s << endl;
 //
 list<C_SV1>::iterator i;
 for(i=sv.evt.begin(); i != sv.evt.end(); ++i) {
 output << (*i) << endl;  
 }
 
 return output;
 }
 */

void C_SV::print(string & outfilename)
{
    // format: optimized to loadFragments.m matlab script  
    // open output binary file. bomb if unable to open
    fstream output(outfilename.c_str(), ios::out);
    output << *this;
}

//------------------------------------------------------------------------------
//  F & R clustered event class (retro)
//------------------------------------------------------------------------------ 
C_SVR1::C_SVR1() {
    pos=0;
    anchor=0;
    length=0;
    q=0;
    p5[0]=0;
    p5[1]=0;
    p3[0]=0;
    p3[1]=0;
    copynumber=0;
    // retro event info
    retro5.clear();
    retro3.clear();
    NfragCluster[0]=0;
    NfragCluster[1]=0;
    NfragCov=0;
    NfragCovOut[0]=0;
    NfragCovOut[1]=0;
    NfragCovExp=0;
    pmedian=0;
    Nconstrain[0]=0;
    Nconstrain[1]=0;
    a5=0;
    a3=0;
    posU=0;
    lenU=0;
    merge5=false;
    merge3=false;
    SampleMap.clear();
    type="";
    subtype="";
    Qsubtype=0;
    gap=0;
    Esense=0; 
    QsubtypeI=0;
    Qsubtype5=0;
    Qsubtype3=0;
    id=0;
}

/*
 ostream &operator<<(ostream &output, C_SVR1 & e1)
 {
 //----------------------------------------------------------------------------
 //  text file position coordinates are 1 based (not c convention)
 //----------------------------------------------------------------------------  
 char b [200];
 sprintf(b,"%10d\t%10d\t%4d\t%4d", e1.pos+1,e1.length,e1.anchor,e1.q);
 string s = b;
 output << s;
 
 C_cluster2d_element1 c1 = e1.cls5;
 sprintf(b,"\t%10d\t%10.0f\t%10.1f\t%10.0f\t%10.0f", int(c1.inp.size()),c1.mean[0],c1.std[0],c1.low[0],c1.high[0]);
 s = b;
 output << s;
 c1 = e1.cls3;
 sprintf(b,"\t%10d\t%10.0f\t%10.1f\t%10.0f\t%10.0f", int(c1.inp.size()),c1.mean[0],c1.std[0],c1.low[0],c1.high[0]);
 s = b;
 output << s;
 
 sprintf(b,"\t%10d\t%10d\t%10d\t%10d", e1.NfragCov,e1.NfragCovOut[0],e1.NfragCovOut[1],e1.NfragCovExp);
 s = b;
 output << s;
 
 sprintf(b,"\t%10d\t%10d\t%10d\t%10d", e1.p5[0]+1,e1.p5[1]+1,e1.p3[0]+1,e1.p3[1]+1);
 s = b;
 output << s;  
 
 C_SVcoverage1 q1=e1.cov;
 sprintf(b,"\t%10d\t%10d\t%10d\t%10d\t%10d\t%10.1f\t%10d", q1.N,q1.p0+1,q1.p1+1,q1.Nsite,q1.Nrepeat,q1.eN,e1.copynumber);
 s = b;
 output << s;
 
 sprintf(b,"\t%10d\t%10d\t%10d",e1.pmedian+1,e1.Nconstrain[0],e1.Nconstrain[1]);
 s = b;
 output << s;  
 
 sprintf(b,"\t%10d\t%10d",e1.posU,e1.lenU);
 s = b;
 output << s;
 
 //----------------------------------------------------------------------------
 // Read Groups
 //----------------------------------------------------------------------------
 std::map<unsigned int,unsigned int, less<unsigned int> >::iterator it;
 
 // iterate over read groups in 5' end;
 c1 = e1.cls5;
 sprintf(b,"\t%d", c1.N);
 s = b;
 for ( it=e1.ReadGroupMap5.begin() ; it != e1.ReadGroupMap5.end(); it++ ) {
 unsigned int RGC1 = (*it).first;
 unsigned int NF1  = e1.ReadGroupMap5[RGC1];
 sprintf(b,"|%u:%u", NF1,RGC1);
 s=s+b;
 }
 output << s;
 
 // iterate over read groups in 3' end;
 c1 = e1.cls3;
 sprintf(b,"\t%d", c1.N);
 s = b;
 for ( it=e1.ReadGroupMap3.begin() ; it != e1.ReadGroupMap3.end(); it++ ) {
 unsigned int RGC1 = (*it).first;
 unsigned int NF1  = e1.ReadGroupMap3[RGC1];
 sprintf(b,"|%u:%u", NF1,RGC1);
 s=s+b;
 }
 output << s;
 
 // merge
 sprintf(b,"\t%10d\t%10d", (e1.merge5)?1:0,(e1.merge3)?1:0);
 s = b;
 output << s;
 
 // end alignability
 sprintf(b,"\t%10.3f\t%10.3f",e1.a5,e1.a3);
 s = b;
 output << s;
 
 return output;
 }
 */

ostream &operator<<(ostream &output, C_SVR1 & e1)
{
    //----------------------------------------------------------------------------
    //  text file position coordinates are 1 based (not c convention)
    //----------------------------------------------------------------------------  
    char b [200];
    string s="INS:"+e1.type;
    string sf="PASS";
    if ((e1.NfragCluster[0]==0)|(e1.NfragCluster[1]==0)) {
		sf="1_SIDE";
	}
    sprintf(b,"%d\t%d\t%d\t.\t<%s>\t%d\t%s\t", e1.anchor, e1.pos+1,e1.id, s.c_str(),e1.q,sf.c_str());
    s = b;
    output << s;
    
    sprintf(b,"SVLEN=%d;NFF=%d;NFR=%d;CIPOS=%d,%d;", e1.length, e1.NfragCluster[0],e1.NfragCluster[1],-e1.posU,e1.lenU);
    s = b;
    output << s;
    
    int PA=(e1.NfragCluster[0]>0? e1.p5[0]-e1.pos:0);
    int PB=(e1.NfragCluster[0]>0? e1.p5[1]-e1.pos:0);
    int PC=(e1.NfragCluster[1]>0? e1.p3[0]-e1.pos:0);
    int PD=(e1.NfragCluster[1]>0? e1.p3[1]-e1.pos:0);
    
    sprintf(b,"PC1=%d,%d;PC2=%d,%d;",PA,PB,PC,PD);
    s = b;
    output << s;
    
    /*
	 sprintf(b,"NR=%d;ER=%.1f;MR=%d;", e1.cov.N,e1.cov.eN,(e1.merge5|e1.merge3)?1:0);
     s = b;
     output << s;  
     
     sprintf(b,"A5=%.2f;A3=%.2f;", e1.a5,e1.a3);
     s = b;
     output << s;  
     
     */
	
	
	s=e1.subtype;  
	size_t t=s.find(".");
	if (t!=string::npos) {
        s=s.substr(t+1);
	}
	
    //sprintf(b,"EL=%s;V0=%d;V1=%d;V5=%d;V3=%d;", s.c_str(),int(e1.Qsubtype),int(e1.QsubtypeI),int(e1.Qsubtype5),int(e1.Qsubtype3));
    sprintf(b,"EL=%s;", s.c_str());
    s = b;
    output << s;  
    
    sprintf(b,"GP=%d;ES=%.1f", e1.gap,e1.Esense);
    s = b;
    output << s;  
    
    return output;
}

C_SVR1& C_SVR1::operator=(const C_SVR1 &rhs)
{
    this->pos = rhs.pos;
    this->anchor = rhs.anchor;
    this->length = rhs.length;
    this->q = rhs.q;
    this->p5[0] = rhs.p5[0];
    this->p5[1] = rhs.p5[1];
    this->p3[0] = rhs.p3[0];
    this->p3[1] = rhs.p3[1];
    this->copynumber = rhs.copynumber;
    this->cov    = rhs.cov;  
    this->cls5   = rhs.cls5;  
    this->cls3   = rhs.cls3;  
    this->retro5 = rhs.retro5;  
    this->retro3 = rhs.retro3;  
    this->NfragCluster[0]= rhs.NfragCluster[0];
    this->NfragCluster[1]= rhs.NfragCluster[1];
    this->NfragCov       = rhs.NfragCov;
    this->NfragCovOut[0] = rhs.NfragCovOut[0];
    this->NfragCovOut[1] = rhs.NfragCovOut[1];
    this->NfragCovExp    = rhs.NfragCovExp;
    this->pmedian        = rhs.pmedian;
    this->Nconstrain[0]  = rhs.Nconstrain[0];
    this->Nconstrain[1]  = rhs.Nconstrain[1];
    this->a5             = rhs.a5;
    this->a3             = rhs.a3;
    this->posU           = rhs.posU;
    this->lenU           = rhs.lenU;
    this->ReadGroupMap5  = rhs.ReadGroupMap5;
    this->ReadGroupMap3  = rhs.ReadGroupMap3;
    this->merge5         = rhs.merge5;
    this->merge3         = rhs.merge3;
    this->gap            = rhs.gap;
    this->type           = rhs.type;
    this->subtype        = rhs.subtype;
    this->Qsubtype       = rhs.Qsubtype;
    this->QsubtypeI      = rhs.QsubtypeI;
    this->Qsubtype5      = rhs.Qsubtype5;
    this->Qsubtype3      = rhs.Qsubtype3;
    this->Esense         = rhs.Esense; 
    this->id             = rhs.id;
    this->SampleMap      = rhs.SampleMap;
    
    return *this;
}

int C_SVR1::operator==(const C_SVR1 &rhs) const
{
    if( this->pos != rhs.pos) return 0;
    if( this->anchor != rhs.anchor) return 0;
    if( this->length != rhs.length) return 0;
    if( this->p5[0]  != rhs.p5[0]) return 0;
    if( this->p5[1]  != rhs.p5[1]) return 0;
    if( this->p3[0]  != rhs.p3[0]) return 0;
    if( this->p3[1]  != rhs.p3[1]) return 0;
    if( this->q != rhs.q) return 0;
    if (this->copynumber != rhs.copynumber) return 0;
    if (this->NfragCov   != rhs.NfragCov) return 0;
    if (this->cls5.inp.size()!= rhs.cls5.inp.size()) return 0;
    if (this->cls3.inp.size()!= rhs.cls3.inp.size()) return 0;
    if (this->pmedian        != rhs.pmedian) return 0;
    if (this->gap            != rhs.gap) return 0;
    
    if( this->posU != rhs.posU) return 0;
    if( this->lenU != rhs.lenU) return 0;
    if (this->merge5 != rhs.merge5) return 0;
    if (this->merge3 != rhs.merge3) return 0;
    if (this->ReadGroupMap5.size() != rhs.ReadGroupMap5.size()) return 0;
    if (this->ReadGroupMap3.size() != rhs.ReadGroupMap3.size()) return 0;
    if( this->a5 != rhs.a5) return 0;
    if( this->a3 != rhs.a3) return 0;
    if( this->id != rhs.id) return 0;
    
    return 1;
}

int C_SVR1::operator<(const C_SVR1 &rhs) const
{
    if( this->anchor < rhs.anchor ) return 1;
    if( this->anchor > rhs.anchor ) return 0;
    if( this->pos < rhs.pos ) return 1;
    if( this->pos > rhs.pos ) return 0;
    if( this->cls5.inp.size() < rhs.cls5.inp.size() ) return 1;
    if( this->cls5.inp.size() > rhs.cls5.inp.size() ) return 0;
    if( this->cls3.inp.size() < rhs.cls3.inp.size()) return 1;
    if( this->cls3.inp.size() > rhs.cls3.inp.size()) return 0;
    if( this->id < rhs.id ) return 1;
    return 0;
}


//------------------------------------------------------------------------------
//  F & R clustered event class (inversions and retro)
//------------------------------------------------------------------------------ 
C_SVR::C_SVR() {
    typeName="z";
    contigName="";
    setName="";
    evt.clear(); 
}

C_SVR::C_SVR(C_contig  & c1,  RunControlParameters & pars) {
    typeName="z";
    contigName =c1.getContigName();
    setName =c1.setName;
    evt.clear();    
}

void C_SVR::finalize(C_contig  & contig, C_libraries & libraries, RunControlParameters & pars) {
    
    //----------------------------------------------------------------------------
    // sample by sample genotype info
    //----------------------------------------------------------------------------
    genotype(contig,libraries,pars);
    
    //----------------------------------------------------------------------------
    // loop over events / clusters within events to calculate:
    //    insertion length 
    //    most likely subtype
    //    element polarity
    //----------------------------------------------------------------------------
    
    // check for events, libraries
    if (evt.size()==0) { return;}
    if (libraries.libmap.size()==0) { return;}
    
    // number of characters in insertion element name 
    int LT=typeName.size();
    
    // declare maps for accumulating votes for element subclass 
    map<string, double, less<string> >  aMap, aMap5,aMap3;
    
    list<C_SVR1>::iterator ie,ie1,ie2;
    // first event
    ie1=evt.begin();
    // last event
    ie2=evt.end();  
    // loop over events
    int ne=0;
    
    // max a5 or a3 alignability metrics
    double aMax=0;
    
    for ( ie=ie1 ; ie != ie2; ++ie ) {
        ne++;
        aMap.clear();
        
        // id
        (*ie).id=ne;
        
        int n5=(*ie).retro5.size();
        int n3=(*ie).retro3.size();
        
        // find max alignability scale
        if ((*ie).a5>aMax) { aMax=(*ie).a5;}
        if ((*ie).a3>aMax) { aMax=(*ie).a3;}
        
        // loop over end
        for (int s=0; s<2; s++) {
            
            // pointer to umpair list for this end
            vector<C_umpair> *um;
            um= & (*ie).retro5;
            
            // pointer to accumulating map for this end
            std::map<string,double, std::less<string> > *am;
            am= & aMap5;
            
            // re-point pointers for 3' end
            if (s==1) {
                um= & (*ie).retro3;
                am= & aMap3;
            } 
            
            // reset map for this end
            (*am).clear();
            
            // number of supporing fragments at this end
            int Np=um->size();
            
            if (Np<1) {continue;};
            
            // loop over frags
            for (int i=0; i<Np; i++)  {
                
                // umpair for this fragment
                C_umpair p1=(*um)[i];
                
                // number of element mapppings
                short nM=p1.nmap;
                
                // element index for this fragment
                short k=p1.read[1].anchor;
                // element name 
                string nameM=libraries.anchorinfo.names[k];
                
                // accumulate votes (weighted by 1/nM) for this element
                aMap[nameM]+=1.0/nM;
                (*am)[nameM]+=1.0/nM;
            }
            
        }    
        
        // remove generic elements from voting (must be an element subclass)
        if (aMap.size()>1) {
            string genericName="L1.L1";
            if (aMap[genericName]>0) { aMap[genericName]=0; }
            if ((aMap5.size()>1)&(aMap5[genericName]>0)) { aMap5[genericName]=0;} 
            if ((aMap3.size()>1)&(aMap3[genericName]>0)) { aMap3[genericName]=0;} 
            genericName="ALU.ALUY";
            if (aMap[genericName]>0) { aMap[genericName]=0; }
            if ((aMap5.size()>1)&(aMap5[genericName]>0)) { aMap5[genericName]=0;} 
            if ((aMap3.size()>1)&(aMap3[genericName]>0)) { aMap3[genericName]=0;} 
        }
        
        // most popular (weighted 1/nM voting) element name
        vector<string> ename;
        ename = sortKeysByValue(aMap,true);
        
        
        string stype=ename[0];    
        
        (*ie).subtype=stype;
        
        // most popular votes
        double v1=aMap[ename[0]];
        string ename1=ename[0];
        
        // next popular votes
        double  v2=0;
        if (ename.size()>1) { v2=aMap[ename[1]];};
        
        
        
        double  vB=0;
        char SC1=ename[0][2*LT+1];
        
        // different subclass level 1 next popular votes
        double vtot=0;
        for (int i=0; i<int(ename.size()); i++) {
            double  v=aMap[ename[i]];
            vtot+=v;
            char SCB=ename[i][2*LT+1];
            if (SCB!=SC1) {
                if (v>vB) {vB=v;}
            }
        }
        
        // confidence metrics for subclasses
        (*ie).Qsubtype=99.*v1/vtot; //10.0*v1*(v1-v2);
        (*ie).QsubtypeI=99.*(v1-vB)/vtot; //10.0*v1*(v1-vB);
        
        // 5' most popular (weighted 1/nM voting) element name
        if (aMap5.size()>0) {
            ename= sortKeysByValue(aMap5,true);
            (*ie).subtype5=ename[0];
            // double v5=aMap5[ename[0]];
            vtot=1.e-10;
            for (int i=0; i<int(ename.size()); i++) {
                vtot+=aMap5[ename[i]];
            }      
            (*ie).Qsubtype5=99.*aMap5[ename1]/vtot; //10.0*v1*(v1-v2);
        }
        
        // 3' most popular (weighted 1/nM voting) element name
        if (aMap3.size()>0) {
            ename= sortKeysByValue(aMap3,true);
            (*ie).subtype3=ename[0];
            //double v3 =aMap3[ename[0]];
            vtot=1.e-10;
            for (int i=0; i<int(ename.size()); i++) {
                vtot+=aMap3[ename[i]];
            }      
            (*ie).Qsubtype3=99.*aMap3[ename1]/vtot; //10.0*v1*(v1-v2);
        }
        
        int pMin=1000000;
        int pMax=0;
        double  sen1[2];
        
        // loop over ends
        for (int s=0; s<2; s++) {
            
            // pointer to umpair list for this end
            vector<C_umpair> *um;
            um= & (*ie).retro5;
            
            // re-point pointers for 3' end
            if (s==1) {
                um= & (*ie).retro3;
            } 
            
            // reset sense sensor
            sen1[s]=0.0;
            int n1=0;
            
            // number of supporing fragments at this end
            int Np=um->size();
            
            // loop over frags
            for (int i=0; i<Np; i++)  {
                C_umpair p1=(*um)[i];
                short aM=p1.read[1].anchor;
                string nameM=libraries.anchorinfo.names[aM];
                
                // only when this is the top element
                if (nameM==(*ie).subtype) {
                    int pM1=p1.read[1].pos;
                    int pM2=p1.read[1].pos+p1.read[1].len;
                    if (pM1<pMin) {pMin=pM1;};
                    if (pM2>pMax) {pMax=pM2;};
                    sen1[s]+=(p1.read[1].sense=='F'?1.0:-1.0);
                    n1++;
                }
                
            }
            
            // normalize sense
            if (n1>0) { sen1[s]=sen1[s]/n1; }
            
        }
        if ((n5>0)&(n3>0)) { 
            
            // element sense // +1 for F, -1 for R
            (*ie).Esense=double(n5*sen1[1]-n3*sen1[0])/double(n5+n3);  
            
            // element insertion length 
            (*ie).length=pMax-pMin;     
        } else if (n5>0) { 
            
            // 5' only
            (*ie).Esense=sen1[0];
        } else {
            
            // 3' only
            (*ie).Esense=sen1[1];
        }
        
        
    }
    
    // normalize alignability to <=1.0
    for ( ie=evt.begin(); ie != evt.end(); ++ie ) {
        (*ie).a5=(*ie).a5/aMax;
        (*ie).a3=(*ie).a3/aMax;
    }
    
} 




void C_SVR::genotype(C_contig  & contig, C_libraries & libraries, RunControlParameters & pars)
{
    // check for events, libraries
    if (evt.size()==0) { return;}
    if (libraries.libmap.size()==0) { return;}
    
    map<string, double, less<string> >  rX = libraries.readFractionSamples();
    
    list<C_SVR1>::iterator ie,ie1,ie2;
    //C_SVR1 e1;
    unsigned int ReadGroupCode1;
    string SAM;
    int NF1;
    
    ie1=evt.begin();
    
    // know when to stop
    ie2=evt.end();  
    
    list<C_localpair>::iterator i;
    // get Minimum Q mapping value for unique read
    int Qmin  = pars.getQmin();
    // loop over pairs
    int np=0;
    for(i=contig.localpairs.begin(); i != contig.localpairs.end(); ++i) {
        np++;
        // skip low mapping quality fragments
        if (((*i).q1<Qmin)&((*i).q2<Qmin)) {
            continue;
        }
        // fragment length limits for this library
        ReadGroupCode1 = (*i).ReadGroupCode;  
        //int LM1=libraries.libmap[ReadGroupCode1].LM;
        int LMhigh1=libraries.libmap[ReadGroupCode1].LMhigh;
        int LMlow1=libraries.libmap[ReadGroupCode1].LMlow;
        
        // fragment length
        int lm = (*i).lm;
        // demand usual fragment
        if ((lm<LMlow1)|(lm>LMhigh1)) {continue; }
        
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        
        // add only the non-read part of the fragment
        int p1 = (*i).pos+(*i).len1;
        int p2 = p1+(*i).lm-(*i).len2;
        
        // loop over events
        int ne=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            
            ne++;
            //e1 = *ie;
            // check for reads starts inside window (pos+/- 100bp)
            if (( abs(int((*i).pos)-int((*ie).pos)) <=100)&((*i).q1>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            // other end
            if (( abs(int((*i).pos+(*i).lm-(*i).len2)-int((*ie).pos)) <=100)&((*i).q2>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            
            int P1=(*ie).pos-(*ie).posU/2;       
            if ((P1-1000)>p2) { break;}         
            if (P1<p1) { continue;}       
            int P2=(*ie).pos+(*ie).posU/2;            
            if (P2>p2) { continue;}
            
            // select spanning fragments - both above Qmin 
            if (((*i).q1<Qmin)|((*i).q2<Qmin)) {
                continue;
            } 
            
            (*ie).SampleMap[SAM].NN++;
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
        
    }
    
    ie1=evt.begin();
    
    std::map<unsigned int,unsigned int, less<unsigned int> >::iterator it;  
    // iterate over read groups in 5' end;
    for ( ie=ie1 ; ie != ie2; ++ie ) {
        //e1 = *ie;
        for ( it=(*ie).ReadGroupMap5.begin() ; it != (*ie).ReadGroupMap5.end(); it++ ) {
            ReadGroupCode1 = (*it).first;
            NF1  = (*ie).ReadGroupMap5[ReadGroupCode1];
            SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
            (*ie).SampleMap[SAM].N5+=NF1;
            (*ie).SampleMap[SAM].N+=NF1;
        }     
        for ( it=(*ie).ReadGroupMap3.begin() ; it != (*ie).ReadGroupMap3.end(); it++ ) {
            ReadGroupCode1 = (*it).first;
            NF1  = (*ie).ReadGroupMap3[ReadGroupCode1];
            SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
            (*ie).SampleMap[SAM].N3+=NF1;
            (*ie).SampleMap[SAM].N+=NF1;
        }      
        
        // loop over libraries for estimated number of reads 
        //map<string, double, less<string> > rX = libraries.readFractionSamples();
        map<string, double, less<string> >::iterator irX;
        for ( irX=rX.begin() ; irX != rX.end(); irX++ ) {
            SAM=(*irX).first;
            (*ie).SampleMap[SAM].ER=double((*ie).cov.eN)*rX[SAM];
        }
    }
    
    // loop over single ends
    
    // fill depth from cross pairs (not so many...)
    int N = contig.crosspairs.size();
    list<C_crosspair>::iterator ix;
    ie1=evt.begin();
    for(ix=contig.crosspairs.begin(); ix != contig.crosspairs.end(); ++ix) {
        // skip low mapping quality fragments
        if ((*ix).read[0].q<Qmin) {continue;}
        // fragment length limits for this library
        ReadGroupCode1 = (*ix).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        // loop over events
        int nx=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nx++;
            // check for reads starts inside window (pos+/- 100bp)
            if ( abs(int((*ix).read[0].pos)-int((*ie).pos)) <=100) {       
                (*ie).SampleMap[SAM].NR++;
            } else if (int((*ix).read[0].pos)<int((*ie).pos-100)) {
                break;
            }
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    // fill depth from dangle starts 
    N = contig.dangle.size();
    list<C_singleEnd>::iterator is; 
    for(is=contig.dangle.begin(); is != contig.dangle.end(); ++is) {
        // skip low mapping quality fragments
        if ((*is).q<Qmin) {continue;}
        ReadGroupCode1 = (*is).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        // loop over events
        int nd=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nd++;
            // check for reads starts inside window (pos+/- 100bp)
            if ( abs(int((*is).pos)-int((*ie).pos)) <=100) {       
                (*ie).SampleMap[SAM].NR++;
            } else if (int((*is).pos)<int((*ie).pos-100)) {
                break;
            }
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    // fill depth from umpairs 
    N = contig.umpairs.size();
    list<C_umpair>::iterator iu;
    
    for(iu=contig.umpairs.begin(); iu != contig.umpairs.end(); ++iu) {
        if ((*iu).read[0].q<Qmin) {continue;}
        ReadGroupCode1 = (*iu).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        // loop over events
        int nd=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nd++;
            // check for reads starts inside window (pos+/- 100bp)
            if ( abs(int((*iu).read[0].pos)-int((*ie).pos)) <=100) {       
                (*ie).SampleMap[SAM].NR++;
            } else if (int((*iu).read[0].pos)<int((*ie).pos-100)) {
                break;
            }
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    // fill depth from unique ends 
    N = contig.singleton.size();
    for(is=contig.singleton.begin(); is != contig.singleton.end(); ++is) {
        if ((*is).q<Qmin) {continue;}
        ReadGroupCode1 = (*is).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        // loop over events
        int nd=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nd++;
            // check for reads starts inside window (pos+/- 100bp)
            if ( abs(int((*is).pos)-int((*ie).pos)) <=100) {       
                (*ie).SampleMap[SAM].NR++;
            } else if (int((*is).pos)<int((*ie).pos-100)) {
                break;
            }
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    
    
}

/*
 ostream &operator<<(ostream &output, C_SVR & sv)
 {
 
 output << sv.typeName << "\t " << sv.setName << "\t " << sv.contigName << endl;
 int Nevt = sv.evt.size();
 output << Nevt << endl;
 char b[200];
 sprintf(b,"%10s\t%10s\t%4s\t%4s", "pos","length","anc","q");
 string s = b;
 output << s;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s\t%10s","N5","p5","std5","low5","high5");
 s = b;
 output << s;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s\t%10s","N3","p3","std3","low3","high3");
 s = b;
 output << s;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s", "Nfrag","NfragPre","NFragPost","NfragExp");
 s = b;
 output << s;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s", "p5a","p5b","p3a","p3b");
 s = b;
 output << s; 
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s", "Nread","p0","p1","Nsite","Nrepeat","eN","CN");
 s = b;
 output << s; 
 sprintf(b,"\t%10s\t%10s\t%10s","pmedian","N5con","N3con");
 s = b;
 output << s;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s", "posU","lenU","Source5","Source3","Merge5","Merge3");
 s = b;
 output << s ;
 sprintf(b,"\t%10s\t%10s", "a5","a3");
 s = b;
 output << s ;
 // samples
 size_t NS=sv.samples.size();
 sprintf(b,"\t%10d",int(NS));
 s = b;
 output << s ;
 for (int ns=0; ns<NS; ns++) {
 s=sv.samples[ns];
 sprintf(b,"\t%s",s.c_str());
 s = b;
 output << s ;
 }  
 output << endl;
 
 
 //
 list<C_SVR1>::iterator i;
 for(i=sv.evt.begin(); i != sv.evt.end(); ++i) {
 output << (*i);
 sprintf(b,"\t%10d",int(NS));
 s = b;
 output << s ;
 for (int ns=0; ns<NS; ns++) {
 s=sv.samples[ns];
 sprintf(b,"\t%d:%d:%d:%d",(*i).SampleMap[s].N,(*i).SampleMap[s].NN,(*i).SampleMap[s].N5,(*i).SampleMap[s].N3);
 s = b;
 output << s;
 }
 output  << endl;  
 }  
 return output;
 }
 */

ostream &operator<<(ostream &output, C_SVR & sv)
{
    
    // SVCF header
    output  << sv.SVCF;
    
    size_t NS = sv.samples.size();
    
    char b[200];
    string s;
    
    // format 
    size_t NFMT=sv.SVCF.FMT.size();
    string FMT=sv.SVCF.FMT[0].Id;
    for (int f=1; f<int(NFMT); f++) {
        FMT=FMT+":"+sv.SVCF.FMT[f].Id;
    }
    // events
    list<C_SVR1>::iterator i;
    for(i=sv.evt.begin(); i != sv.evt.end(); ++i) {
        output << (*i);
        
        output << "\t" << FMT ;
        
        for (int ns=0; ns<int(NS); ns++) {
            s=sv.samples[ns];
            // double cn = 2.0*double((*i).SampleMap[s].NR)/((*i).SampleMap[s].ER+0.01);
            //sprintf(b,"\t%d:%d:%d:%d:%d:%.1f",(*i).SampleMap[s].N,(*i).SampleMap[s].NN,(*i).SampleMap[s].N5,(*i).SampleMap[s].N3,(*i).SampleMap[s].NR,(*i).SampleMap[s].ER); //cn);
            sprintf(b,"\t%d:%d",(*i).SampleMap[s].N5,(*i).SampleMap[s].N3); //cn);
            s = b;
            output << s;
        }
        output  << endl;  
    }  
    return output;
}


void C_SVR::print(string & outfilename)
{
    // format: optimized to loadFragments.m matlab script  
    // open output binary file. bomb if unable to open
    fstream output(outfilename.c_str(), ios::out);
    output << *this;
}

//------------------------------------------------------------------------------
// Inversion F & R event class
//------------------------------------------------------------------------------
C_SVV1::C_SVV1() {
    type="inversion";
    pos=0;
    anchor=0;
    length=0;
    q=0;
    p5[0]=0;
    p5[1]=0;
    p3[0]=0;
    p3[1]=0;
    copynumber=0;
    // retro event info
    pair5.clear();
    pair3.clear();
    NfragCluster[0]=0;
    NfragCluster[1]=0;
    NfragCov=0;
    NfragCovOut[0]=0;
    NfragCovOut[1]=0;
    NfragCovExp=0;
    
    posU=0;
    lenU=0;
    qOutlier=0;
    qAberrantLM=0;
    merge5=false;
    merge3=false;
    a5[0]=0;
    a3[0]=0;
    a5[1]=0;
    a3[1]=0;
    
    SampleMap.clear();
    id=0;
    
}

ostream &operator<<(ostream &output, C_SVV1 & e1)
{
    //----------------------------------------------------------------------------
    //  text file position coordinates are 1 based (not c convention)
    //----------------------------------------------------------------------------  
    char b [200];
    string s=e1.type;
    sprintf(b,"%d\t%d\t%d\t.\t%s\t%d\t0\t", e1.anchor, e1.pos+1,e1.id, s.c_str(),e1.q);
    s = b;
    output << s;
    
    sprintf(b,"L=%d;NF=%d;UP=%d;UL=%d;", e1.length, int(e1.cls5.inp.size())+int(e1.cls3.inp.size()), e1.posU,e1.lenU);
    s = b;
    output << s;
    
    int PA=(e1.NfragCluster[0]>0? e1.p5[0]-e1.pos:0);
    int PB=(e1.NfragCluster[0]>0? e1.p5[1]-e1.pos:0);
    int PC=(e1.NfragCluster[1]>0? e1.p3[0]-e1.pos:0);
    int PD=(e1.NfragCluster[1]>0? e1.p3[1]-e1.pos:0);
    
    sprintf(b,"N5=%d;N3=%d;PA=%d;PB=%d;PC=%d;PD=%d;",e1.NfragCluster[0],e1.NfragCluster[1],PA,PB,PC,PD);
    s = b;
    output << s;
    
    sprintf(b,"NR=%d;ER=%.1f;MR=%d;", e1.cov.N,e1.cov.eN,(e1.merge5|e1.merge3)?1:0);
    s = b;
    output << s;  
    
    sprintf(b,"A5=%.2f;A3=%.2f;", e1.a5[0],e1.a3[0]);
    s = b;
    output << s;  
    
    return output;
}


/*
 ostream &operator<<(ostream &output, C_SVV1 & e1)
 {
 //----------------------------------------------------------------------------
 //  text file position coordinates are 1 based (not c convention)
 //----------------------------------------------------------------------------  
 char b [200];
 sprintf(b,"%10d\t%10d\t%4d\t%4d", e1.pos+1,e1.length,e1.anchor,e1.q);
 string s = b;
 output << s;
 
 C_cluster2d_element1 c1 = e1.cls5;
 sprintf(b,"\t%10d\t%10.0f\t%10.1f\t%10.0f\t%10.0f", int(c1.inp.size()),c1.mean[0],c1.std[0],c1.low[0],c1.high[0]);
 s = b;
 output << s;
 sprintf(b,"\t%10.0f\t%10.1f\t%10.0f\t%10.0f", c1.mean[1],c1.std[1],c1.low[1],c1.high[1]);
 s = b;
 output << s;
 c1 = e1.cls3;
 sprintf(b,"\t%10d\t%10.0f\t%10.1f\t%10.0f\t%10.0f", int(c1.inp.size()),c1.mean[0],c1.std[0],c1.low[0],c1.high[0]);
 s = b;
 output << s;
 sprintf(b,"\t%10.0f\t%10.1f\t%10.0f\t%10.0f", c1.mean[1],c1.std[1],c1.low[1],c1.high[1]);
 s = b;
 output << s;
 
 sprintf(b,"\t%10d\t%10d\t%10d\t%10d", e1.NfragCov,e1.NfragCovOut[0],e1.NfragCovOut[1],e1.NfragCovExp);
 s = b;
 output << s;
 
 sprintf(b,"\t%10d\t%10d\t%10d\t%10d", e1.p5[0]+1,e1.p5[1]+1,e1.p3[0]+1,e1.p3[1]+1);
 s = b;
 output << s;  
 
 C_SVcoverage1 q1=e1.cov;
 sprintf(b,"\t%10d\t%10d\t%10d\t%10d\t%10d\t%10.1f\t%10d", q1.N,q1.p0+1,q1.p1+1,q1.Nsite,q1.Nrepeat,q1.eN,e1.copynumber);
 s = b;
 output << s;
 
 sprintf(b,"\t%10d\t%10d\t%10d\t%10d", e1.posU,e1.lenU,e1.qOutlier,e1.qAberrantLM);
 s = b;
 output << s;
 
 // iterate over read groups in 5' end;
 c1 = e1.cls5;
 sprintf(b,"\t%d", c1.N);
 s = b;
 std::map<unsigned int,unsigned int, less<unsigned int> >::iterator it;
 for ( it=e1.ReadGroupMap5.begin() ; it != e1.ReadGroupMap5.end(); it++ ) {
 unsigned int RGC1 = (*it).first;
 unsigned int NF1  = e1.ReadGroupMap5[RGC1];
 sprintf(b,"|%u:%u", NF1,RGC1);
 s=s+b;
 }
 output << s;
 
 // iterate over read groups in 3' end;
 c1 = e1.cls3;
 sprintf(b,"\t%d", c1.N);
 s = b;
 for ( it=e1.ReadGroupMap3.begin() ; it != e1.ReadGroupMap3.end(); it++ ) {
 unsigned int RGC1 = (*it).first;
 unsigned int NF1  = e1.ReadGroupMap3[RGC1];
 sprintf(b,"|%u:%u", NF1,RGC1);
 s=s+b;
 }
 output << s;
 
 // merge
 sprintf(b,"\t%10d\t%10d", (e1.merge5)?1:0,(e1.merge3)?1:0);
 s = b;
 output << s;
 
 // end alignability
 float a5=(e1.a5[0]>e1.a5[1]? e1.a5[0]: e1.a5[1] );
 float a3=(e1.a3[0]>e1.a3[1]? e1.a3[0]: e1.a3[1] );
 sprintf(b,"\t%10.3f\t%10.3f\t%10d\t%10d\t%10d\t%10d",a5,a3,e1.cov5.Nrepeat,e1.cov3.Nrepeat,e1.cov5.N,e1.cov3.N);
 s = b;
 output << s;
 
 
 
 
 
 return output;
 }
 
 */

C_SVV1& C_SVV1::operator=(const C_SVV1 &rhs)
{
    this->pos = rhs.pos;
    this->anchor = rhs.anchor;
    this->length = rhs.length;
    this->q = rhs.q;
    this->p5[0] = rhs.p5[0];
    this->p5[1] = rhs.p5[1];
    this->p3[0] = rhs.p3[0];
    this->p3[1] = rhs.p3[1];
    this->copynumber = rhs.copynumber;
    this->cov    = rhs.cov;  
    this->cls5   = rhs.cls5;  
    this->cls3   = rhs.cls3;  
    this->pair5 = rhs.pair5;  
    this->pair3 = rhs.pair3;  
    this->NfragCluster[0]= rhs.NfragCluster[0];
    this->NfragCluster[1]= rhs.NfragCluster[1];
    this->NfragCov       = rhs.NfragCov;
    this->NfragCovOut[0] = rhs.NfragCovOut[0];
    this->NfragCovOut[1] = rhs.NfragCovOut[1];
    this->NfragCovExp    = rhs.NfragCovExp;
    //
    this->cov5  = rhs.cov5;
    this->cov3  = rhs.cov3;
    this->a5[0] = rhs.a5[0];
    this->a5[1] = rhs.a5[1];
    this->a3[0] = rhs.a3[0];
    this->a3[1] = rhs.a3[1];
    this->posU  = rhs.posU;
    this->lenU  = rhs.lenU;
    this->qOutlier     = rhs.qOutlier;
    this->qAberrantLM  = rhs.qAberrantLM;
    this->ReadGroupMap5= rhs.ReadGroupMap5;
    this->ReadGroupMap3= rhs.ReadGroupMap3;
    this->merge5 = rhs.merge5;
    this->merge3 = rhs.merge3;
    this->id             = rhs.id;
    this->type           = rhs.type;
    this->SampleMap      = rhs.SampleMap;
    return *this;
}

int C_SVV1::operator==(const C_SVV1 &rhs) const
{
    if( this->pos != rhs.pos) return 0;
    if( this->anchor != rhs.anchor) return 0;
    if( this->length != rhs.length) return 0;
    if( this->p5[0] != rhs.p5[0]) return 0;
    if( this->p5[1] != rhs.p5[1]) return 0;
    if( this->p3[0] != rhs.p3[0]) return 0;
    if( this->p3[1] != rhs.p3[1]) return 0;
    if( this->q != rhs.q) return 0;
    if (this->copynumber != rhs.copynumber) return 0;
    if (this->NfragCov   != rhs.NfragCov) return 0;
    if (this->cls5.inp.size()!= rhs.cls5.inp.size()) return 0;
    if (this->cls3.inp.size()!= rhs.cls3.inp.size()) return 0;
    
    if( this->posU != rhs.posU) return 0;
    if( this->lenU != rhs.lenU) return 0;
    if( this->qOutlier != rhs.qOutlier) return 0;
    if( this->qAberrantLM != rhs.qAberrantLM) return 0;
    if (this->merge5 != rhs.merge5) return 0;
    if (this->merge3 != rhs.merge3) return 0;
    if (this->ReadGroupMap5.size() != rhs.ReadGroupMap5.size()) return 0;
    if (this->ReadGroupMap3.size() != rhs.ReadGroupMap3.size()) return 0;
    if( this->a5[0] != rhs.a5[0]) return 0;
    if( this->a3[0] != rhs.a3[0]) return 0;
    if( this->a5[1] != rhs.a5[1]) return 0;
    if( this->a3[1] != rhs.a3[1]) return 0;
    
    return 1;
}

int C_SVV1::operator<(const C_SVV1 &rhs) const
{
    if( this->anchor < rhs.anchor ) return 1;
    if( this->anchor == rhs.anchor && this->pos < rhs.pos ) return 1;
    if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->cls5.inp.size() < rhs.cls5.inp.size() ) return 1;
    if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->cls5.inp.size() == rhs.cls5.inp.size()  
       && this->cls3.inp.size() < rhs.cls3.inp.size()) return 1;
    return 0;
}


//------------------------------------------------------------------------------
//  F & R clustered event class (inversions)
//------------------------------------------------------------------------------ 
C_SVV::C_SVV() {
    typeName="v";
    contigName="";
    setName="";
    evt.clear(); 
}

C_SVV::C_SVV(C_contig  & c1,  RunControlParameters & par) {
    typeName="v";
    contigName =c1.getContigName();
    setName =c1.setName;
    evt.clear(); 
}

void C_SVV::finalize(C_contig  & contig, C_libraries & libraries, RunControlParameters & pars) {
    
    //----------------------------------------------------------------------------
    // sample by sample genotype info
    //----------------------------------------------------------------------------
    genotype(contig,libraries,pars);
    
    //----------------------------------------------------------------------------
    // loop over events / clusters within events to calculate:
    //----------------------------------------------------------------------------
    
    // check for events, libraries
    if (evt.size()==0) { return;}
    if (libraries.libmap.size()==0) { return;}
    
    list<C_SVV1>::iterator ie,ie1,ie2;
    // first event
    ie1=evt.begin();
    // last event
    ie2=evt.end();  
    // loop over events
    int ne=0;
    
    // max a5 or a3 alignability metrics
    double aMax=0;
    
    for ( ie=ie1 ; ie != ie2; ++ie ) {
        ne++;
        
        // id
        (*ie).id=ne;
        
        // find max alignability scale
        if ((*ie).a5[0]>aMax) { aMax=(*ie).a5[0];}
        if ((*ie).a3[0]>aMax) { aMax=(*ie).a3[0];}
        if ((*ie).a5[1]>aMax) { aMax=(*ie).a5[1];}
        if ((*ie).a3[1]>aMax) { aMax=(*ie).a3[1];}
        
    }
    
    // normalize alignability to <=1.0
    for ( ie=evt.begin(); ie != evt.end(); ++ie ) {
        (*ie).a5[0]=(*ie).a5[0]/aMax;
        (*ie).a3[0]=(*ie).a3[0]/aMax;
        (*ie).a5[1]=(*ie).a5[1]/aMax;
        (*ie).a3[1]=(*ie).a3[1]/aMax;
    }
    
} 




void C_SVV::genotype(C_contig  & contig, C_libraries & libraries, RunControlParameters & pars)
{
    // check for events, libraries
    if (evt.size()==0) { return;}
    if (libraries.libmap.size()==0) { return;}
    
    map<string, double, less<string> >  rX = libraries.readFractionSamples();
    
    list<C_SVV1>::iterator ie,ie1,ie2;
    //C_SVR1 e1;
    unsigned int ReadGroupCode1;
    string SAM;
    int NF1;
    
    ie1=evt.begin();
    
    // know when to stop
    ie2=evt.end();  
    
    list<C_localpair>::iterator i;
    // get Minimum Q mapping value for unique read
    int Qmin  = pars.getQmin();
    // loop over pairs
    int np=0;
    for(i=contig.localpairs.begin(); i != contig.localpairs.end(); ++i) {
        np++;
        // skip low mapping quality fragments
        if (((*i).q1<Qmin)&((*i).q2<Qmin)) {
            continue;
        }
        // fragment length limits for this library
        ReadGroupCode1 = (*i).ReadGroupCode;  
        //int LM1=libraries.libmap[ReadGroupCode1].LM;
        int LMhigh1=libraries.libmap[ReadGroupCode1].LMhigh;
        int LMlow1=libraries.libmap[ReadGroupCode1].LMlow;
        
        // fragment length
        int lm = (*i).lm;
        // demand usual fragment
        if ((lm<LMlow1)|(lm>LMhigh1)) {continue; }
        
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        
        // add only the non-read part of the fragment
        int p0 = (*i).pos;
        int p1 = (*i).pos+(*i).len1;
        int p2 = p1+(*i).lm-(*i).len2;
        
        // loop over events
        int ne=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            
            ne++;
            
            // two spanning windows eP1a-eP1b and eP2b-eP2b  
            int eP1a=(*ie).pos-(*ie).posU/2;     
            int eP1b=int((*ie).pos)+int((*ie).posU/2);
            int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
            int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
            
            // check for reads starts inside event window eP1b-eP2a
            if ((p0>=eP1b)&(p0<=eP2a)&((*i).q1>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            // other end of fragment
            if ((p2>=eP1a)&(p2<=eP2a)&((*i).q1>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            
            // beginning of event well beyond this fragment
            if ((eP1a-1000)>p2) { break;}        
            
            // end of event not yet reaching this fragment
            if (eP2b<p1) { continue;}       
            
            // select spanning fragments - both above Qmin 
            if (((*i).q1<Qmin)|((*i).q2<Qmin)) {
                continue;
            } 
            
            // spanning 5' break
            if ((p1<eP1a)&(p2>eP1b)) {
                (*ie).SampleMap[SAM].NN++;
            }
            // spanning 3' break
            if ((p1<eP2a)&(p2>eP2b)) {
                (*ie).SampleMap[SAM].NN++;
            }
            
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
        
    }
    
    ie1=evt.begin();
    
    std::map<unsigned int,unsigned int, less<unsigned int> >::iterator it;  
    // iterate over read groups in 5' end;
    for ( ie=ie1 ; ie != ie2; ++ie ) {
        //e1 = *ie;
        for ( it=(*ie).ReadGroupMap5.begin() ; it != (*ie).ReadGroupMap5.end(); it++ ) {
            ReadGroupCode1 = (*it).first;
            NF1  = (*ie).ReadGroupMap5[ReadGroupCode1];
            SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
            (*ie).SampleMap[SAM].N5+=NF1;
            (*ie).SampleMap[SAM].N+=NF1;
        }     
        for ( it=(*ie).ReadGroupMap3.begin() ; it != (*ie).ReadGroupMap3.end(); it++ ) {
            ReadGroupCode1 = (*it).first;
            NF1  = (*ie).ReadGroupMap3[ReadGroupCode1];
            SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
            (*ie).SampleMap[SAM].N3+=NF1;
            (*ie).SampleMap[SAM].N+=NF1;
        }      
        
        // loop over libraries for estimated number of reads 
        //map<string, double, less<string> > rX = libraries.readFractionSamples();
        map<string, double, less<string> >::iterator irX;
        for ( irX=rX.begin() ; irX != rX.end(); irX++ ) {
            SAM=(*irX).first;
            (*ie).SampleMap[SAM].ER=double((*ie).cov.eN)*rX[SAM];
        }
    }
    
    // loop over single ends
    
    // fill depth from cross pairs (not so many...)
    int N = contig.crosspairs.size();
    list<C_crosspair>::iterator ix;
    ie1=evt.begin();
    for(ix=contig.crosspairs.begin(); ix != contig.crosspairs.end(); ++ix) {
        // skip low mapping quality fragments
        if ((*ix).read[0].q<Qmin) {continue;}
        // fragment length limits for this library
        ReadGroupCode1 = (*ix).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        
        // read start
        int p1=int((*ix).read[0].pos);
        
        // loop over events
        int nx=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nx++;
            
            int eP1a=(*ie).pos-(*ie).posU/2;     
            int eP1b=int((*ie).pos)+int((*ie).posU/2);
            int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
            int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
            
            // beginning of event well beyond this fragment
            if ((eP1a-1000)>p1) { break;}        
            
            // end of event not yet reaching this fragment
            if (eP2b<p1) { continue;}       
            
            // check for reads starts inside event window eP1b-eP2a
            if ((p1>=eP1b)&(p1<=eP2a)&((*ix).read[0].q>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    // fill depth from dangle starts 
    N = contig.dangle.size();
    list<C_singleEnd>::iterator is; 
    for(is=contig.dangle.begin(); is != contig.dangle.end(); ++is) {
        // skip low mapping quality fragments
        if ((*is).q<Qmin) {continue;}
        ReadGroupCode1 = (*is).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        
        // read start
        int p1=int((*is).pos);
        
        // loop over events
        int nd=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nd++;
            
            int eP1a=(*ie).pos-(*ie).posU/2;     
            int eP1b=int((*ie).pos)+int((*ie).posU/2);
            int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
            int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
            
            // beginning of event well beyond this fragment
            if ((eP1a-1000)>p1) { break;}        
            
            // end of event not yet reaching this fragment
            if (eP2b<p1) { continue;}       
            
            // check for reads starts inside event window eP1b-eP2a
            if ((p1>=eP1b)&(p1<=eP2a)&((*is).q>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    // fill depth from umpairs 
    N = contig.umpairs.size();
    list<C_umpair>::iterator iu;
    
    for(iu=contig.umpairs.begin(); iu != contig.umpairs.end(); ++iu) {
        if ((*iu).read[0].q<Qmin) {continue;}
        ReadGroupCode1 = (*iu).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        
        // read start
        int p1=int((*iu).read[0].pos);
        
        // loop over events
        int nd=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nd++;
            
            int eP1a=(*ie).pos-(*ie).posU/2;     
            int eP1b=int((*ie).pos)+int((*ie).posU/2);
            int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
            int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
            
            // beginning of event well beyond this fragment
            if ((eP1a-1000)>p1) { break;}        
            
            // end of event not yet reaching this fragment
            if (eP2b<p1) { continue;}       
            
            // check for reads starts inside event window eP1b-eP2a
            if ((p1>=eP1b)&(p1<=eP2a)&((*iu).read[0].q>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    // fill depth from unique ends 
    N = contig.singleton.size();
    for(is=contig.singleton.begin(); is != contig.singleton.end(); ++is) {
        if ((*is).q<Qmin) {continue;}
        ReadGroupCode1 = (*is).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        
        // read start
        int p1=int((*is).pos);
        
        // loop over events
        int nd=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nd++;
            int eP1a=(*ie).pos-(*ie).posU/2;     
            int eP1b=int((*ie).pos)+int((*ie).posU/2);
            int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
            int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
            
            // beginning of event well beyond this fragment
            if ((eP1a-1000)>p1) { break;}        
            
            // end of event not yet reaching this fragment
            if (eP2b<p1) { continue;}       
            
            // check for reads starts inside event window eP1b-eP2a
            if ((p1>=eP1b)&(p1<=eP2a)&((*is).q>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    
}

ostream &operator<<(ostream &output, C_SVV & sv)
{
    
    // SVCF header
    output  << sv.SVCF;
    
    size_t NS = sv.samples.size();
    
    char b[200];
    string s;
    
    // format 
	size_t NFMT=sv.SVCF.FMT.size();
    string FMT=sv.SVCF.FMT[0].Id;
    for (int f=1; f<int(NFMT); f++) {
		FMT=FMT+":"+sv.SVCF.FMT[f].Id;
    }
	
    // events
    list<C_SVV1>::iterator i;
    for(i=sv.evt.begin(); i != sv.evt.end(); ++i) {
        output << (*i);
        
        output << "\t" << FMT ;
        
        for (int ns=0; ns<int(NS); ns++) {
            s=sv.samples[ns];
            // double cn = 2.0*double((*i).SampleMap[s].NR)/((*i).SampleMap[s].ER+0.01);
            sprintf(b,"\t%d:%d:%d:%d:%d:%.1f",(*i).SampleMap[s].N,(*i).SampleMap[s].NN,(*i).SampleMap[s].N5,(*i).SampleMap[s].N3,(*i).SampleMap[s].NR,(*i).SampleMap[s].ER); //cn);
            s = b;
            output << s;
        }
        output  << endl;  
    }  
    return output;
}


/*
 ostream &operator<<(ostream &output, C_SVV & sv)
 {
 output << sv.typeName << "\t " << sv.setName << "\t " << sv.contigName << endl;
 int Nevt = sv.evt.size();
 output << Nevt << endl;
 char b[200];
 sprintf(b,"%10s %10s %4s %4s", "pos","length","anc","q");
 string s = b;
 output << s;
 sprintf(b," %10s %10s %10s %10s %10s","N5","p5","std5","low5","high5");
 s = b;
 output << s;
 sprintf(b," %10s %10s %10s %10s","lm5","lmstd5","lmlow5","lmhigh5");
 s = b;
 output << s;
 sprintf(b," %10s %10s %10s %10s %10s","N3","p3","std3","low3","high3");
 s = b;
 output << s;
 sprintf(b," %10s %10s %10s %10s","lm3","lmstd3","lmlow3","lmhigh3");
 s = b;
 output << s;
 sprintf(b,"%10s %10s %10s %10s", "Nfrag","NfragPre","NFragPost","NfragExp");
 s = b;
 output << s;
 sprintf(b," %10s %10s %10s %10s", "p5a","p5b","p3a","p3b");
 s = b;
 output << s; 
 sprintf(b," %10s %10s %10s %10s %10s %10s %10s", "Nread","p0","p1","Nsite","Nrepeat","eN","CN");
 s = b;
 output << s;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s", "posU","lenU","Qoutlier","Qaberr","Source5","Source3","Merge5","Merge3");
 s = b;
 output << s ;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s", "a5","a3","Nr5","Nr3","N5","N3");
 s = b;
 output << s << endl;
 
 
 
 //
 list<C_SVV1>::iterator i;
 for(i=sv.evt.begin(); i != sv.evt.end(); ++i) {
 output << (*i) << endl;  
 }
 
 
 return output;
 }
 */

void C_SVV::print(string & outfilename)
{
    // format: optimized to loadFragments.m matlab script  
    // open output binary file. bomb if unable to open
    fstream output(outfilename.c_str(), ios::out);
    output << *this;
}


//------------------------------------------------------------------------------
// Cross link  F & R event class
//------------------------------------------------------------------------------
C_SVX1::C_SVX1() {
    type="cross";
    pos=0;
    anchor=0;
    length=0;
    q=0;
    p5[0]=0;
    p5[1]=0;
    p3[0]=0;
    p3[1]=0;
    //
    anchor2=0;    
    p5a2[0]=0;
    p3a2[0]=0;  
    p5a2[1]=0;
    p3a2[1]=0; 
    // 
    copynumber=0;
    cross5.clear();
    cross3.clear();
    NfragCluster[0]=0;
    NfragCluster[1]=0;
    NfragCov=0;
    NfragCovOut[0]=0;
    NfragCovOut[1]=0;
    NfragCovExp=0;
    //
    posU=0;
    lenU=0;
    qOutlier=0;
    merge5=false;
    merge3=false;
    a5[0]=0;
    a3[0]=0;
    a5[1]=0;
    a3[1]=0;
    
    SampleMap.clear();
    id=0;
    
}


ostream &operator<<(ostream &output, C_SVX1 & e1)
{
    //----------------------------------------------------------------------------
    //  text file position coordinates are 1 based (not c convention)
    //----------------------------------------------------------------------------  
    char b [200];
    string s=e1.type;
    sprintf(b,"%d\t%d\t%d\t.\t%s\t%d\t0\t", e1.anchor, e1.pos+1,e1.id, s.c_str(),e1.q);
    s = b;
    output << s;
    
    sprintf(b,"L=%d;NF=%d;UP=%d;UL=%d;", e1.length, int(e1.cls5.inp.size())+int(e1.cls3.inp.size()), e1.posU,e1.lenU);
    s = b;
    output << s;
    
    sprintf(b,"P2=%d;L2=%d;A2=%d;", e1.pos2+1,e1.length2,e1.anchor2);
    s = b;
    output << s;
    
    int PA=(e1.NfragCluster[0]>0? e1.p5[0]-e1.pos:0);
    int PB=(e1.NfragCluster[0]>0? e1.p5[1]-e1.pos:0);
    int PC=(e1.NfragCluster[1]>0? e1.p3[0]-e1.pos:0);
    int PD=(e1.NfragCluster[1]>0? e1.p3[1]-e1.pos:0);
    
    sprintf(b,"N5=%d;N3=%d;PA=%d;PB=%d;PC=%d;PD=%d;",e1.NfragCluster[0],e1.NfragCluster[1],PA,PB,PC,PD);
    s = b;
    output << s;
    
    sprintf(b,"NR=%d;ER=%.1f;MR=%d;", e1.cov.N,e1.cov.eN,(e1.merge5|e1.merge3)?1:0);
    s = b;
    output << s;  
    
    sprintf(b,"A5=%.2f;A3=%.2f;", e1.a5[0],e1.a3[0]);
    s = b;
    output << s;  
    
    return output;
}

/*
 ostream &operator<<(ostream &output, C_SVX1 & e1)
 {
 //----------------------------------------------------------------------------
 //  text file position coordinates are 1 based (not c convention)
 //----------------------------------------------------------------------------  
 char b [200];
 sprintf(b,"%10d\t%10d\t%4d\t%4d", e1.pos+1,e1.length,e1.anchor,e1.q);
 string s = b;
 output << s;
 
 sprintf(b,"%10d\t%10d\t%4d", e1.pos2+1,e1.length2,e1.anchor2);
 s = b;
 output << s;
 
 C_cluster2d_element1 c1 = e1.cls5;
 sprintf(b,"\t%10d\t%10.0f\t%10.1f\t%10.0f\t%10.0f", int(c1.inp.size()),c1.mean[0],c1.std[0],c1.low[0],c1.high[0]);
 s = b;
 output << s;
 sprintf(b,"\t%10.0f\t%10.1f\t%10.0f\t%10.0f", c1.mean[1],c1.std[1],c1.low[1],c1.high[1]);
 s = b;
 output << s;
 c1 = e1.cls3;
 sprintf(b,"\t%10d\t%10.0f\t%10.1f\t%10.0f\t%10.0f", int(c1.inp.size()),c1.mean[0],c1.std[0],c1.low[0],c1.high[0]);
 s = b;
 output << s;
 sprintf(b,"\t%10.0f\t%10.1f\t%10.0f\t%10.0f", c1.mean[1],c1.std[1],c1.low[1],c1.high[1]);
 s = b;
 output << s;
 
 sprintf(b,"\t%10d\t%10d\t%10d\t%10d", e1.NfragCov,e1.NfragCovOut[0],e1.NfragCovOut[1],e1.NfragCovExp);
 s = b;
 output << s;
 
 sprintf(b,"\t%10d\t%10d\t%10d\t%10d", e1.p5[0]+1,e1.p5[1]+1,e1.p3[0]+1,e1.p3[1]+1);
 s = b;
 output << s;  
 
 sprintf(b,"\t%10d\t%10d\t%10d\t%10d", e1.p5a2[0]+1,e1.p5a2[1]+1,e1.p3a2[0]+1,e1.p3a2[1]+1);
 s = b;
 output << s;  
 
 C_SVcoverage1 q1=e1.cov;
 sprintf(b,"\t%10d\t%10d\t%10d\t%10d\t%10d\t%10.1f\t%10d", q1.N,q1.p0+1,q1.p1+1,q1.Nsite,q1.Nrepeat,q1.eN,e1.copynumber);
 s = b;
 output << s;
 
 sprintf(b,"\t%10d\t%10d\t%10d", e1.posU,e1.lenU,e1.qOutlier);
 s = b;
 output << s;
 
 // iterate over read groups in 5' end;
 c1 = e1.cls5;
 sprintf(b,"\t%d", c1.N);
 s = b;
 std::map<unsigned int,unsigned int, less<unsigned int> >::iterator it;
 for ( it=e1.ReadGroupMap5.begin() ; it != e1.ReadGroupMap5.end(); it++ ) {
 unsigned int RGC1 = (*it).first;
 unsigned int NF1  = e1.ReadGroupMap5[RGC1];
 sprintf(b,"|%u:%u", NF1,RGC1);
 s=s+b;
 }
 output << s;
 
 // iterate over read groups in 3' end;
 c1 = e1.cls3;
 sprintf(b,"\t%d", c1.N);
 s = b;
 for ( it=e1.ReadGroupMap3.begin() ; it != e1.ReadGroupMap3.end(); it++ ) {
 unsigned int RGC1 = (*it).first;
 unsigned int NF1  = e1.ReadGroupMap3[RGC1];
 sprintf(b,"|%u:%u", NF1,RGC1);
 s=s+b;
 }
 output << s;
 
 // merge
 sprintf(b,"\t%10d\t%10d", (e1.merge5)?1:0,(e1.merge3)?1:0);
 s = b;
 output << s;
 
 // end alignability
 float a5=(e1.a5[0]>e1.a5[1]? e1.a5[0]: e1.a5[1] );
 float a3=(e1.a3[0]>e1.a3[1]? e1.a3[0]: e1.a3[1] );
 sprintf(b,"\t%10.3f\t%10.3f\t%10d\t%10d\t%10d\t%10d",a5,a3,e1.cov5.Nrepeat,e1.cov3.Nrepeat,e1.cov5.N,e1.cov3.N);
 s = b;
 output << s;
 
 return output;
 }
 */

C_SVX1& C_SVX1::operator=(const C_SVX1 &rhs)
{
    this->pos = rhs.pos;
    this->anchor = rhs.anchor;
    this->length = rhs.length;
    this->pos2 = rhs.pos2;
    this->anchor2 = rhs.anchor2;
    this->length2 = rhs.length2;
    this->q = rhs.q;
    this->p5[0] = rhs.p5[0];
    this->p5[1] = rhs.p5[1];
    this->p3[0] = rhs.p3[0];
    this->p3[1] = rhs.p3[1];
    this->p5a2[0] = rhs.p5a2[0];
    this->p5a2[1] = rhs.p5a2[1];
    this->p3a2[0] = rhs.p3a2[0];
    this->p3a2[1] = rhs.p3a2[1];
    this->copynumber = rhs.copynumber;
    this->cov    = rhs.cov;  
    this->cls5   = rhs.cls5;  
    this->cls3   = rhs.cls3;  
    this->cross5 = rhs.cross5;  
    this->cross3 = rhs.cross3;  
    this->NfragCluster[0]= rhs.NfragCluster[0];
    this->NfragCluster[1]= rhs.NfragCluster[1];
    this->NfragCov       = rhs.NfragCov;
    this->NfragCovOut[0] = rhs.NfragCovOut[0];
    this->NfragCovOut[1] = rhs.NfragCovOut[1];
    this->NfragCovExp    = rhs.NfragCovExp;
    //
    this->cov5  = rhs.cov5;
    this->cov3  = rhs.cov3;
    this->a5[0] = rhs.a5[0];
    this->a5[1] = rhs.a5[1];
    this->a3[0] = rhs.a3[0];
    this->a3[1] = rhs.a3[1];
    this->posU  = rhs.posU;
    this->lenU  = rhs.lenU;
    this->qOutlier     = rhs.qOutlier;
    this->ReadGroupMap5= rhs.ReadGroupMap5;
    this->ReadGroupMap3= rhs.ReadGroupMap3;
    this->merge5 = rhs.merge5;
    this->merge3 = rhs.merge3;
    this->id             = rhs.id;
    this->type           = rhs.type;
    this->SampleMap      = rhs.SampleMap;
    
    return *this;
}

int C_SVX1::operator==(const C_SVX1 &rhs) const
{
    if( this->pos != rhs.pos) return 0;
    if( this->anchor != rhs.anchor) return 0;
    if( this->length != rhs.length) return 0;
    if( this->pos2 != rhs.pos2) return 0;
    if( this->anchor2 != rhs.anchor2) return 0;
    if( this->length2 != rhs.length2) return 0;
    if( this->p5[0] != rhs.p5[0]) return 0;
    if( this->p5[1] != rhs.p5[1]) return 0;
    if( this->p3[0] != rhs.p3[0]) return 0;
    if( this->p3[1] != rhs.p3[1]) return 0;
    if( this->p5a2[0] != rhs.p5a2[0]) return 0;
    if( this->p5a2[1] != rhs.p5a2[1]) return 0;
    if( this->p3a2[0] != rhs.p3a2[0]) return 0;
    if( this->p3a2[1] != rhs.p3a2[1]) return 0;
    if( this->q != rhs.q) return 0;
    if (this->copynumber != rhs.copynumber) return 0;
    if (this->NfragCov   != rhs.NfragCov) return 0;
    if (this->cls5.inp.size()!= rhs.cls5.inp.size()) return 0;
    if (this->cls3.inp.size()!= rhs.cls3.inp.size()) return 0;
    
    if( this->posU != rhs.posU) return 0;
    if( this->lenU != rhs.lenU) return 0;
    if( this->qOutlier != rhs.qOutlier) return 0;
    if (this->merge5 != rhs.merge5) return 0;
    if (this->merge3 != rhs.merge3) return 0;
    if (this->ReadGroupMap5.size() != rhs.ReadGroupMap5.size()) return 0;
    if (this->ReadGroupMap3.size() != rhs.ReadGroupMap3.size()) return 0;
    if( this->a5[0] != rhs.a5[0]) return 0;
    if( this->a3[0] != rhs.a3[0]) return 0;
    if( this->a5[1] != rhs.a5[1]) return 0;
    if( this->a3[1] != rhs.a3[1]) return 0;
    
    return 1;
}

int C_SVX1::operator<(const C_SVX1 &rhs) const
{
    if( this->anchor < rhs.anchor ) return 1;
    if( this->anchor > rhs.anchor ) return 0;
    if( this->pos < rhs.pos ) return 1;
    if( this->pos > rhs.pos ) return 0;
    if( this->anchor2 < rhs.anchor2 ) return 1;
    if( this->anchor2 > rhs.anchor2 ) return 0;
    if( this->pos2 < rhs.pos2 ) return 1;
    if( this->pos2 > rhs.pos2 ) return 0;
    if( this->length < rhs.length ) return 1;
    if( this->length > rhs.length ) return 0;
    if( this->length2 < rhs.length2 ) return 1;
    if( this->length2 > rhs.length2 ) return 0;
    if( this->cls5.inp.size() < rhs.cls5.inp.size() ) return 1;
    if( this->cls5.inp.size() > rhs.cls5.inp.size() ) return 0;
    if( this->cls3.inp.size() < rhs.cls3.inp.size() ) return 1;
    if( this->cls3.inp.size() > rhs.cls3.inp.size() ) return 0;
    return 0;
}


//------------------------------------------------------------------------------
//  F & R clustered event class (inversions)
//------------------------------------------------------------------------------ 
C_SVX::C_SVX() {
    typeName="x";
    contigName="";
    setName="";
    evt.clear(); 
}

C_SVX::C_SVX(C_contig  & c1,  RunControlParameters & par) {
    typeName="x";
    contigName =c1.getContigName();
    setName =c1.setName;
    evt.clear(); 
}

void C_SVX::finalize(C_contig  & contig, C_libraries & libraries, RunControlParameters & pars) {
    
    //----------------------------------------------------------------------------
    // sample by sample genotype info
    //----------------------------------------------------------------------------
    genotype(contig,libraries,pars);
    
    //----------------------------------------------------------------------------
    // loop over events / clusters within events to calculate:
    //----------------------------------------------------------------------------
    
    // check for events, libraries
    if (evt.size()==0) { return;}
    if (libraries.libmap.size()==0) { return;}
    
    list<C_SVX1>::iterator ie,ie1,ie2;
    // first event
    ie1=evt.begin();
    // last event
    ie2=evt.end();  
    // loop over events
    int ne=0;
    
    // max a5 or a3 alignability metrics
    double aMax=0;
    
    for ( ie=ie1 ; ie != ie2; ++ie ) {
        ne++;
        
        // id
        (*ie).id=ne;
        
        // find max alignability scale
        if ((*ie).a5[0]>aMax) { aMax=(*ie).a5[0];}
        if ((*ie).a3[0]>aMax) { aMax=(*ie).a3[0];}
        if ((*ie).a5[1]>aMax) { aMax=(*ie).a5[1];}
        if ((*ie).a3[1]>aMax) { aMax=(*ie).a3[1];}
        
    }
    
    // normalize alignability to <=1.0
    for ( ie=evt.begin(); ie != evt.end(); ++ie ) {
        (*ie).a5[0]=(*ie).a5[0]/aMax;
        (*ie).a3[0]=(*ie).a3[0]/aMax;
        (*ie).a5[1]=(*ie).a5[1]/aMax;
        (*ie).a3[1]=(*ie).a3[1]/aMax;
    }
    
} 


void C_SVX::genotype(C_contig  & contig, C_libraries & libraries, RunControlParameters & pars)
{
    // check for events, libraries
    if (evt.size()==0) { return;}
    if (libraries.libmap.size()==0) { return;}
    
    map<string, double, less<string> >  rX = libraries.readFractionSamples();
    
    list<C_SVX1>::iterator ie,ie1,ie2;
    //C_SVR1 e1;
    unsigned int ReadGroupCode1;
    string SAM;
    int NF1;
    
    ie1=evt.begin();
    
    // know when to stop
    ie2=evt.end();  
    
    list<C_localpair>::iterator i;
    // get Minimum Q mapping value for unique read
    int Qmin  = pars.getQmin();
    // loop over pairs
    int np=0;
    for(i=contig.localpairs.begin(); i != contig.localpairs.end(); ++i) {
        np++;
        // skip low mapping quality fragments
        if (((*i).q1<Qmin)&((*i).q2<Qmin)) {
            continue;
        }
        // fragment length limits for this library
        ReadGroupCode1 = (*i).ReadGroupCode;  
        //int LM1=libraries.libmap[ReadGroupCode1].LM;
        int LMhigh1=libraries.libmap[ReadGroupCode1].LMhigh;
        int LMlow1=libraries.libmap[ReadGroupCode1].LMlow;
        
        // fragment length
        int lm = (*i).lm;
        // demand usual fragment
        if ((lm<LMlow1)|(lm>LMhigh1)) {continue; }
        
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        
        // add only the non-read part of the fragment
        int p0 = (*i).pos;
        int p1 = (*i).pos+(*i).len1;
        int p2 = p1+(*i).lm-(*i).len2;
        
        // loop over events
        int ne=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            
            ne++;
            
            // two spanning windows eP1a-eP1b and eP2b-eP2b  
            int eP1a=(*ie).pos-(*ie).posU/2;     
            int eP1b=int((*ie).pos)+int((*ie).posU/2);
            int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
            int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
            
            // check for reads starts inside event window eP1b-eP2a
            if ((p0>=eP1b)&(p0<=eP2a)&((*i).q1>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            // other end of fragment
            if ((p2>=eP1a)&(p2<=eP2a)&((*i).q1>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            
            // beginning of event well beyond this fragment
            if ((eP1a-1000)>p2) { break;}        
            
            // end of event not yet reaching this fragment
            if (eP2b<p1) { continue;}       
            
            // select spanning fragments - both above Qmin 
            if (((*i).q1<Qmin)|((*i).q2<Qmin)) {
                continue;
            } 
            
            // spanning 5' break
            if ((p1<eP1a)&(p2>eP1b)) {
                (*ie).SampleMap[SAM].NN++;
            }
            // spanning 3' break
            if ((p1<eP2a)&(p2>eP2b)) {
                (*ie).SampleMap[SAM].NN++;
            }
            
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
        
    }
    
    ie1=evt.begin();
    
    std::map<unsigned int,unsigned int, less<unsigned int> >::iterator it;  
    // iterate over read groups in 5' end;
    for ( ie=ie1 ; ie != ie2; ++ie ) {
        //e1 = *ie;
        for ( it=(*ie).ReadGroupMap5.begin() ; it != (*ie).ReadGroupMap5.end(); it++ ) {
            ReadGroupCode1 = (*it).first;
            NF1  = (*ie).ReadGroupMap5[ReadGroupCode1];
            SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
            (*ie).SampleMap[SAM].N5+=NF1;
            (*ie).SampleMap[SAM].N+=NF1;
        }     
        for ( it=(*ie).ReadGroupMap3.begin() ; it != (*ie).ReadGroupMap3.end(); it++ ) {
            ReadGroupCode1 = (*it).first;
            NF1  = (*ie).ReadGroupMap3[ReadGroupCode1];
            SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
            (*ie).SampleMap[SAM].N3+=NF1;
            (*ie).SampleMap[SAM].N+=NF1;
        }      
        
        // loop over libraries for estimated number of reads 
        //map<string, double, less<string> > rX = libraries.readFractionSamples();
        map<string, double, less<string> >::iterator irX;
        for ( irX=rX.begin() ; irX != rX.end(); irX++ ) {
            SAM=(*irX).first;
            (*ie).SampleMap[SAM].ER=double((*ie).cov.eN)*rX[SAM];
        }
    }
    
    // loop over single ends
    
    // fill depth from cross pairs (not so many...)
    int N = contig.crosspairs.size();
    list<C_crosspair>::iterator ix;
    ie1=evt.begin();
    for(ix=contig.crosspairs.begin(); ix != contig.crosspairs.end(); ++ix) {
        // skip low mapping quality fragments
        if ((*ix).read[0].q<Qmin) {continue;}
        // fragment length limits for this library
        ReadGroupCode1 = (*ix).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        
        // read start
        int p1=int((*ix).read[0].pos);
        
        // loop over events
        int nx=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nx++;
            
            int eP1a=(*ie).pos-(*ie).posU/2;     
            int eP1b=int((*ie).pos)+int((*ie).posU/2);
            int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
            int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
            
            // beginning of event well beyond this fragment
            if ((eP1a-1000)>p1) { break;}        
            
            // end of event not yet reaching this fragment
            if (eP2b<p1) { continue;}       
            
            // check for reads starts inside event window eP1b-eP2a
            if ((p1>=eP1b)&(p1<=eP2a)&((*ix).read[0].q>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    // fill depth from dangle starts 
    N = contig.dangle.size();
    list<C_singleEnd>::iterator is; 
    for(is=contig.dangle.begin(); is != contig.dangle.end(); ++is) {
        // skip low mapping quality fragments
        if ((*is).q<Qmin) {continue;}
        ReadGroupCode1 = (*is).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        
        // read start
        int p1=int((*is).pos);
        
        // loop over events
        int nd=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nd++;
            
            int eP1a=(*ie).pos-(*ie).posU/2;     
            int eP1b=int((*ie).pos)+int((*ie).posU/2);
            int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
            int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
            
            // beginning of event well beyond this fragment
            if ((eP1a-1000)>p1) { break;}        
            
            // end of event not yet reaching this fragment
            if (eP2b<p1) { continue;}       
            
            // check for reads starts inside event window eP1b-eP2a
            if ((p1>=eP1b)&(p1<=eP2a)&((*is).q>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    // fill depth from umpairs 
    N = contig.umpairs.size();
    list<C_umpair>::iterator iu;
    
    for(iu=contig.umpairs.begin(); iu != contig.umpairs.end(); ++iu) {
        if ((*iu).read[0].q<Qmin) {continue;}
        ReadGroupCode1 = (*iu).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        
        // read start
        int p1=int((*iu).read[0].pos);
        
        // loop over events
        int nd=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nd++;
            
            int eP1a=(*ie).pos-(*ie).posU/2;     
            int eP1b=int((*ie).pos)+int((*ie).posU/2);
            int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
            int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
            
            // beginning of event well beyond this fragment
            if ((eP1a-1000)>p1) { break;}        
            
            // end of event not yet reaching this fragment
            if (eP2b<p1) { continue;}       
            
            // check for reads starts inside event window eP1b-eP2a
            if ((p1>=eP1b)&(p1<=eP2a)&((*iu).read[0].q>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
            
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    // fill depth from unique ends 
    N = contig.singleton.size();
    for(is=contig.singleton.begin(); is != contig.singleton.end(); ++is) {
        if ((*is).q<Qmin) {continue;}
        ReadGroupCode1 = (*is).ReadGroupCode;  
        SAM = libraries.libmap[ReadGroupCode1].Info.SampleName;
        
        // read start
        int p1=int((*is).pos);
        
        // loop over events
        int nd=0;
        for ( ie=ie1 ; ie != ie2; ++ie ) {
            nd++;
            int eP1a=(*ie).pos-(*ie).posU/2;     
            int eP1b=int((*ie).pos)+int((*ie).posU/2);
            int eP2a=int((*ie).pos)+int((*ie).length)-int((*ie).posU/2);
            int eP2b=int((*ie).pos)+int((*ie).length)+int((*ie).posU/2);
            
            // beginning of event well beyond this fragment
            if ((eP1a-1000)>p1) { break;}        
            
            // end of event not yet reaching this fragment
            if (eP2b<p1) { continue;}       
            
            // check for reads starts inside event window eP1b-eP2a
            if ((p1>=eP1b)&(p1<=eP2a)&((*is).q>=Qmin)) {       
                (*ie).SampleMap[SAM].NR++;
            }
        }
        // start next loop over events where this one left off - 3
        if (ie!=ie1)  ie1=ie--;
        if (ie1!=evt.begin())  ie1--;
        if (ie1!=evt.begin())  ie1--;
    }
    
    
}

ostream &operator<<(ostream &output, C_SVX & sv)
{
    
    // SVCF header
    output  << sv.SVCF;
    
    size_t NS = sv.samples.size();
    
    char b[200];
    string s;
    
    // format 
	size_t NFMT=sv.SVCF.FMT.size();
    string FMT=sv.SVCF.FMT[0].Id;
    for (int f=1; f<int(NFMT); f++) {
		FMT=FMT+":"+sv.SVCF.FMT[f].Id;
    }
	
    // events
    list<C_SVX1>::iterator i;
    for(i=sv.evt.begin(); i != sv.evt.end(); ++i) {
        output << (*i);
        
        output << "\t" << FMT ;
        
        for (int ns=0; ns<int(NS); ns++) {
            s=sv.samples[ns];
            // double cn = 2.0*double((*i).SampleMap[s].NR)/((*i).SampleMap[s].ER+0.01);
            sprintf(b,"\t%d:%d:%d:%d:%d:%.1f",(*i).SampleMap[s].N,(*i).SampleMap[s].NN,(*i).SampleMap[s].N5,(*i).SampleMap[s].N3,(*i).SampleMap[s].NR,(*i).SampleMap[s].ER); //cn);
            s = b;
            output << s;
        }
        output  << endl;  
    }  
    return output;
}

/*
 ostream &operator<<(ostream &output, C_SVX & sv)
 {
 output << sv.typeName << "\t " << sv.setName << "\t " << sv.contigName << endl;
 int Nevt = sv.evt.size();
 output << Nevt << endl;
 char b[200];
 sprintf(b,"%10s %10s %4s %4s", "pos","length","anc","q");
 string s = b;
 output << s;
 sprintf(b,"%10s %10s %4s", "pos2","length2","anc2");
 s = b;
 output << s;
 sprintf(b," %10s %10s %10s %10s %10s","N5","p5","std5","low5","high5");
 s = b;
 output << s;
 sprintf(b," %10s %10s %10s %10s","p5a2","std5a2","low5a2","high5a2");
 s = b;
 output << s;
 sprintf(b," %10s %10s %10s %10s %10s","N3","p3","std3","low3","high3");
 s = b;
 output << s;
 sprintf(b," %10s %10s %10s %10s","p3a2","std3a2","low3a2","high3a2");
 s = b;
 output << s;
 sprintf(b,"%10s %10s %10s %10s", "Nfrag","NfragPre","NFragPost","NfragExp");
 s = b;
 output << s;
 sprintf(b," %10s %10s %10s %10s", "p5a","p5b","p3a","p3b");
 s = b;
 output << s; 
 sprintf(b," %10s %10s %10s %10s", "p5a2","p5b2","p3a2","p3b2");
 s = b;
 output << s; 
 sprintf(b," %10s %10s %10s %10s %10s %10s %10s", "Nread","p0","p1","Nsite","Nrepeat","eN","CN");
 s = b;
 output << s;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s", "posU","lenU","Qoutlier","Source5","Source3","Merge5","Merge3");
 s = b;
 output << s ;
 sprintf(b,"\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s", "a5","a3","Nr5","Nr3","N5","N3");
 s = b;
 output << s << endl;
 
 //
 list<C_SVX1>::iterator i;
 for(i=sv.evt.begin(); i != sv.evt.end(); ++i) {
 output << (*i) << endl;  
 }
 
 return output;
 }
 */

void C_SVX::print(string & outfilename)
{
    // format: optimized to loadFragments.m matlab script  
    // open output binary file. bomb if unable to open
    fstream output(outfilename.c_str(), ios::out);
    output << *this;
}



//------------------------------------------------------------------------------
//  Clustering class (not for Retro mob insertions) 
//------------------------------------------------------------------------------ 
C_SpannerCluster::C_SpannerCluster(C_contig & c1, C_libraries & libs1, RunControlParameters & pars1) {
    //
    pars = pars1;
    contigName=c1.getContigName();
    setName=c1.setName; 
    int Qmin = pars.getQmin();
    
    //library info
    libraries = libs1;
    
    L = c1.Length;   
    double Npair = c1.localpairs.size();
    PairDensity = Npair/L;
    double Nend = c1.dangle.size();
    EndDensity = 0.5*Nend/L;
    //-------------------------------------------------------------------------
    // select abberant pairs  -> longpairs, shortpairs, invert5, invert3
    // length too big, too small, or flipped orientation
    //-------------------------------------------------------------------------
    selectPairs(c1,Qmin);
    //-------------------------------------------------------------------------
    //inter-chromosomal linked fragments
    //-------------------------------------------------------------------------
    selectCross(c1,Qmin);
    //-------------------------------------------------------------------------
    // set neighborhood window to some scale of fragment length width 
    // for breakpoint spanning pairs
    //-------------------------------------------------------------------------
    // HistObj h = pars.getFragHist(); // 12 July 2009 - libraries not FragHist
    // double window = h.std/2.0; // original
    // double window = 3*h.std;   // X Oct 2008
    // double window = 4.0*h.std;    // 5 Nov 2008
    // double window = 4.0*h.std;    // 5 Nov 2008
    double fwindow = 1.0;    // 12 July 2009 - relative to 
    pars.setClusteringLength(int(fwindow));
    //-------------------------------------------------------------------------
    // long /short LM clusters for inversion detection
    //-------------------------------------------------------------------------
    // double fwindow2[2] = {fwindow,fwindow};
    // double fwindow2[2] = {fwindow*,fwindow/2.0};
    double fwindow2[2] = {fwindow*1.25,fwindow/2.0};
    //-------------------------------------------------------------------------
    // both forward mapped ~ '>' reads 5'
    //-------------------------------------------------------------------------
    vector<vector<double> > x2;
    x2.clear();    
    vector<vector<double> > wx2;
    wx2.clear();    
    vector<int> ip2;
    ip2.clear();    
    
    //cout << "clustering " << endl;
    
    int N =  makepairX(cross5,'F',x2,wx2,ip2);    
    if (N>0) {
        C_NNcluster2d q(x2, wx2, ip2, fwindow2,"cross5",c1.setName,c1.getContigName());  
        //q.mergeClusters();   
        q.cleanClusters();     
        cross5c = q; 
    }
    //cout << "\t" << setw(15) << cross5c.typeName << "\t " << cross5c.NC << endl;
    
    x2.clear();    
    wx2.clear();    
    ip2.clear();    
    
    N =  makepairX(cross3,'R',x2,wx2,ip2);    
    if (N>0) {
        C_NNcluster2d q(x2, wx2, ip2, fwindow2,"cross3",c1.setName,c1.getContigName());  
        //q.mergeClusters();   
        q.cleanClusters();     
        cross3c = q; 
    }
    //cout << "\t" << setw(15) << cross3c.typeName << "\t " << cross3c.NC << endl;
    
    x2.clear();    
    wx2.clear();    
    ip2.clear();    
    
    N =  makepairP(longpair,'-',x2,wx2,ip2);    
    if (N>0) {
        C_NNcluster2d q(x2, wx2, ip2, fwindow2,"longpair",c1.setName,c1.getContigName());  
        //q.mergeClusters();   
        q.cleanClusters();     
        longpairc = q; 
    }
    //cout << "\t" << setw(15) << longpairc.typeName << "\t " << longpairc.NC << endl;
    x2.clear();    
    wx2.clear();    
    ip2.clear();    
    
    N =  makepairP(shortpair,'-',x2,wx2,ip2);
    if (N>0) {
        C_NNcluster2d q(x2, wx2, ip2, fwindow2,"shortpair",c1.setName,c1.getContigName());  
        q.cleanClusters();
        shortpairc = q; 
    }
    //cout << "\t" << setw(15) << shortpairc.typeName << "\t " << shortpairc.NC << endl;
    x2.clear();
    wx2.clear();
    ip2.clear();    
    
    //-------------------------------------------------------------------------
    // inverted clusters for inversion detection
    //-------------------------------------------------------------------------
    
    // double window = h.std/2.0; // original
    // double window = 3*h.std;   // X Oct 2008
    // double window = 4.0*h.std;    // 5 Nov 2008
    //window = h.mean;    // 23 April  Nov 2009 - thanks to EEE et al. 
    //pars.setClusteringLength(int(window));
    
    
    N =  makepairP(invert5,'>',x2,wx2,ip2);
    if (N>0) {
        C_NNcluster2d q(x2, wx2, ip2, fwindow2,"invert5",c1.setName,c1.getContigName());  
        q.cleanClusters();
        invert5c = q; 
    }
    //cout << "\t" << setw(15)  << "invert5" << "\t " << invert5c.NC << endl;
    //-------------------------------------------------------------------------
    // both reverse mapped ~ '<' reads 3'
    //-------------------------------------------------------------------------
    x2.clear();
    wx2.clear();
    ip2.clear();    
    
    N =  makepairP(invert3,'<',x2,wx2,ip2);
    if (N>0) {
        C_NNcluster2d q(x2, wx2, ip2, fwindow2,"invert3",c1.setName,c1.getContigName());  
        q.cleanClusters();
        invert3c = q; 
    }
    //cout << "\t" << setw(15) << "invert3" << "\t " << invert3c.NC << endl;
    
    //-------------------------------------------------------------------------
    // set neighborhood window to some scale of fragment length  
    // for fragments dangling into insert dna
    //-------------------------------------------------------------------------
    double fwindow1 = 2.0;  //  optimize ? 
    pars.setClusteringLengthRetro(int(fwindow1));
    //-------------------------------------------------------------------------
    // dangling end clusters for insertion detection
    //-------------------------------------------------------------------------    
    /* 
     vector<double> x;
     N =  makeDangleX(c1,'5',x);  
     cout << "dangle5  " << N << "\t " << x.size() << endl;   
     if (N>0) {
     C_NNcluster1d q(x, window1,"dangle5",c1.setName,c1.getContigName());  
     dangle5c = q; 
     }
     cout << dangle5c.typeName << "\t " << dangle5c.NC << endl;
     //-------------------------------------------------------------------------
     // other end of the insertion
     //-------------------------------------------------------------------------
     x.clear();
     N =  makeDangleX(c1,'3',x);    
     if (N>0) {
     C_NNcluster1d q(x, window1,"dangle3",c1.setName,c1.getContigName());  
     dangle3c = q; 
     }
     cout << dangle3c.typeName << "\t " << dangle3c.NC << endl;
     */
}



int C_SpannerCluster::makeDangleX(C_contig & c1, char e, vector<double> & x1) {
    //list<C_readmap>* d1;
    char S='F';
    if (e=='5') {
        S='F';
        //d1 = & c1.dangleEnd;
    } else if (e=='3') {
        S='R';
        //d1 = & c1.dangleStart;
    } else { 
        cerr << " bad arg to makeDangleX: " << e << endl;
        exit(1);
    }
    int N =  c1.dangle.size();    
    x1.resize(N,0);        
    int j=0;
    if (N>0) {
        list<C_singleEnd>::iterator i;
        for(i=c1.dangle.begin(); i != c1.dangle.end(); ++i) {
            char s= (*i).sense;
            if (s==S) {   
                x1[j] = (*i).pos;
                if (e=='5') {  // add the read length for the 5' end 
                    x1[j] = x1[j]+(*i).len;
                }
                j++;
            }
        }
    }
    x1.resize(j+1);     
    return int(x1.size());
}

void C_SpannerCluster::selectPairs(C_contig & c1, int Qmin) {
    int N0 =  c1.localpairs.size();    
    invert5.clear();        
    invert3.clear();        
    longpair.clear();        
    shortpair.clear();
    list<C_localpair>::iterator i;
    if (N0>0) {
        for(i=c1.localpairs.begin(); i != c1.localpairs.end(); ++i) {
            
            //----------------------------------------------------------------------
            // only high quality mapped read pairs
            //----------------------------------------------------------------------        
            //if ((*i).q1>100) (*i).q1=0;
            //if ((*i).q2>100) (*i).q2=0;          
            if (int((*i).q1)<Qmin) continue;
            if (int((*i).q2)<Qmin) continue;
            
            char o1=(*i).orient;  
            //double p1= double((*i).pos);
            int lm1=(*i).lm;  
            
            // skip flipped libs (ie. unfixed jumping libs)
            int LF1 = libraries.libmap[(*i).ReadGroupCode].LM;
            if (LF1<0) continue;
            
            if (o1=='>') {
                invert5.push_back((*i));
            } else if (o1=='<') {
                invert3.push_back((*i));
            } else if (o1=='-') {
                // library based selection (7/19/2009)
                int LFhigh = libraries.libmap[(*i).ReadGroupCode].LMhigh;        
                int LFlow = libraries.libmap[(*i).ReadGroupCode].LMlow;
                

                //--------------------------------------------------------------------
                // many libraries have significant tail extending down to LR
                // use min of LR & LMlow; 
                //--------------------------------------------------------------------
                // 1000G library fix for too many dups 9/11/2009
                // int LRmin = libraries.libmap[(*i).ReadGroupCode].LRmin;
                // LFlow=(LFlow<LRmin? LFlow: LRmin);
                //
                if (lm1>LFhigh) {  //pars.getFragmentLengthHi()) {
                    longpair.push_back((*i));
                } else if (lm1< LFlow) {  //pars.getFragmentLengthLo()) {
                    shortpair.push_back((*i));
                } // only normal pairs remain...
            } else {
                cerr << " bad orientation in selectpairs " << (*i) << endl;
                exit(1);
            } 
        }
    }
}

void C_SpannerCluster::selectCross(C_contig & c1, int Qmin) {
    int N0 =  c1.crosspairs.size();    
    cross5.clear();        
    cross3.clear();        
    list<C_crosspair>::iterator i;
    if (N0>0) {
        for(i=c1.crosspairs.begin(); i != c1.crosspairs.end(); ++i) {
            
            //----------------------------------------------------------------------
            // only high quality mapped read pairs
            //----------------------------------------------------------------------        
            if ((*i).read[0].q>100) (*i).read[0].q=0;
            if ((*i).read[1].q>100) (*i).read[1].q=0;
            if (int((*i).read[0].q)<Qmin) continue;
            if (int((*i).read[1].q)<Qmin) continue;
            
            char s1=(*i).read[0].sense;  
            
            if (s1=='F') {
                cross5.push_back((*i));
            } else if (s1=='R') {
                cross3.push_back((*i));
            } else {
                cerr << " bad orientation in selectCross " << (*i) << endl;
                exit(1);
            } 
        }
    }
}


bool compare_vector (vector<double> & v1 , vector<double> & v2)
{
    size_t N=4;
    if (v1.size()!=v2.size()) N=(v1.size()>v2.size()?v2.size(): v1.size());
    for (int n=0; n<int(N); n++) {
        if (v1[n]!=v2[n]) {
            return (v1[n]<v2[n]);
        }
    }
    // should not get here
    return false;
}


//------------------------------------------------------------------------------
// cross span clusters
//------------------------------------------------------------------------------
int C_SpannerCluster::makepairX(vector<C_crosspair> & p1,  char sense
                                , vector<vector<double> >& x1, vector<vector<double> >& wx1, vector<int> & inp1) {
    
    int N0 =  p1.size();    
    x1.clear();          
    wx1.clear();          
    inp1.clear();
    
    C_vectorDouble x4(5,0);
    vector<double> x2(2,0);
    
    // use a list (for sorting...)  
    list<C_vectorDouble> xlist; 
    list<C_vectorDouble>::iterator it;
    
    if (N0>0) {
        for(int i=0; i< N0; i++) {
            char s1=p1[i].read[0].sense;  
            if (s1!=sense) {
                cerr << " wrong orientation in makepairP" << p1[i] << endl;
            }
            // index to p1 order
            x4[4]=i;
            
            // pos at  5' end of fragment
            x4[0]=double(p1[i].read[0].pos);
            x4[1]=double(p1[i].read[1].anchor*1e10+p1[i].read[1].pos);  
            
            // skip flipped libs (i.e. unfixed Jumping libs)
            int LF1=libraries.libmap[p1[i].ReadGroupCode].LMhigh;
            if (LF1<0) continue;
            
            // get clustering width from library info
            int LFlow = libraries.libmap[p1[i].ReadGroupCode].LMlow;
            int LFhigh = libraries.libmap[p1[i].ReadGroupCode].LMhigh;
            double W=double(LFhigh-LFlow);
            
            x4[2]=W;
            x4[3]=W;
            xlist.push_back(x4);
        }
    }
    
    //sort
    xlist.sort(compare_vector);
    
    // dump into x1 and wx1
    for (it=xlist.begin(); it!=xlist.end(); ++it) {
        x4=*it;
        x2[0]=x4[0];
        x2[1]=x4[1];
        x1.push_back(x2);
        x2[0]=x4[2];
        x2[1]=x4[3];
        wx1.push_back(x2);
        inp1.push_back(int(round(x4[4])));
        
    }
    return int(x1.size());
}


int C_SpannerCluster::makepairP(vector<C_localpair> & p1,  char orient
                                , vector<vector<double> >& x1, vector<vector<double> >& wx1, vector<int> & inp1) {
    
    int N0 =  p1.size();    
    x1.clear();          
    wx1.clear();          
    inp1.clear();
    
    C_vectorDouble x4(5,0);
    vector<double> x2(2,0);
    
    // use a list (for sorting...)  
    list<C_vectorDouble> xlist; 
    list<C_vectorDouble>::iterator it;
    
    if (N0>0) {
        for(int i=0; i< N0; i++) {
            char o1=p1[i].orient;  
            if (o1!=orient) {
                cerr << " wrong orientation in makepairP" << p1[i] << endl;
            }
            // index to p1 order
            x4[4]=i;
            
            // position at center of fragment
            // x2[0]=double(p1[i].pos);  5' end of fragment
            x4[0]=double(p1[i].pos+(p1[i].lm/2));
            
            // get average LF from library info record
            int LF1 = libraries.libmap[p1[i].ReadGroupCode].LM;
            
            //  skip libraries with negative fragment length 
            if (LF1<0) continue;
            
            x4[1]=double(p1[i].lm - LF1);  
            //x1.push_back(x);
            // get clustering width from library info
            double W1=LF1;  // deletions, duplications, inversions
            double W2=2*LF1;  // inversions
            
            if (orient=='-') { // deletions, insertions, duplications
                int LFlow = libraries.libmap[p1[i].ReadGroupCode].LMlow;
                int LFhigh = libraries.libmap[p1[i].ReadGroupCode].LMhigh;
                W2=double(LFhigh-LFlow);
            }
            
            x4[2]=W1;
            x4[3]=W2;
            //wx1.push_back(x);
            xlist.push_back(x4);
        }
    }
    
    //sort
    xlist.sort(compare_vector);
    
    // dump into x1 and wx1
    for (it=xlist.begin(); it!=xlist.end(); ++it) {
        x4=*it;
        x2[0]=x4[0];
        x2[1]=x4[1];
        x1.push_back(x2);
        x2[0]=x4[2];
        x2[1]=x4[3];
        wx1.push_back(x2);
        inp1.push_back(int(round(x4[4])));
        
    }
    return int(x1.size());
}




void C_SpannerCluster::writeall() {
    string area = pars.getOutputDir();
    int k = 1+area.find_last_of("/");
    string subdir = area.substr(k);
    if (area.size()>0) { area = area+"/";}
    string prefix = pars.getPrefix();
    if (prefix.size()>0) { prefix = prefix+".";}
    string sn1 = setName;
    // dont need redundant setName if subdirectory is already the setName
    // if (sn1==subdir) sn1="";
    if (sn1.size()>0) { sn1 = sn1+".";}
    string basename = area+prefix+sn1+contigName;
    if (setName==contigName) basename = area+prefix+contigName;
    
    
    // replace evil character "|" with benign "_" 
    size_t found=basename.find("|");
    while (found!=string::npos) {
        basename.replace(found,1,"_");
        found=basename.find("|");
    }
    
    if (dangle5c.cluster.size()>0) { 
        // cout << "write clusters for dangle 5' : " << sn1+contigName << endl;
        string fname = basename+".dangle5.cls.span";
        cout << "write " << fname << "\t " << dangle5c.cluster.size() << endl;
        dangle5c.write(fname);
        //write(fname,dangle5c,dangle5);
    }
    if (dangle3c.cluster.size()>0) { 
        // cout << "write clusters for dangle 3': " << sn1+contigName << endl;
        string fname = basename+".dangle3.cls.span";
        cout << "write " << fname << "\t " << dangle3c.cluster.size() << endl;
        dangle3c.write(fname);
        //write(fname,dangle3c,dangle3);
    }
    // invert clusters
    if (invert5c.cluster.size()>0) { 
        // cout << "write clusters for invert 5': " << sn1+contigName << endl;
        string fname = basename+".invert5.cls.span";
        cout << "write " << fname << "\t " << invert5c.cluster.size() << endl;
        //invert5c.write(fname);
        write(fname,invert5c,invert5);
    }
    if (invert3c.cluster.size()>0) { 
        //cout << "write clusters for invert 3': " << sn1+contigName << endl;
        string fname = basename+".invert3.cls.span";
        cout << "write " << fname << "\t " << invert3c.cluster.size() << endl;
        //invert3c.write(fname);
        write(fname,invert3c,invert3);
    }
    // long LM clusters
    if (longpairc.cluster.size()>0) { 
        //cout << "write clusters for long LM: " << sn1+contigName << endl;
        string fname = basename+".long.cls.span";
        cout << "write " << fname << "\t " << longpairc.cluster.size() << endl;
        //longpairc.write(fname);
        write(fname,longpairc,longpair);
    }
    if (shortpairc.cluster.size()>0) { 
        //cout << "write clusters for short LM: " << sn1+contigName << endl;
        string fname = basename+".short.cls.span";
        cout << "write " << fname << "\t " << shortpairc.cluster.size() << endl;
        //shortpairc.write(fname);
        write(fname,shortpairc,shortpair);
    }
    if (cross5c.cluster.size()>0) { 
        // cout << "write clusters for cross 5' : " << sn1+contigName << endl;
        string fname = basename+".cross5.cls.span";
        cout << "write " << fname << "\t " << cross5c.cluster.size() << endl;
        //cross5c.write(fname);
        write(fname,cross5c,cross5);
    }
    if (cross3c.cluster.size()>0) { 
        //cout << "write clusters for cross 3': " << sn1+contigName << endl;
        string fname = basename+".cross3.cls.span";
        cout << "write " << fname << "\t " << cross3c.cluster.size() << endl;
        //cross3c.write(fname);
        write(fname,cross3c,cross3);
    }
    
    
}

void C_SpannerCluster::write(string & outfilename, C_NNcluster2d & cls1, 
                             vector<C_localpair> & lpair1) //const
{
    // format: optimized to loadFragments.m matlab script  
    // open output binary file. bomb if unable to open
    fstream output(outfilename.c_str(), ios::out | ios::binary);
    if (!output) {
        cerr << "Unable to open file: " << outfilename << endl;
        return;
    }
    C_headerSpan h;
    h.V = 2203;
    h.setName = setName;
    h.contigName = contigName;
    int t =0;
    while(outfilename.find(h.spanextc[t])==string::npos) t++;
    h.typeName = h.spanextc[t];
    h.reclen =   8*sizeof(double)+sizeof(int);
    h.N = cls1.cluster.size();
    h.write(output);
    int Ncluster = 0;
    vector<int> icluster(lpair1.size(),0);
    C_cluster2d_element1 c1;
    C_cluster2d_elements::iterator it;
    for ( it=cls1.cluster.begin() ; it != cls1.cluster.end(); it++ ) {
        int i = (*it).first;
        c1 = cls1.cluster[i];
        output.write(reinterpret_cast<const char *>(&c1.N), sizeof(int));
        for (int j=0; j<2; j++) {
            output.write(reinterpret_cast<const char *>(&c1.mean[j]), sizeof(double));
            output.write(reinterpret_cast<const char *>(&c1.std[j]), sizeof(double));
            output.write(reinterpret_cast<const char *>(&c1.low[j]), sizeof(double));
            output.write(reinterpret_cast<const char *>(&c1.high[j]), sizeof(double));
        }
        Ncluster++;
        for (int j=0; j<int(c1.inp.size()); j++) {
            icluster[c1.inp[j]]=Ncluster;
        }
    }
    h.typeName = "localpair";
    h.reclen =   4*sizeof(int)+2*sizeof(short)+5*sizeof(char);
    h.N = lpair1.size();
    h.write(output);
    for(int i=0; i<int(lpair1.size()); i++) {
        unsigned int pos = lpair1[i].pos;  
        int lm = lpair1[i].lm;  
        char o = lpair1[i].orient;  
        char q1 = lpair1[i].q1;  
        char q2 = lpair1[i].q2;  
        char mm1 = lpair1[i].mm1;  
        char mm2 = lpair1[i].mm2;  
        short len1 = lpair1[i].len1;  
        short len2 = lpair1[i].len2; 
        unsigned int ReadGroupCode = lpair1[i].ReadGroupCode;  
        int icls = icluster[i]; 
        //position
        output.write(reinterpret_cast<const char *>(&icls), sizeof(int));
        output.write(reinterpret_cast<const char *>(&pos), sizeof(int));
        output.write(reinterpret_cast<const char *>(&lm), sizeof(int));
        output.write(reinterpret_cast<const char *>(&o), sizeof(char));
        output.write(reinterpret_cast<const char *>(&len1), sizeof(short));
        output.write(reinterpret_cast<const char *>(&len2), sizeof(short));
        output.write(reinterpret_cast<const char *>(&q1), sizeof(char));    
        output.write(reinterpret_cast<const char *>(&q2), sizeof(char));    
        output.write(reinterpret_cast<const char *>(&mm1), sizeof(char));    
        output.write(reinterpret_cast<const char *>(&mm2), sizeof(char));    
        output.write(reinterpret_cast<const char *>(&ReadGroupCode), sizeof(int));    
    }
    output.close();
}

//-------------------------------------------------------------
// write clusters for crosspairs
//-------------------------------------------------------------
void C_SpannerCluster::write(string & outfilename, C_NNcluster2d & cls1, 
                             vector<C_crosspair> & xpair1) //const
{
    // format: optimized to loadFragments.m matlab script  
    // open output binary file. bomb if unable to open
    fstream output(outfilename.c_str(), ios::out | ios::binary);
    if (!output) {
        cerr << "Unable to open file: " << outfilename << endl;
        return;
    }
    C_headerSpan h;
    h.V = 2203;
    h.setName = setName;
    h.contigName = contigName;
    int t =0;
    while(outfilename.find(h.spanextc[t])==string::npos) t++;
    h.typeName = h.spanextc[t];
    h.reclen =   8*sizeof(double)+sizeof(int);
    h.N = cls1.cluster.size();
    h.write(output);
    int Ncluster = 0;
    vector<int> icluster(xpair1.size(),0);
    C_cluster2d_element1 c1;
    C_cluster2d_elements::iterator it;
    for ( it=cls1.cluster.begin() ; it != cls1.cluster.end(); it++ ) {
        int i = (*it).first;
        c1 = cls1.cluster[i];
        output.write(reinterpret_cast<const char *>(&c1.N), sizeof(int));
        for (int j=0; j<2; j++) {
            output.write(reinterpret_cast<const char *>(&c1.mean[j]), sizeof(double));
            output.write(reinterpret_cast<const char *>(&c1.std[j]), sizeof(double));
            output.write(reinterpret_cast<const char *>(&c1.low[j]), sizeof(double));
            output.write(reinterpret_cast<const char *>(&c1.high[j]), sizeof(double));
        }
        Ncluster++;
        for (int j=0; j<int(c1.inp.size()); j++) {
            icluster[c1.inp[j]]=Ncluster;
        }
    }
    
    h.typeName = "crosspair";
    h.reclen =    2*(sizeof(int)+2*sizeof(short)+3*sizeof(char))+sizeof(int);
    h.N = xpair1.size();
    h.write(output);
    
    for(int i=0; i<int(xpair1.size()); i++) {
        for (int e=0; e<2; e++) {
            unsigned int p0 = xpair1[i].read[e].pos;
            unsigned short l0 = xpair1[i].read[e].len;          // length of this read aligment 
            unsigned short a0 = xpair1[i].read[e].anchor;    // anchor index  
            char  s0 = xpair1[i].read[e].sense;    // anchor index  
            char  q0 = xpair1[i].read[e].q;    // anchor index  
            char  mm0 = xpair1[i].read[e].mm;    // anchor index  
            //position
            output.write(reinterpret_cast<const char *>(&p0), sizeof(int));
            //length
            output.write(reinterpret_cast<const char *>(&l0), sizeof(short));
            output.write(reinterpret_cast<const char *>(&a0), sizeof(short));
            output.write(reinterpret_cast<const char *>(&s0), sizeof(char));
            output.write(reinterpret_cast<const char *>(&q0), sizeof(char));
            output.write(reinterpret_cast<const char *>(&mm0), sizeof(char));
        }
        unsigned int ReadGroupCode0 = xpair1[i].ReadGroupCode;
        output.write(reinterpret_cast<const char *>(&ReadGroupCode0), sizeof(int));    
    }
    output.close();
}



//------------------------------------------------------------------------------
//  Clustering class ( for Retro mob insertions) 
//------------------------------------------------------------------------------ 
C_SpannerRetroCluster::C_SpannerRetroCluster(C_contig & c1,  C_libraries & libs1, RunControlParameters & pars1, 
                                             int e, string & retrotype) {
    
    //library info
    libraries = libs1;
    //  parameters
    pars  = pars1;
    // retro element bit index
    etype = e;
    // retro element type name
    typeName=retrotype;
    contigName=c1.getContigName();
    setName=c1.setName; 
    L = c1.Length;   
    double Npair = c1.localpairs.size();
    PairDensity = Npair/L;
    double Nretro = c1.umpairs.size();
    RetroDensity = Nretro/L;
    //-------------------------------------------------------------------------
    // select retro free element frags  -> alu5, alu3
    // UM with M hitting alu,  and no satisfied frag length criteria
    //-------------------------------------------------------------------------
    int Qmin = 0; //pars.getQmin();
    selectRetro(c1, Qmin);
    //-------------------------------------------------------------------------
    // set neighborhood window to some scale of fragment length width 
    // for breakpoint spanning pairs
    //-------------------------------------------------------------------------
    // HistObj h = pars.getFragHist();
    // double window = 1.0*h.median;    // 24 April 2009
    
    double fwindow = 1.0;    // 12 July 2009 - relative to lib frag width
    pars.setClusteringLength(int(fwindow));
    //-------------------------------------------------------------------------
    // long /short LM clusters for inversion detection
    //-------------------------------------------------------------------------
    double fwindow2[2] = {fwindow*1.25,1e-3};
    //-------------------------------------------------------------------------
    // both forward mapped ~ '>' reads 5'
    //-------------------------------------------------------------------------
    vector<vector<double> > x2;
    x2.clear();    
    vector<vector<double> > wx2;
    wx2.clear();    
    vector<int> ip2;
    ip2.clear();    
    
    int N =  makepairP(e5,'F',x2,wx2,ip2);    
    /*
     double window = pars.getFragmentLengthHi();     // 27 April 2009
     pars.setClusteringLength(int(window));
     //-------------------------------------------------------------------------
     // cluster r5's  (F UM fragments with element)
     //-------------------------------------------------------------------------
     vector<double>  x1;
     x1.clear();    
     int N =  makeRetroX(e5,'F',x1);    
     */
    if (N>0) {
        C_NNcluster2d q(x2, wx2, ip2, fwindow2,"e5pair",c1.setName,c1.getContigName());  
        //q.mergeClusters();   
        q.cleanClusters();     
        
        e5c = q; 
    }
    // cout << e5c.typeName << "\t " << e5c.NC << endl;
    //-------------------------------------------------------------------------
    // cluster r3's (F UM fragments with element)
    //-------------------------------------------------------------------------
    x2.clear();    
    wx2.clear();    
    ip2.clear();    
    N =  makepairP(e3,'R',x2,wx2,ip2);    
    //N =  makeRetroX(e3,'R',x1);    
    if (N>0) {
        C_NNcluster2d q(x2, wx2, ip2, fwindow2,"e3pair",c1.setName,c1.getContigName());  
        //C_NNcluster1d q(x1, window,"e5pair",c1.setName,c1.getContigName());  
        q.cleanClusters();     
        e3c = q; 
    }
    // cout << e3c.typeName << "\t " << e3c.NC << endl;
}

//------------------------------------------------------------------------------
// select retro fragments for clustering...
//------------------------------------------------------------------------------
void C_SpannerRetroCluster::selectRetro(C_contig & c1,int Qmin) {
    // init
    int N0 =  c1.umpairs.size();    
    e5.clear();        
    e3.clear();          
    
    /*
     int aly elementmask  = 1;
     int alu elementmask  = 2;
     int l1  elementmask  = 4;
     int sva elementmask  = 8;
     */
    int elementmask = 1<<etype;
    
    // frag length parameters
    // double lmLow = pars.getFragmentLengthLo();
    // double lmHigh = pars.getFragmentLengthHi();    
    //double fraglen = pars.getFragmentLength();        
    
    list<C_umpair>::iterator i;
    if (N0>0) {
        for(i=c1.umpairs.begin(); i != c1.umpairs.end(); ++i) {
            
            //----------------------------------------------------------------------
            // only high quality mapped read pairs
            //----------------------------------------------------------------------        
            // kludge for subzero q getting rolled over to 100+ something...  
            if ((*i).read[0].q>100) (*i).read[0].q=0;
            if (int((*i).read[0].q)<Qmin) continue;
            
            // U side
            char su=(*i).read[0].sense;                   // forward ('F') or reverse complement ('R')
            unsigned int pu=(*i).read[0].pos;             // position in contig (unpadded)
            unsigned short lu=(*i).read[0].len;           // read length
            unsigned short au=(*i).read[0].anchor;        // anchor index  
            //char qu=(*i).read[0].q;                     // mapping quality
            // M side
            // char sm=(*i).read[1].sense;                // forward ('F') or reverse complement ('R')
            unsigned int pm=(*i).read[1].pos;             // position in contig (unpadded)
            unsigned short lm=(*i).read[1].len;           // read length
            unsigned short am=(*i).read[1].anchor;        // anchor index  
            //char qm=(*i).read[1].q;                     // mapping quality
            // element
            int element = (*i).elements;
            // M mappings
            //int nmap = (*i).nmap;
            
            unsigned int RGC1=(*i).ReadGroupCode;
            double lm1=libraries.libmap[RGC1].LM;
            
            // skip flipped libs (i.e. unfixed Jumping libs)
            if (lm1<0) continue;

            
            double lmHigh=libraries.libmap[RGC1].LMhigh;
            double lmLow=libraries.libmap[RGC1].LMlow;
            
            //----------------------------------------------------------------------
            // widen resolve window by 2x to clean up clusters
            //----------------------------------------------------------------------
            lmHigh=lmHigh+(lmHigh-lm1);
            lmLow=lmLow-(lm1-lmLow);
            
            // frag constraint
            bool constrain=false;
            if (au==am) {
                double mapumdist = double(pm)+lm-double(pu);
                if (su=='R') {
                    mapumdist = double(pu)+lu-double(pm);
                } 
                constrain = (mapumdist>lmLow)&(mapumdist<lmHigh);
            }      
            
            // element bit set ?
            bool isae = (element & elementmask) > 0;
            
            // sva should be only bit 
            // ALU can be any alu or aluy bit 
            if (etype==3) { 
                isae = (element == elementmask); // sva alone
            }
            
            // if  (isae)  {
            // include only um frags with element and not resolved 
            if ( (isae) & (!constrain) ) {
                if (su=='F') {
                    e5.push_back((*i));
                } else {
                    e3.push_back((*i));
                }
            }        
        }
    }
}

//------------------------------------------------------------------------------
//  umpair to simple vector function 
//------------------------------------------------------------------------------
int C_SpannerRetroCluster::makepairP(vector<C_umpair> & p1,  char orient
                                     , vector<vector<double> >& x1, vector<vector<double> >& wx1, vector<int> & inp1) {
    
    int N0 =  p1.size();    
    x1.clear();          
    wx1.clear();          
    inp1.clear();
    
    C_vectorDouble x4(5,0);
    vector<double> x2(2,0);
    
    // use a list (for sorting...)  
    list<C_vectorDouble> xlist; 
    list<C_vectorDouble>::iterator it;
    
    if (N0>0) {
        for(int i=0; i< N0; i++) {
            
            if (p1[i].read[1].pos>6000) {
                //--------------------------------------------------------------------
                // skip poly-A tail in L1.HS 
                //--------------------------------------------------------------------
                continue;
            }
            char o1=p1[i].read[0].sense;  
            if (o1!=orient) {
                cerr << " wrong orientation in makepairP" << p1[i] << endl;
            }
            // index to p1 order
            x4[4]=i;
            
            // position at center of fragment
            // x2[0]=double(p1[i].pos);  5' end of fragment
            // p1[i].read[0].pos;
            // get average LF from library info record
            int LF1 = libraries.libmap[p1[i].ReadGroupCode].LM;
            
            // skip flipped libs (i.e. unfixed Jumping libs)
            if (LF1<0) continue;

            // get library fragment width info
            int LFlow = libraries.libmap[p1[i].ReadGroupCode].LMlow;
            int LFhigh = libraries.libmap[p1[i].ReadGroupCode].LMhigh;
            double LFW=double(LFhigh-LFlow);
            
            //----------------------------------------------------------------------
            // x4[0] is position to cluster
            // position to middle of frag ~ edge of element
            //----------------------------------------------------------------------
            if (orient=='F') { 
                x4[0]=double(p1[i].read[0].pos)+double(LF1)/2.0;  
            } else {
                x4[0]=double(p1[i].read[0].pos)+double(p1[i].read[0].len)-double(LF1)/2.0; 
            }
            // x4[1] is a dummy var for the 2d clustering
            x4[1]=0;  
            // neighborhood for this frag
            x4[2]=LFW;
            // dummy small neighborhood
            x4[3]=1e-3;
            //wx1.push_back(x);
            xlist.push_back(x4);
        }
    }
    
    //sort
    xlist.sort(compare_vector);
    
    // dump into x1 and wx1
    for (it=xlist.begin(); it!=xlist.end(); ++it) {
        x4=*it;
        x2[0]=x4[0];
        x2[1]=x4[1];
        x1.push_back(x2);
        x2[0]=x4[2];
        x2[1]=x4[3];
        wx1.push_back(x2);
        inp1.push_back(int(round(x4[4])));
        
    }
    return int(x1.size());
}

//------------------------------------------------------------------------------
//  obsolete 8/2009 umpair to simple vector function 
//------------------------------------------------------------------------------
int C_SpannerRetroCluster::makeRetroX(vector<C_umpair> & r1,  char orient, vector<double> & x1) {
    int N0 =  r1.size();    
    x1.clear();          
    double x0;
    if (N0>0) {
        for(int i=0; i< int(r1.size()); i++) {
            char o1=r1[i].read[0].sense;  
            if (o1!=orient) {
                cerr << " wrong orientation in makepairP" << r1[i] << endl;
            }
            x0=double(r1[i].read[0].pos);  
            x1.push_back(x0);
        }
        
    }
    return int(x1.size());
}

//------------------------------------------------------------------------------
// select retro fragments for clustering...
//------------------------------------------------------------------------------
void C_SpannerRetroCluster::writeall() {
    string area = pars.getOutputDir();
    int k = 1+area.find_last_of("/");
    string subdir = area.substr(k);
    if (area.size()>0) { area = area+"/";}
    string prefix = pars.getPrefix();
    if (prefix.size()>0) { prefix = prefix+".";}
    string sn1 = setName;
    // dont need redundant setName if subdirectory is already the setName
    // if (sn1==subdir) sn1="";
    if (sn1.size()>0) { sn1 = sn1+".";}
    string basename = area+prefix+sn1+contigName;
    if (setName==contigName) basename = area+prefix+contigName;
    
    // replace evil character "|" with benign "_" 
    size_t found=basename.find("|");
    while (found!=string::npos) {
        basename.replace(found,1,"_");
        found=basename.find("|");
    }
    
    if (e5c.cluster.size()>0) { 
        // cout << "write clusters for 5': " << typeName << " " << sn1+contigName << endl;
        string fname = basename+"."+typeName+"5.cls.span";
        cout << "write " << fname << "\t " << e5c.cluster.size() << endl;
        write(fname,e5c,e5);
    }
    if (e3c.cluster.size()>0) { 
        //cout << "write clusters for 3': " << typeName << " " << sn1+contigName << endl;
        string fname = basename+"."+typeName+"3.cls.span";
        cout << "write" << fname << "\t " << e3c.cluster.size() << endl;
        write(fname,e3c,e3);
    }
}

void C_SpannerRetroCluster::write(string & outfilename, C_NNcluster2d & cls1, 
                                  vector<C_umpair> & retro1) //const
{
    // format: optimized to loadFragments.m matlab script  
    // open output binary file. bomb if unable to open
    fstream output(outfilename.c_str(), ios::out | ios::binary);
    if (!output) {
        cerr << "Unable to open file: " << outfilename << endl;
        return;
    }
    C_headerSpan h;
    h.V = 2203;
    h.setName = setName;
    h.contigName = contigName;
    /*
     int t =0;  
     int  nt=h.spanextc.size();
     while( (outfilename.find(h.spanextc[t])==string::npos)&&(t<nt)) t++;
     h.typeName = h.spanextc[t];
     */
    h.typeName = typeName;
    h.reclen =   4*sizeof(double)+sizeof(int);
    h.N = cls1.cluster.size();
    h.write(output);
    int Ncluster = 0;
    vector<int> icluster(retro1.size(),0);
    C_cluster2d_element1 c1;
    C_cluster2d_elements::iterator it;
    for ( it=cls1.cluster.begin() ; it != cls1.cluster.end(); it++ ) {
        int i = (*it).first;
        c1 = cls1.cluster[i];
        output.write(reinterpret_cast<const char *>(&c1.N), sizeof(int));
        output.write(reinterpret_cast<const char *>(&c1.mean[0]), sizeof(double));
        output.write(reinterpret_cast<const char *>(&c1.std[0]), sizeof(double));
        output.write(reinterpret_cast<const char *>(&c1.low[0]), sizeof(double));
        output.write(reinterpret_cast<const char *>(&c1.high[0]), sizeof(double));
        Ncluster++;
        for (int j=0; j<int(c1.inp.size()); j++) {
            icluster[c1.inp[j]]=Ncluster;
        }
    }
    h.typeName = "umpair";
    
    // two reads + nmap
    h.reclen =6*sizeof(int)+4*sizeof(short)+6*sizeof(char);
    h.N = retro1.size();
    h.write(output);
    // loop
    //vector<C_umpair>::iterator i;
    for(int i=0; i<int(retro1.size()); i++) {
        //for(i=retro1.begin(); i != retro1.end(); ++i) {
        // cluster
        //unsigned int ReadGroupCode = (*i).ReadGroupCode;  
        unsigned int ReadGroupCode = retro1[i].ReadGroupCode;  
        int icls = icluster[i]; 
        output.write(reinterpret_cast<const char *>(&icls), sizeof(int));
        output.write(reinterpret_cast<const char *>(&ReadGroupCode), sizeof(int)); 
        for (int e=0; e<2; e++) {
            unsigned int p0 = retro1[i].read[e].pos;
            unsigned short l0 = retro1[i].read[e].len;       // length of this read aligment 
            unsigned short a0 = retro1[i].read[e].anchor;    // anchor index  
            char  s0 = retro1[i].read[e].sense;              // read orientation 
            char  q0 = retro1[i].read[e].q;                  // mapping q
            char  mm0 = retro1[i].read[e].mm;                // mismatches
            output.write(reinterpret_cast<const char *>(&p0), sizeof(int));
            output.write(reinterpret_cast<const char *>(&l0), sizeof(short));
            output.write(reinterpret_cast<const char *>(&a0), sizeof(short));
            output.write(reinterpret_cast<const char *>(&s0), sizeof(char));
            output.write(reinterpret_cast<const char *>(&q0), sizeof(char));
            output.write(reinterpret_cast<const char *>(&mm0), sizeof(char));    
            
        }
        int nmap= retro1[i].nmap;   //nmapped positions
        output.write(reinterpret_cast<const char *>(&nmap), sizeof(int));
        int elements= retro1[i].elements;   //nmapped positions
        output.write(reinterpret_cast<const char *>(&elements), sizeof(int));
    }  
    
    output.close();
}


int C_SpannerRetroCluster::Mask(C_BedChr & mask, int LMmax) {
    
    // set start of mask loop 
    int im0=0;
    // 
    C_cluster2d_element1 c1;
    C_cluster2d_elements::iterator it;
    vector<int> itoss;
    
    
    // 5' clusters
    for ( it=e5c.cluster.begin() ; it != e5c.cluster.end(); it++ ) {
        int i = (*it).first;
        c1 = e5c.cluster[i];
        //int pc1 = int(e5c.cluster[i].high[0]);
        int pc1 = int(e5c.cluster[i].mean[0]);
        int pc2 = pc1+LMmax;
        int im; 
        bool hit=false;
        for (im=im0; im<int(mask.p.size()); im++) {
            int pm1=mask.p[im];
            int pm2=mask.p[im]+int(mask.l[im]);        
            hit= ( (pc1<=pm1) & (pm1<=pc2) );
            hit= hit | ( (pm1<=pc1) & (pc1<=pm2) );
            if (hit) { 
                itoss.push_back(i);
                break;
            }
            if (pm1>pc2) { break;}
        }
        if (im>im0) {
            im0=im; 
        }
    }
    
    // remove masked 5' clusters
    int NT = itoss.size();
    for (int i=0; i<NT; i++) {
        e5c.cluster.erase(itoss[i]);
    }
    int NC5 = e5c.cluster.size();
    
    
    // 3' clusters
    im0=0;
    itoss.clear();
    for ( it=e3c.cluster.begin() ; it != e3c.cluster.end(); it++ ) {
        int i = (*it).first;
        c1 = e3c.cluster[i];
        //int pc2 = int(e3c.cluster[i].low[0]);
        int pc2 = int(e3c.cluster[i].mean[0]);
        int pc1 = pc2-LMmax;
        int im;
        bool hit=false;
        for (im=im0; im<int(mask.p.size()); im++) {
            int pm1=mask.p[im];
            int pm2=mask.p[im]+int(mask.l[im]);        
            hit= ( (pc1<=pm1) & (pm1<=pc2) );
            hit= hit | ( (pm1<=pc1) & (pc1<=pm2) );
            if (hit) { 
                itoss.push_back(i);
                break;
            }
            if (pm1>pc2) { break;}
        }
        if (im>im0) {
            im0=im; 
        }
    }
    NT = itoss.size();
    for (int i=0; i<NT; i++) {
        e3c.cluster.erase(itoss[i]);
    }
    int NC3 = e5c.cluster.size();
    
    return (NC3+NC5);
    
}


//------------------------------------------------------------------------------
//  Coverage class
//------------------------------------------------------------------------------ 
C_SVcoverage1::C_SVcoverage1() {
    p0=0;               // start base position of region
    p1=0;               // end base position
    Nsite=0;            // number of available start sites in region
    Nrepeat=0;          // number of masked start sites in region
    N=0;                // number of reads 
    eN=0;            // expected number of reads  
    p=0;   
}

C_SVcoverage1::C_SVcoverage1(C_contig & contig, int P0, int P1,C_NominalCov & nc) {
    p0 = P0;               // start base position of region
    p1 = P1;               // end base position
    Nrepeat=0;             // number of masked start sites in region
    
    N = 0;    
    size_t Lrn= contig.repeat.n.size(); 
    
    for (int i=p0; i<=p1; i++) {
        N +=int(contig.read_start.n[i]); 
        if (i<int(Lrn)) Nrepeat+=int(contig.repeat.n[i]>0); 
    }
    Nsite=p1-p0+1-Nrepeat;    // number of available start sites in region
    //int totRepeats = int(contig.repeat.Stats.N*contig.repeat.Stats.mean);
    int totSites = contig.Length - int(contig.totalRepeatBases)-int(contig.totalNoCovBases);
    eN = float(Nsite*contig.totalUniqueReads/totSites);  // expected number of reads  
    float r1 = float(N)/eN;
    p = nc.getplow(Nsite,r1);
}

ostream &operator<<(ostream &output, C_SVcoverage1 & q1)
{
    char b [100];
    sprintf(b,"%10d %10d %10d %10d %10d %10.1f %10.5f", q1.N,q1.p0,q1.p1,q1.Nsite,q1.Nrepeat,q1.eN,q1.p);
    string s = b;
    output << s;
    return output;
}

C_SVcoverage1& C_SVcoverage1::operator=(const C_SVcoverage1 &rhs)
{
    this->p0 = rhs.p0;
    this->p1 = rhs.p1;
    this->Nsite = rhs.Nsite;
    this->Nrepeat = rhs.Nrepeat;
    this->N= rhs.N;
    this->eN = rhs.eN;
    this->p = rhs.p;
    return *this;
}

int C_SVcoverage1::operator==(const C_SVcoverage1 &rhs) const
{
    if( this->p0 != rhs.p0) return 0;
    if( this->p1 != rhs.p1) return 0;
    if( this->Nsite != rhs.Nsite) return 0;
    if( this->Nrepeat != rhs.Nrepeat) return 0;
    if( this->N != rhs.N) return 0;
    if( this->eN != rhs.eN) return 0;
    if( this->p != rhs.p) return 0;
    return 1;
}

int C_SVcoverage1::operator<(const C_SVcoverage1 &rhs) const
{
    if( this->p0 < rhs.p0 ) return 1;
    if( this->p0 == rhs.p0 && this->p1 < rhs.p1 ) return 1;
    return 0;
}



C_NominalCov::C_NominalCov()  {
    N=0;
}
//------------------------------------------------------------------------------
// create class for Nominal coverage histograms given set (data) and 
// nuv (vector of region lengths) 
//------------------------------------------------------------------------------
C_NominalCov::C_NominalCov(C_set & set, vector<int> & nuv)  {
    nu = nuv;
    N = nu.size();
    // loop over contigs to add up fractions 
    vector<unsigned int> Lt;
    vector<int> L;
    vector<string> cn;
    C_contigs::const_iterator iterContig;   
    unsigned int Ltot=0; 
    for (iterContig = set.contig.begin();	iterContig != set.contig.end(); iterContig++) {
        string name = iterContig->first;
        if (set.contig[name].Length>0) {
            int L1 =  set.contig[name].Length;
            Ltot+=L1;
            L.push_back(L1);
            Lt.push_back(Ltot);
            cn.push_back(name);
        }
    }
    int Nc=L.size();
    /* initialize random seed: */
    srand ( 347714 );
    //----------------------------------------------------------------------------
    // number of samples & number of bins parameters
    //----------------------------------------------------------------------------
    int imax=1000;
    int nbin=500;
    //
    for (int n=0; n<N; n++) {
        int nu1 = nu[n];
        HistObj h;
        //h.Initialize(nu1+1,-0.5,nu1+0.5);  
        h.Initialize(nbin,0,5.0);  
        char buffer [100];
        sprintf(buffer," Nu %d:  reads / expected",nu1); 
        string title = buffer;
        h.setTitle(title);
        h.setXlabel("Nr/<Nr>");
        for (int i=0; i<imax; i++) {
            int r1  = rand() % Ltot;
            int c = 0;
            for (c=0; c<Nc; c++) {
                if (r1>int(Lt[c])) {break;}
            }
            c--;
            int dp = L[c];
            int p0=0;
            int p1=0;
            //-----------------------------------------------------------------------
            // don't sample coverage from repeat wastelands with < 10% unique left...
            //-----------------------------------------------------------------------
            while (dp>(10*nu1)) {
                p0  = rand() % (L[c]-2*nu1);
                p1  = set.contig[cn[c]].howFar(p0,nu1);
                dp = p1-p0;
            }
            bool wacky=((p0<0)|(p0>=int(L[c]))|(p1<0)|(p1>=int(L[c])));
            if (wacky) {
                cout << " wacky limits in Nominal_Cov " << p0 << "\t to " << p1 << endl;
                cout << " wacky ... L[c] " << L[c] << "\t nu1 " << nu1 << endl;
                exit(0);
            }      
            //------------------------------------------------------------------------
            // this had problems with unitialized hist for this 
            //C_SVcoverage1 sc1(set.contig[cn[c]], p0, p1, *this);
            // do the count and Ecount by hand here....
            //------------------------------------------------------------------------
            double count = 0;
            double repeat =0;
            size_t Lrn= set.contig[cn[c]].repeat.n.size(); 
            for (int i=p0; i<=p1; i++) {
                count +=int(set.contig[cn[c]].read_start.n[i]); 
                if (i>=int(Lrn)) continue; 
                repeat+=int(set.contig[cn[c]].repeat.n[i]>0); 
            }      
            double Nsite=p1-p0+1-repeat;    // number of available start sites in region
            //int totRepeats = int(contig.repeat.Stats.N*contig.repeat.Stats.mean);
            int totSites = set.contig[cn[c]].Length - int(set.contig[cn[c]].totalRepeatBases)-int(set.contig[cn[c]].totalNoCovBases);
            float Ecount = float(Nsite*set.contig[cn[c]].totalUniqueReads/totSites);  // expected number of reads  
            
            //int Nr = sc1.N;
            //h.Fill1(Nr);  
            //h.Fill1(sc1.N/sc1.eN);  
            h.Fill1(count/Ecount);  
        }
        h.Finalize();
        hist[nu1]=h;
        // cout << hist[nu1] << endl;
    }
}

double C_NominalCov::getplow(int nu1, double r1) {
    int n;
    double p;
    int nu0=nu[0];
    if (hist[nu0].Nbin==0) {
        return 0;
    }
    if (nu1<=nu0) {
        p=hist[nu0].x2pTrim(r1);
        return p;
    } else if (nu1>=nu[N-1]) {
        p=hist[nu[N-1]].x2pTrim(r1);
        return p;
    }
    for (n=1; n<N; n++) {
        if ( (nu1<nu[n]) & ( nu1>=nu[n-1]) ) {break;}
    }
    if (n==N) {
        cout << " bug in getplow " << nu1 << endl;
        exit(0);
    }
    int n0= nu[n-1];
    double p0=hist[n0].x2pTrim(r1); 
    int n1= nu[n];
    double p1=hist[n1].x2pTrim(r1);
    p = interp(n0,p0,n1,p1,nu1);
    return p;
}

double C_NominalCov::interp(int i0,double x0,int i1,double x1,int i) {
    double x = x0+double(i-i0)*(x1-x0)/double(i1-i0);
    return x;
}


ostream &operator<<(ostream &output, C_NominalCov & sc)
{
    output << "Nominal coverage samples:\t " << sc.N << endl;
    output << "\t ";
    char b[100];
    sprintf(b,"%10s %10s %10s", "Nub","<N/eN>","stdev");
    string s = b;
    output << s << endl;
    for (int i=0; i<sc.N; i++) {
        int nu1 = sc.nu[i];
        sprintf(b,"%10d %10.3f %10.3f", nu1,sc.hist[nu1].mean,sc.hist[nu1].std);
        s = b;
        output << s << endl;
    }
    for (int i=0; i<sc.N; i++) {
        int nu1 = sc.nu[i];
        output << sc.hist[nu1];
    }
    return output;
}

//------------------------------------------------------------------------------
// create list of retro elements from the anchor info
//------------------------------------------------------------------------------
C_retroElements::C_retroElements(C_anchorinfo & a1, string & prefix) {
    
	// element name pattern
	string patternElement("^"+prefix+"(\\S+)\\.\\S+");
	string match;
	
	
    
    // fill elements in anchor object
    a1.anchorElements();
    int n= a1.element.size();
    e.clear();
    name.clear();
    N=0;
    std::list<int>  elist;
    std::list<int>::iterator it;
    std::map<int, string, std::less<int> >  nmap;
    for (int i = 0; i<n; i++) {
        if (a1.element[i]>=0) {
            bool gotit=false;
            for (it=elist.begin() ; it != elist.end(); it++ ) {
                int e1=*it;
                gotit=(e1==a1.element[i]);
                if (gotit) {break;}
            }
            if (!gotit) {
                N++;
                elist.push_back(a1.element[i]);
                nmap[a1.element[i]]=a1.names[i].substr(0,3);
				
                // check pattern for template
                if (RE2::FullMatch(a1.names[i].c_str(),patternElement.c_str(),&match) ) {
				  	nmap[a1.element[i]] = match;
				}
				
            }
        }
    }
    elist.sort();
    for (it=elist.begin() ; it != elist.end(); it++ ) {
        int e1=*it;
        e.push_back(e1);
        name.push_back(nmap[e1]);
    }
}
