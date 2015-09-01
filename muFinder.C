/*
	muFinder.C
	Clint Wiseman, USC/Majorana
	August 2015.
	
 => This program can be run on PDSF or locally.  It takes a .txt file of run numbers
	as an input argument, and uses the name of the text file to generate output
	in the folder ./output
 
 => More info goes here .........

 => May want to put the "processing functions" into a separate file, and include a header,
 	like this: http://stackoverflow.com/questions/25274312/defining-c-functions-inside-header-files

	Usage:
	CINT: root[0] .X muFinder.C  <-- uses "Debug" file list
		  root[0] .X muFinder.C ("Filename_without_extension")  <--- NO .TXT extension.
	
	Compiled mode: root[0] .X muFinder.C++ ("Filename_without_extension")
	
	bash: root -b -q -l muFinder.C ("Filename_without_extension")
		  rootMJ -b  muFinder.C++ \(\"Debug1\"\)

*/

// Resolving problem with OS X stdint.h header
// stdint.h on >10.9 is "too complex" for CINT to parse, needs CLANG.
// Workaround from https://github.com/psi46/pxar/issues/20
#if defined(__CINT__) && defined (__APPLE__)
#undef GNUC 
typedef char __signed; 
typedef char int8_t; 
#endif

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>

// MUST load MDGOMJ classes in ROOT before calling these headers
#include "MJTRun.hh"
#include "MJTVetoData.hh"
#include "MGTBasicEvent.hh"

using namespace std;

const int numPanels = 32;
const int vd_size = 16;	
const Int_t numPlanes = 7;

// Structure to hold a complete matched veto event 
struct VEvent {
	Int_t isGood;		// flag a bad VEvent.
    Int_t run;			// run number
    Int_t vEnt; 		// qdc eventCount (i.e. entry)
    Int_t sEnt; 		// scaler entry
    Int_t card[3];		// 0: scaler ID, 1: qdc card 1 (panels 1-16), 2: qdc card 2 (panels 17-32)
    Float_t xTime;		// "combination time"
    Float_t sTime;		// scaler time
    Bool_t sTimeBad;	// bad scaler entry
    Float_t lTime;		// LED time
    Float_t eTime;		// entry time
    Int_t QDCTotal;		// summed QDC from all panels
    Int_t QDC[numPanels];	// store qdc values for both cards in the same array 
    Int_t IsUnderThreshold[numPanels];
    Int_t IsOverflow[numPanels];
	Bool_t coinType[32]; 	// coincidence type
	Bool_t cutType[32];		// bools used in cuts 
	Bool_t plane[numPlanes];	// which planes were hit: 0:L-Bot 1:U-Bot 2:Top 3:North 4:South 5:West 6:East
	Int_t planeHitCount[numPlanes];  // how many hits in planes

    // constructor
    VEvent(){
    	isGood=1; run=0; vEnt=0; sEnt=0;
    	card[0]=0;card[1]=0;card[2]=0;
    	xTime=0; sTime=0; sTimeBad=0; lTime=0; eTime=0;
    	QDCTotal = 0;
    	for(int i=0;i<numPanels;i++) {
    		QDC[i]=0; IsUnderThreshold[i]=0; IsOverflow[i]=0;
    	}
    	for(int j=0;j<10;j++) { coinType[j]=0; cutType[j]=0; }
    	for(int k=0;k<numPlanes;k++) { plane[k]=0; planeHitCount[k]=0; }
    }
};

// Structure to hold a qdc entry
struct QEvent{
    Int_t runNumber;
    UInt_t crate;
    UInt_t card;
    UInt_t EventCount;
	Int_t QDC[16];
    Int_t IsUnderThreshold[vd_size];
    Int_t IsOverflow[vd_size];

    // constructor
    QEvent(){
    	runNumber=0; crate=0; card=0; EventCount=0;
    	for(int i=0;i<vd_size;i++){
    		QDC[i]=0; IsUnderThreshold[i]=0; IsOverflow[i]=0;	
    	}
    }
};

// Structure to hold info about the veto panel's LED peak
// Read these in from a calibration table made by builtVetoCal.
struct PeakInfo {
	Int_t panel;
	Float_t mean;
	Float_t meanErr;
	Float_t sigma;
	Float_t sigmaErr;
	Float_t chi2_NDF;
	Bool_t badFit;

	// constructor
	PeakInfo() {
		panel=0; mean=0; meanErr=0;
		sigma=0; sigmaErr=0; chi2_NDF=0;
		badFit=false;
	}
};

//==============PROCESSING FUNCTIONS=========================

// move data into QEvent structure
QEvent WriteQEvent(Int_t verbose, Int_t run, Int_t i, MJTVetoData *vd[vd_size]) {
	QEvent qdc;
	qdc.runNumber=run;
	qdc.crate=vd[0]->GetCrate();
	qdc.card=vd[0]->GetCard();
	qdc.EventCount=vd[0]->GetEventCount();
	int k = 0;
	for (int j = 0; j<vd_size; j++)	{
	k = vd[j]->GetChannel();
		qdc.QDC[k]=vd[j]->GetAmplitude();
		qdc.IsUnderThreshold[k]=vd[j]->IsUnderThreshold();
		qdc.IsOverflow[k]=vd[j]->IsOverflow();
	}

	if (verbose==2) {
	// check qdc data after moving into QEvent structure
		printf("QDC -- run: %i  Entry: %i  crate: %i  card: %i  EventCount: %i \n",qdc.runNumber,i,qdc.crate,qdc.card,qdc.EventCount);
		cout << "QDC: "; for (int k = 0; k<16; ++k) { cout << qdc.QDC[k] << "  "; } cout << endl;
		cout << "UTh: "; for (int k = 0; k<16; ++k) { cout << qdc.IsUnderThreshold[k] << "  "; } cout << endl;
		cout << "Ovr: "; for (int k = 0; k<16; ++k) { cout << qdc.IsOverflow[k] << "  "; } cout << endl;	
	}
	if (verbose==1) {
	// check entry numbers
		printf("Entry: %i  card: %i  EventCount: %i",i,qdc.card,qdc.EventCount);
		printf("  ||  ScalerCount: %i  TimeStamp: %.5f  IsBadTs: %i \n",
			vd[0]->GetScalerCount(),vd[0]->GetTimeStamp()/1E8,vd[0]->IsBadTS());
	}

	return qdc;
}

// set flag if current qdc entry has same EventCount as previous one.
Int_t EventMatch(Int_t i, Int_t run, Long64_t nentries, QEvent qdc, QEvent prevqdc) {
	if (qdc.EventCount == prevqdc.EventCount) { 
		//printf("EventMatch true.  qdc.EventCount:%i  prevqdc.EventCount:%i \n",qdc.EventCount,prevqdc.EventCount);
		return 1;
	}
	else if (abs((Int_t)qdc.EventCount - (Int_t)prevqdc.EventCount) > 1 && i > 2) {
		printf(" EventCount mismatch! Run:%i  Entry:%i  Card:%i  Prev.card:%i  Breaking at %.0f %% through file.\n",
			run,i,qdc.card,prevqdc.card,((Float_t)i/nentries)*100);
		return 0;
	}
	else return 0;
}

// move matching events into VEvent structure and incorporate scaler info.
VEvent WriteVEvent(Int_t verbose, Int_t i, Int_t run, Long64_t nentries, Float_t duration, Int_t card1, 
	Int_t card2, QEvent qdc, QEvent prevqdc, MJTVetoData *vd[vd_size]) {

	VEvent veto;
	veto.run = qdc.runNumber;
	veto.vEnt = qdc.EventCount;
	veto.sEnt = vd[0]->GetScalerCount();
	veto.sTime = vd[0]->GetTimeStamp()/1E8;
	veto.sTimeBad = vd[0]->IsBadTS();
	veto.eTime = ((Float_t)i/nentries)*duration;	
	// case 1
	if ((Int_t)prevqdc.card==card1 && (Int_t)qdc.card==card2) {
		veto.card[0]=vd[0]->GetScalerID();
		veto.card[1]=prevqdc.card;
		veto.card[2]=qdc.card;
		for (int k = 0; k<16; k++) {
			veto.QDC[k]=prevqdc.QDC[k];
			veto.QDC[16+k]=qdc.QDC[k];
			veto.QDCTotal += veto.QDC[k];
			veto.QDCTotal += veto.QDC[16+k];
			veto.IsUnderThreshold[k]=prevqdc.IsUnderThreshold[k];
			veto.IsUnderThreshold[16+k]=qdc.IsUnderThreshold[k];
			veto.IsOverflow[k]=prevqdc.IsOverflow[k];
			veto.IsOverflow[16+k]=qdc.IsOverflow[k];
		}
	}
	// case 2
	else if ((Int_t)prevqdc.card==card2 && (Int_t)qdc.card==card1) {
		veto.card[0]=vd[0]->GetScalerID();
		veto.card[1]=qdc.card;
		veto.card[2]=prevqdc.card;
		for (int k = 0; k<16; k++) {
			veto.QDC[k]=qdc.QDC[k];
			veto.QDC[16+k]=prevqdc.QDC[k];
			veto.QDCTotal += veto.QDC[k];
			veto.QDCTotal += veto.QDC[16+k];
			veto.IsUnderThreshold[k]=qdc.IsUnderThreshold[k];
			veto.IsUnderThreshold[16+k]=prevqdc.IsUnderThreshold[k];
			veto.IsOverflow[k]=qdc.IsOverflow[k];
			veto.IsOverflow[16+k]=prevqdc.IsOverflow[k];
		}
	}
	else if ((Int_t)prevqdc.card==-1) { 
		cout << "Previous Card was 0, EventMatch: " << EventMatch(i,run,nentries,qdc,prevqdc) << endl;
		veto.isGood=0;  
	}
	else { 
		printf("Failed to match!  Run: %i  VetoTree entry: %i  Card:%i  Prev.Card:%i  Breaking at %.0f%% through file.\n"
			,run,i,qdc.card,prevqdc.card,((Float_t)i/nentries)*100); 
		veto.isGood=-1;
	}

	if (verbose==1) {
	// check VEvent data
		printf("run:%i  vEnt:%i  sEnt:%i  card0:%i  card1:%i  card2:%i  sTime:%.5f  sTimeBad:%i  QDCTotal:%i\n",
			veto.run,veto.vEnt,veto.sEnt,veto.card[0],veto.card[1],veto.card[2],veto.sTime,veto.sTimeBad,veto.QDCTotal);
		cout << "QDC: "; for (int k = 0; k<numPanels; ++k) { cout << veto.QDC[k] << "  "; } cout << endl;
		cout << "UTh: "; for (int k = 0; k<numPanels; ++k) { cout << veto.IsUnderThreshold[k] << "  "; } cout << endl;
		cout << "Ovr: "; for (int k = 0; k<numPanels; ++k) { cout << veto.IsOverflow[k] << "  "; } cout << endl;
	}

	return veto;
}

// For tagging plane-based coincidences.
int PanelMap(int i){

	//0:L-Bot 1:U-Bot 2:Top 3:North 4:South 5:West 6:East

	if 		(i == 0) return 0;  // L-bot 1
	else if (i == 1) return 0;  
	else if (i == 2) return 0;
	else if (i == 3) return 0;
	else if (i == 4) return 0;
	else if (i == 5) return 0;  // L-bot 6

	else if (i == 6) return 1;  // U-bot 1
	else if (i == 7) return 1;
	else if (i == 8) return 1;
	else if (i == 9) return 1;
	else if (i == 10) return 1;
	else if (i == 11) return 1; // U-bot 6

	else if (i == 17) return 2; // Top 1
	else if (i == 18) return 2; // Top 2
	else if (i == 20) return 2; // Top 3
	else if (i == 21) return 2; // Top 4

	else if (i == 15) return 3; // North 1
	else if (i == 16) return 3; // North 2
	else if (i == 19) return 3; // North 3
	else if (i == 23) return 3; // North 4

	else if (i == 24) return 4; // South 1
	else if (i == 25) return 4; // South 2
	else if (i == 26) return 4; // South 3
	else if (i == 27) return 4; // South 4
	
	else if (i == 12) return 5; // West 1
	else if (i == 13) return 5; // West 2
	else if (i == 14) return 5; // West 3
	else if (i == 22) return 5; // West 4
	
	else if (i == 28) return 6; // East 1
	else if (i == 29) return 6; // East 2
	else if (i == 30) return 6; // East 3
	else if (i == 31) return 6; // East 4

	else return -1;
}

void ReadOutCoin(Bool_t ApproxTime, Float_t xTime, Float_t x_deltaT, Int_t planeHitCount[numPlanes], 
	vector<Int_t> candidate0, vector<Int_t> candidate1, vector<Int_t> candidate2) {

	if (ApproxTime) printf("Bad scaler time! Approximating ... event time:%.2f, time since last LED:%.2f\n"
		,xTime,x_deltaT);

	printf("Coincidence at %.2f sec, %.2f sec after last LED, hit %i panels.\n"
		,xTime,x_deltaT,(Int_t)candidate0.size());

	// read out: plane,panel,qdc,withinGaussian.
	printf("   Counts in planes: ");
	for (int q = 0; q < numPlanes; q++)
		printf("(%i):%i  ",q,planeHitCount[q]);
	cout << endl;
	vector<Int_t>::iterator it_c0 = candidate0.begin(); 
	vector<Int_t>::iterator it_c1 = candidate1.begin(); 
	vector<Int_t>::iterator it_c2 = candidate2.begin(); 
	while (it_c0 != candidate0.end()) {
		printf("   panel:%i  qdc:%i  within Gaussian:%i \n",*it_c0,*it_c1,*it_c2);
		++it_c0; ++it_c1; ++it_c2;
	}
}

//==============END PROCESSING FUNCTIONS=========================


// global pointers for qdc histograms.
TH1F *hRawQDC[numPanels];  
TH1F *hLEDCutQDC[numPanels];

void muFinder(string Input = ""){

	int mode = 0;		// switch: 0 for local files, 1 for pdsf files
	UInt_t card1 = 13;	// 11: prototype, 13: module 1
	UInt_t card2 = 18;
	Bool_t useThresh = true; // if true, also enables fitting LED peaks

	// "low" qdc threshold values
	//Int_t thresh[numPanels] = {123,115,95,93,152,115,105,103,119,91,108,103,94,107,95,167,
	//	53,150,89,120,65,85,132,62,130,101,80,108,145,164,119,82};

	// "high" qdc threshold values
	Int_t thresh[numPanels] = {124,117,96,93,155,115,112,105,120,91,109,108,95,112,96,168,
		63,157,100,127,72,100,140,65,145,125,82,112,151,168,122,94};

	// Input a list of run numbers
	Char_t InputName[200];
	string def = "Debug"; // default, if you don't specify an input at run time
	if (Input == "") strcpy(InputName,def.c_str());
	else strcpy(InputName,Input.c_str());
	Char_t InputFile[200];
	sprintf(InputFile,"%s.txt",InputName);
	ifstream InputList;
	InputList.open(InputFile);
	Char_t TheFile[200];

	// Input calibration file(s?) from builtVetoCal
	/* 
		Need to address:
		WHICH RUN RANGES a particular file applies to.
		Does the run range need to be inside the calibration table?  
		Should it be manually added here?
		Should it be added to the PeakInfo structure?

		Possible types of bad fits: 
	 	0. Fit values don't exist for a particular panel.
	 	1. Mean value is negative or greater than 4200
	 	2. Chi2NDF is > 15 ....... but 15 is kind of arbitrary
	 	3. Error in mean is > 4 ...... but 4 is kind of arbitrary
	 	4. builtVetoCal accidentally fits pedestal or gamma peak instead of LED.
	    	Check: Look at the output of BVC.  Does the pedestal in hThreshQDC# survive into hCutQDC#?
	*/
	const int cals = 1;
	string File[cals];
	File[0] = "./output/BVC_Weekend_CalibrationTable.txt";
	PeakInfo peaks[cals][numPanels];
	ifstream CalFiles;
	for (int i = 0; i<cals; i++) {
		CalFiles.open(File[i].c_str());
		cout << "Opened " << File[i] << endl;
		PeakInfo pk; // temporary peak
		while(!CalFiles.eof()){
			CalFiles >> pk.panel >> pk.mean >> pk.meanErr >> pk.sigma >> pk.sigmaErr >> pk.chi2_NDF;
			//printf("%i  %.0f  %.0f  %.0f  %.0f  %.0f\n",pk.panel,pk.mean,pk.meanErr,pk.sigma,pk.sigmaErr,pk.chi2_NDF);
			peaks[i][pk.panel].panel = pk.panel;
			peaks[i][pk.panel].mean = pk.mean;
			peaks[i][pk.panel].meanErr = pk.meanErr;
			peaks[i][pk.panel].sigma = pk.sigma;
			peaks[i][pk.panel].sigmaErr = pk.sigmaErr;
			peaks[i][pk.panel].chi2_NDF = pk.chi2_NDF;
	    }
	    CalFiles.close();
	}
	Int_t badFits[cals] = {0};
	for (int k = 0; k < numPanels; k++) {
		if (peaks[0][k].panel != k) {
			cout << "Missing peak entry! Panel: " << k << endl;
			peaks[0][k].badFit=true; badFits[0]++;
		}
		if ((peaks[0][k].mean < 0 || peaks[0][k].mean > 4200) || peaks[0][k].meanErr > 4 || peaks[0][k].chi2_NDF > 15) {
			peaks[0][k].badFit=true; badFits[0]++;
		}
	}
	if(badFits[0]>0) {
		cout << "Possible bad fits in " << badFits[0] << " panel(s):";
		for (int k = 0; k<numPanels; k++) { if (peaks[0][k].badFit==true) cout << " " << k; }
		cout << endl;
	}
	else cout << endl;

	//=== Global counters / variables / plots ===

		Float_t duration = 0;
		Float_t durationTotal = 0;
		
	  	// muon table variables (same structure as VEvent)
	  	VEvent veto;
    	Int_t run = 0;		// run number
		
		TH1D *TotalMultiplicity = new TH1D("TotalMultiplicity","Events over threshold",32,0,32);
	 	TotalMultiplicity->GetXaxis()->SetTitle("number of panels hit");

	 	const Int_t nqdc_bins=1400;  // this gives 3 qdc / bin
		const Float_t ll_qdc=0.;
		const Float_t ul_qdc=4200.;
		Char_t hname[50];
		for (Int_t i=0; i<numPanels; i++){
			sprintf(hname,"hRawQDC%d",i);
			hRawQDC[i] = new TH1F(hname,hname,nqdc_bins,ll_qdc,ul_qdc);
			sprintf(hname,"hLEDCutQDC%d",i);
			hLEDCutQDC[i] = new TH1F(hname,hname,500,ll_qdc,500);
		}

	 	// get number of files in dataset for any TGraphs needed
		Int_t filesToScan = 0;
		Int_t filesScanned = 0;
	  	while(!InputList.eof()) { InputList >> run; filesToScan++; }
	  	cout << "Scanning " << filesToScan << " files." << endl;	  	
	  	InputList.close();
	  	InputList.open(InputFile);
	  	run=0;

	  	TGraph *LEDFreq = new TGraph(filesToScan);
	 	LEDFreq->GetXaxis()->SetTitle("Frequency (Hz)");

		
	//=== End ===

	// Set up output file 
	Char_t OutputFile[200];
	sprintf(OutputFile,"./output/MF_%s.root",InputName);
	TFile *RootFile = new TFile(OutputFile, "RECREATE"); 	
  	TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."
	TTree *VTree = new TTree("VTree","Veto Panel Events");
	VTree->Branch("isGood",&veto.isGood,"isGood/I");
	VTree->Branch("run",&veto.run,"run/I");
	VTree->Branch("vEnt",&veto.vEnt,"vEnt/I");
	VTree->Branch("sEnt",&veto.sEnt,"sEnt/I");
	VTree->Branch("card",&veto.card,"card[3]/I");
	VTree->Branch("xTime",&veto.xTime,"xTime/F");
  	VTree->Branch("sTime",&veto.sTime,"sTime/F");
  	VTree->Branch("sTimeBad",&veto.sTimeBad,"sTimeBad/O"); 
  	VTree->Branch("lTime",&veto.lTime,"lTime/F");    
  	VTree->Branch("eTime",&veto.eTime,"eTime/F");    
  	VTree->Branch("QDCTotal",&veto.QDCTotal,"QDCTotal/I");
  	VTree->Branch("QDC",&veto.QDC,"QDC[32]/I");
  	VTree->Branch("IsUnderThreshold",&veto.IsUnderThreshold,"IsUnderThreshold[32]/I");
  	VTree->Branch("IsOverflow",&veto.IsOverflow,"IsOverflow[32]/I");
  	VTree->Branch("coinType",&veto.coinType,"coinType[10]/O");
  	VTree->Branch("cutType",&veto.cutType,"cutType[10]/O");
  	VTree->Branch("plane",&veto.plane,"plane[7]/O");
  	VTree->Branch("planeHitCount",&veto.planeHitCount,"planeHitCount[7]/I");

	// Loop over files
	while(!InputList.eof()){

		// initialize 
		InputList >> run;
		if (mode==0) sprintf(TheFile,"~/dev/datasets/muFinder/OR_run%i.root",run);
		else if (mode==1) sprintf(TheFile,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/P3JDY/OR_run%u.root",run); 
		TChain *VetoTree = new TChain("VetoTree");
		VetoTree->AddFile(TheFile);
		TChain *MGTree = new TChain("MGTree");
		MGTree->AddFile(TheFile);
		MJTRun *MyRun = new MJTRun();
		MGTree->SetBranchAddress("run",&MyRun);
		Long64_t nentries = VetoTree->GetEntries();
		MGTBasicEvent *b = new MGTBasicEvent; 
		VetoTree->SetBranchAddress("vetoEvent",&b);
		MJTVetoData *vd[vd_size];
		VetoTree->GetEntry(0);
		for (int i=0; i<vd_size; i++) { vd[i] = dynamic_cast<MJTVetoData*>(b->GetDetectorData()->At(i)); }
        MGTree->GetEntry(0);
        duration = MyRun->GetStopTime() - MyRun->GetStartTime();
        if (duration < 0) {
        	printf("\nRun %i has duration %.0f, skipping file!\n\n",run,duration);
        	continue;
        }
        durationTotal += duration;

    	//=== Single-file counters / variables / plots

		QEvent qdc, prevqdc;
		prevqdc.EventCount=-1;
		VEvent prevLED, prevEvent;

		Int_t corruptCount = 0;

		Int_t multiplicity = 0;
		Int_t qdc_within_Gaussian = 0;
		Int_t underthresh_but_within_Gaussian = 0;
		Int_t qdc_outside_Gaussian = 0;
		Int_t LEDCount = 0;
		Int_t LEDTotal = 0;

		Int_t missed = 0;	// LEDs which are tagged correctly but have the wrong delta-T
		Int_t expected = 0; // expected LEDs

		Float_t firstLEDTime = 0;
		Float_t e_prevTime = 0;
		Float_t e_deltaT = 0;
		Float_t l_prevTime = 0;
		Float_t l_deltaT = 0;
		Float_t s_prevTime = 0;
		Float_t s_deltaT = 0;

		vector<Float_t> s_deltaTs;
		vector<Int_t> underthresh[2]; // 0: panel, 1: qdc entry.  keep track of which panels are underthresh for a given VEvent

		vector<Int_t> candidate[3]; // keep track of: 1: panel number, 2: qdc value, 3: value within Gaussian

		
		sprintf(hname,"LEDDeltaT_run%i",run);
		TH1F *LEDDeltaT = new TH1F(hname,hname,100000,0,100); // 0.001 sec/bin
		Bool_t WriteLEDDeltaT = true; // set true to write histos to ROOT file.
		Int_t maxbin = 0;
		Float_t rms = 0;	// error in frequency estimation due to NATURAL DRIFT

		sprintf(hname,"lTime_vs_sTime_run%i",run);
		TH1F *lTime_vs_sTime = new TH1F(hname,hname,100,-0.05,0.05); // 0.001 sec/bin
		lTime_vs_sTime->GetXaxis()->SetTitle("LEDTime-ScalerTime (sec)");

		sprintf(hname,"eTime_vs_sTime_run%i",run);
		TH1F *eTime_vs_sTime = new TH1F(hname,hname,6000,-30,30); // 0.01 sec/bin
		eTime_vs_sTime->GetXaxis()->SetTitle("EntryTime-ScalerTime (sec)");

		Float_t LEDDeltaT_avg = 0;
		Float_t LEDDeltaT_hist = 0;

		//=== End ===

		// ========= 1st loop over VetoTree entries - Find & count LED's. =========

		printf("====================\nNow scanning run %i: %lli entries, %.2f sec. \n",run,nentries,duration);
		for (int i = 0; i < nentries; i++) {
			VetoTree->GetEntry(i);
			qdc = WriteQEvent(0,run,i,vd);

			if (EventMatch(i,run,nentries,qdc,prevqdc)) {
				veto = WriteVEvent(0,i,run,nentries,duration,card1,card2,qdc,prevqdc,vd);
				if (veto.isGood==0) continue; // previous card was 0
				if (veto.isGood==-1) break;  // failed to match
				
				//====================BEGIN ANALYSIS=================

				// multiplicity is usually 31 for LED's because Card 13, Ch. 15 is broken. (7/2015)
				// requiring the qdc_within variable to be within 5 of multiplicity:
				// this is trying to tag a "pileup" (muon+LED) event.

				// check scaler corruption 
				if (veto.sTimeBad) corruptCount++;

				// loop over the VEvent
				// flag events within +/- 3 sigma of each panel's LED peak 
				for (int k=0;k<numPanels;k++) {
					if (veto.QDC[k]>thresh[k]) {
						multiplicity++;
					}
					else if (veto.QDC[k] < thresh[k]
						&& peaks[0][k].badFit==false
						&& veto.QDC[k] < (peaks[0][k].mean + peaks[0][k].sigma*3) 
						&& veto.QDC[k] > (peaks[0][k].mean - peaks[0][k].sigma*3)) {
						underthresh_but_within_Gaussian++;
						underthresh[0].push_back(k);
						underthresh[1].push_back(veto.QDC[k]);
					}
					if (peaks[0][k].badFit==false
						&& veto.QDC[k] < (peaks[0][k].mean + peaks[0][k].sigma*3) 
						&& veto.QDC[k] > (peaks[0][k].mean - peaks[0][k].sigma*3)) {
						qdc_within_Gaussian++;
					}
				}
				// **** CLINT'S LED CONDITION ****
				// 1. Basic multiplicity cut: this knocks out 2*T's but doesn't use any peak info
				//if (multiplicity >= numPanels-4) { 
				
				// 2. Add underthresh entries back into multiplicity, but don't go over 31 (broken qdc channel)
				if (underthresh_but_within_Gaussian > 0 
					&& multiplicity > numPanels-4
					&& multiplicity < numPanels-1) 
					multiplicity += underthresh_but_within_Gaussian;

				if (multiplicity >= numPanels-4) { // 32,31,30,29,28
				// ******* IT'S SO AWESOME *******
					
					/*
					// Debug adding underthresh values back to multiplicity
					// This is likely to grab a "pedestal" value if the Gaussian extends below the SW threshold.
					if (underthresh_but_within_Gaussian>0) {
						printf("\n  LED Event: final multiplicity:%i  underthresh:%i  original multiplicity:%i  Panels underthresh:"
							,multiplicity,underthresh_but_within_Gaussian,multiplicity-underthresh_but_within_Gaussian);
						for (vector<Int_t>::iterator it1 = underthresh[0].begin(); it1 != underthresh[0].end(); ++it1) {
							printf(" #%i ",*it1);
						}
						for (vector<Int_t>::iterator it2 = underthresh[1].begin(); it2 != underthresh[1].end(); ++it2) {
							printf(" qdc: %i ",*it2);
						}
						//cout << endl;
					}
					*/
					e_deltaT = veto.eTime-e_prevTime;
					s_deltaT = veto.sTime-s_prevTime;
					s_deltaTs.push_back(s_deltaT);
					//printf("\n  LED event: -- multip:%i  gaussian:%i\n  eTime:%.4f  prevETime:%.4f  e_deltaT(don't trust):%.4f \n  sTime:%.4f  prevSTime:%.4f  s_deltaT:%.8f\n"
					//	,multiplicity,qdc_within_Gaussian,veto.eTime,e_prevTime,e_deltaT,veto.sTime,s_prevTime,s_deltaT);
					s_prevTime = veto.sTime;
					e_prevTime = veto.eTime;
					e_deltaT=0; s_deltaT=0;
				}
				else if (multiplicity < numPanels-4 && multiplicity > numPanels-7) { //28,27,26,25
					e_deltaT = veto.eTime-e_prevTime;
					s_deltaT = veto.sTime-s_prevTime;
					printf("\n  Possible LED event (not counted): -- final multip:%i  underthresh:%i\n  eTime:%.4f  prevETime:%.4f  e_deltaT(don't trust):%.4f \n  sTime:%.4f  prevSTime:%.4f  s_deltaT:%.4f \n\n"
						,multiplicity,underthresh_but_within_Gaussian,veto.eTime,e_prevTime,e_deltaT
						,veto.sTime,s_prevTime,s_deltaT);
					s_prevTime = veto.sTime;
					e_prevTime = veto.eTime;
					e_deltaT=0; s_deltaT=0;
				}
				//else 
					//printf("\n Non-LED event: -- multip:%i  gaussian:%i  eTime:%.4f  e_deltaT:%.4f  sTime:%.4f  s_deltaT:%.4f" 
					//	,multiplicity,qdc_within_Gaussian,veto.eTime[0],veto.eTime[0]-e_prevTime
					//	,veto.sTime[0],veto.sTime[0]-s_prevTime);
				
				// Reset & save before getting next VEvent
				underthresh[0].erase(underthresh[0].begin(),underthresh[0].end());
				underthresh[1].erase(underthresh[1].begin(),underthresh[1].end());
				multiplicity=0;
				qdc_within_Gaussian=0;
				underthresh_but_within_Gaussian=0;
				badFits[0]=0;

				//=====================END ANALYSIS===================

			} // end EventMatch condition

			// Save qdc into prevqdc before getting next VetoTree entry.
			prevqdc = qdc;

		}	// ========= End 1st loop over VetoTree entries. =========

		// ============= 1st loop output ===============

		// 1. Assuming the 250MHz scaler is limited by the 100MHz trigger card (10ns steps)
		// 	  the most precision we can get out of sTime is 8 decimal places.
		// 2. Comparing the "average method" with the "histo method" 
		// 	  looks like a good way to tell if we've found all LEDs.

		printf(" Scaler Corruption: %i entries, %.2f%%.\n",corruptCount,((Float_t)corruptCount/nentries)*100);

		LEDTotal = s_deltaTs.size();
		LEDDeltaT_avg = LEDTotal/duration;
		printf(" LED Results: Found %i LEDs.\n  1. Averaging method (only for comparison): %.8fs between LEDs, %.4f Hz. \n"
			,LEDTotal,1/((LEDTotal)/duration),LEDDeltaT_avg);
		
		for (vector<Float_t>::iterator it = s_deltaTs.begin(); it != s_deltaTs.end(); ++it) {
    		LEDDeltaT->Fill(*it);
		}
		sprintf(hname,"LEDDeltaT_run%i",run);
    	if(WriteLEDDeltaT) LEDDeltaT->Write(hname,TObject::kOverwrite); // currently writes a new one for each file.
		maxbin = LEDDeltaT->GetMaximumBin();
		LEDDeltaT->GetXaxis()->SetRange(maxbin-10,maxbin+10);
		rms = LEDDeltaT->GetRMS();
		LEDDeltaT_hist = LEDDeltaT->GetMean();
		printf("  2. Delta-T histo method (used to calculate LED time):\n");
		printf("     Max bin:%i  Max bin entries:%.0f, %.8fs between LEDs (RMS:%.8fs), %.4f Hz. \n"
			,maxbin,LEDDeltaT->GetBinContent(maxbin),LEDDeltaT_hist,rms,1/LEDDeltaT_hist);

		cout << "     Non-corrupted delta-t's outside RMS error: ";
		for (vector<Float_t>::iterator it = s_deltaTs.begin(); it != s_deltaTs.end(); ++it) {
			if (fabs(*it-LEDDeltaT_hist)>rms) { 
				if (*it > 0 && *it < 100) printf(" %.4f ",*it); 
				missed++; 
			}
		}
		printf("\n     %i LEDs (counting corrupted entries) had error larger than RMS, %.2f%%.\n",missed,((Float_t)missed/nentries)*100);
		
		expected = (1/LEDDeltaT_hist)*(Int_t)duration;
		printf("     Expect %i (+/-1) LEDs from hist method, found %i.\n",expected,LEDTotal);
		if (abs(expected-LEDTotal) > 1) 
			printf("  Warning!  Possibly missed %i LED's.\n",abs(expected-LEDTotal));

		LEDFreq->SetPoint(filesScanned,run,1/LEDDeltaT_hist);

    	// ========= end of 1st loop output ============


		// ========= 2nd loop over VetoTree entries - find some damn muons! =========

		prevqdc.EventCount=-1;  // re-initialize
		LEDCount=0;
		printf("\n Looking for muon candidates .... \n Planes are: 0:L-Bot 1:U-Bot 2:Top 3:North 4:South 5:West 6:East\n\n");
		for (int i = 0; i < nentries; i++) {
			VetoTree->GetEntry(i);
			qdc = WriteQEvent(0,run,i,vd);

			if (EventMatch(i,run,nentries,qdc,prevqdc)) {
				veto = WriteVEvent(0,i,run,nentries,duration,card1,card2,qdc,prevqdc,vd);
				if (veto.isGood==0) continue; // previous card was 0
				if (veto.isGood==-1) break;  // failed to match
				
				//=====================BEGIN ANALYSIS=================
				
				// loop over the VEvent 
				// flag events within +/- 3 sigma of each panel's LED peak 
				for (int k=0;k<numPanels;k++) {
					if (veto.QDC[k]>thresh[k]) {
						multiplicity++;
					}
					else if (veto.QDC[k] < thresh[k]
						&& peaks[0][k].badFit==false
						&& veto.QDC[k] < (peaks[0][k].mean + peaks[0][k].sigma*3) 
						&& veto.QDC[k] > (peaks[0][k].mean - peaks[0][k].sigma*3)) {
						underthresh_but_within_Gaussian++;
					}
					if (peaks[0][k].badFit==false
						&& veto.QDC[k] < (peaks[0][k].mean + peaks[0][k].sigma*3) 
						&& veto.QDC[k] > (peaks[0][k].mean - peaks[0][k].sigma*3)) {
						qdc_within_Gaussian++;
					}
				}

				// Use same LED condition as 1st scan, and adjust multiplicity of high-multiplicity events:
				if (underthresh_but_within_Gaussian > 0 
					&& multiplicity > numPanels-4
					&& multiplicity < numPanels-1) 
					multiplicity += underthresh_but_within_Gaussian;
				
				if (multiplicity >= numPanels-4) { // 32,31,30,29,28
					LEDCount++;
					if (LEDCount==1 && !veto.sTimeBad) firstLEDTime = veto.sTime; // save first LED time, it's always less than T
					else if (LEDCount==1 && veto.sTimeBad) printf("Warning! Failed to find first LED time from scaler due to corruption!\n");
				}

				// Find LED time (veto.lTime) of the VEvent (needed if sTime is bad)
				if (LEDCount==1) veto.lTime = firstLEDTime;
				if (LEDCount>1) veto.lTime = firstLEDTime + LEDDeltaT_hist * (LEDCount-1);
				//printf(" M:%-3i   LEDs:%-5i   LEDTime:%-10.4f   sTime:%-10.4f  diff:%-5.4f\n"
				//	,multiplicity,LEDCount,veto.lTime,veto.sTime,veto.lTime-veto.sTime);
				eTime_vs_sTime->Fill(veto.eTime-veto.sTime);
				
				// store previous LED VEvent for time comparison
				if (multiplicity >= numPanels-4) { // 32,31,30,29
					prevLED = veto;  
					lTime_vs_sTime->Fill(veto.lTime-veto.sTime);
				}

				// PanelMap outputs:
				// 0:L-Bot 1:U-Bot 2:Top 3:North 4:South 5:West 6:East
				for (int k = 0; k<numPanels; k++) {
					
					if (veto.QDC[k]>thresh[k]){
						if (PanelMap(k)==0) { veto.plane[0]=true; veto.planeHitCount[0]++; }
						else if (PanelMap(k)==1) { veto.plane[1]=true; veto.planeHitCount[1]++; }
						else if (PanelMap(k)==2) { veto.plane[2]=true; veto.planeHitCount[2]++; }
						else if (PanelMap(k)==3) { veto.plane[3]=true; veto.planeHitCount[3]++; }
						else if (PanelMap(k)==4) { veto.plane[4]=true; veto.planeHitCount[4]++; }
						else if (PanelMap(k)==5) { veto.plane[5]=true; veto.planeHitCount[5]++; }
						else if (PanelMap(k)==6) { veto.plane[6]=true; veto.planeHitCount[6]++; }
						
						candidate[0].push_back(k); // save the panel number
						candidate[1].push_back(veto.QDC[k]);
						if (peaks[0][k].badFit==false
							&& veto.QDC[k] < (peaks[0][k].mean + peaks[0][k].sigma*3) 
							&& veto.QDC[k] > (peaks[0][k].mean - peaks[0][k].sigma*3)) 
							candidate[2].push_back(1);  // value is within 3 sigma of LED peak
						else candidate[2].push_back(0); // value is not within 3 sigma or fit is bad.
					}
				}
				
				// 1. Time cut:
				// Throw away events with a delta-t that matches the period calculated in the last section.
				Bool_t TimeCut = false;

				// get time information of last LED
				l_deltaT = veto.lTime-prevLED.lTime;
				s_deltaT = veto.sTime-prevLED.sTime;
				l_prevTime = prevLED.lTime;
				s_prevTime = prevLED.sTime;

				// get entry time information of last event (this is less than LED interval)
				// used when current veto.sTime is bad.
				e_deltaT = veto.eTime-prevEvent.eTime;
				e_prevTime = prevEvent.eTime;

				// placeholder time if sTime is corrupted.
				Float_t xTime = 0; 
				Float_t x_deltaT = 0; 
				Bool_t ApproxTime = false;
				if (!veto.sTimeBad && !prevLED.sTimeBad) {
					veto.xTime = veto.sTime;
					x_deltaT = s_deltaT;
				}
				else if (!veto.sTimeBad && prevLED.sTimeBad) {
					veto.xTime = veto.sTime;
					x_deltaT = veto.sTime-prevLED.lTime;
					ApproxTime = true;
					//printf("Bad scaler time. xTime = veto.sTime (OK), x_deltaT = veto.sTime - prevLED.lTime : %.2f = %.2f - %.2f\n"
					//	,x_deltaT,veto.sTime,prevLED.lTime);
				}
				else if (veto.sTimeBad && e_deltaT < LEDDeltaT_hist) { 
					veto.xTime = veto.lTime+e_deltaT;
					x_deltaT = e_deltaT;
					ApproxTime = true;
					//printf("Bad scaler time. xTime = lTime + e_deltaT : %.2f = %.2f + %.2f, x_deltaT = e_deltaT : %.2f\n"
					//	,xTime,veto.lTime,e_deltaT,e_deltaT);
				}
				else if (veto.sTimeBad && e_deltaT > LEDDeltaT_hist) {
					veto.xTime = veto.lTime;
					x_deltaT = l_deltaT;
					ApproxTime = true;
					//printf("Bad scaler time! xTime = lTime : %.2f, x_deltaT = l_deltaT (won't pass time cut) : %.2f\n"
					//	,xTime,l_deltaT);
				}
				else{
					printf("Bad scaler time. Warning: Failed all test conditions!\n");
				}
				if (fabs(LEDDeltaT_hist-x_deltaT)>rms) TimeCut = true;


				// 2. "QDC within Gaussian" cut:
				// Will apply only to high-multiplicity events - say, M>20.
				// If more than half the panels are within their gaussian, figure it's an LED.
				// Once we know locations of muon peaks in panels, can use them to tighten this cut
				// to exclude single gammas.
				Bool_t GaussianCut = false;
				if (multiplicity>20 && qdc_within_Gaussian>multiplicity/2) {
					GaussianCut = false;
				}
				else GaussianCut = true; 


				// 3. "QDC is outside Gaussian" cut:
				// AT LEAST ONE of the panels in THIS CANDIDATE event should 
				// have a QDC value outside its LED's Gaussian.
				Bool_t EventOutsideGaussianCut = false;
				vector<Int_t>::iterator it_g0 = candidate[0].begin(); 
				vector<Int_t>::iterator it_g1 = candidate[1].begin(); 
				vector<Int_t>::iterator it_g2 = candidate[2].begin(); 
				while (it_g0 != candidate[0].end()) {
					if (*it_g2==0) {
						//printf(" EOGC: time:%.2f  Approx:%i  panel:%i  qdc:%i  within Gaussian:%i\n"
						//	,xTime,ApproxTime,*it_g0,*it_g1,*it_g2);
						EventOutsideGaussianCut=true;
						qdc_outside_Gaussian++;
					}
					++it_g0; ++it_g1; ++it_g2;
				}
				
				/*
				// debug: check candidate events.
				if (!TimeCut){
					printf("Debug: Candidate event is:\n");
					ReadOutCoin(ApproxTime,xTime,x_deltaT,planeHitCount,candidate[0],candidate[1],candidate[2]);
					printf("Bools: ApproxTime:%i  TimeCut:%i  GaussianCut:%i  EventOutsideGaussianCut:%i\n",ApproxTime,TimeCut,GaussianCut,EventOutsideGaussianCut);
					printf("Times: xTime:%.2f  x_deltaT:%.2f  TimeCut(s):%.2f\n\n",xTime,x_deltaT,fabs(LEDDeltaT_hist-x_deltaT));
				}
				*/

				// applying time cut would exclude a muon that had the LED Delta T.
				if (TimeCut && GaussianCut && EventOutsideGaussianCut) {

					// 0:L-Bot  1:U-Bot  2:Top  3:North  4:South  5:West  6:East

					// 1. (Yuri's favorite)
					if (veto.plane[0] && veto.plane[1] && veto.plane[2]) {
						veto.coinType[1]=true;
						printf("3-Plane coincidence: Top + U-Bot + L-Bot, %i panels outside Gaussian\n",qdc_outside_Gaussian);
						ReadOutCoin(ApproxTime,veto.xTime,x_deltaT,veto.planeHitCount,candidate[0],candidate[1],candidate[2]);
						cout << endl;
					}
					
					// 2. (Side + U-Bot + L-Bot)
					if (veto.plane[0] && veto.plane[1] && (veto.plane[3]||veto.plane[4]||veto.plane[5]||veto.plane[6])) {
						veto.coinType[2]=true;
						printf("3-Plane coincidence: Side + U-Bot + L-Bot, %i panels outside Gaussian\n",qdc_outside_Gaussian);
						ReadOutCoin(ApproxTime,veto.xTime,x_deltaT,veto.planeHitCount,candidate[0],candidate[1],candidate[2]);
						cout << endl;
					}

					// 3. "2-plane" (any two sides) + qdc_outside_Gaussian >= 1
					if (((veto.plane[0] && veto.plane[1]) || (veto.plane[0] && veto.plane[2]) || (veto.plane[0] && veto.plane[3])
						|| (veto.plane[0] && veto.plane[4]) || (veto.plane[0] && veto.plane[5]) || (veto.plane[0] && veto.plane[6])
						|| (veto.plane[1] && veto.plane[2]) || (veto.plane[1] && veto.plane[3]) || (veto.plane[1] && veto.plane[4])
						|| (veto.plane[1] && veto.plane[5]) || (veto.plane[1] && veto.plane[6]) || (veto.plane[2] && veto.plane[3])
						|| (veto.plane[2] && veto.plane[4]) || (veto.plane[2] && veto.plane[5]) || (veto.plane[2] && veto.plane[6])
						|| (veto.plane[3] && veto.plane[4]) || (veto.plane[3] && veto.plane[5]) || (veto.plane[3] && veto.plane[6])
						|| (veto.plane[4] && veto.plane[5]) || (veto.plane[4] && veto.plane[6]) || (veto.plane[5] && veto.plane[6]))) {
						veto.coinType[3]=true;
						//printf("2-side coincidence, %i panels outside Gaussian\n",qdc_outside_Gaussian);
						//ReadOutCoin(ApproxTime,veto.xTime,x_deltaT,planeHitCount,candidate[0],candidate[1],candidate[2]);
						//cout << endl;
					}

				}

				// Fill VTree
				if (TimeCut) veto.cutType[0]=true;
				if (GaussianCut) veto.cutType[1]=true;
				if (EventOutsideGaussianCut) veto.cutType[2]=true;
				if (ApproxTime) veto.cutType[3]=true;
				VTree->Fill();

				//=====================END ANALYSIS===================

				// reset & save before getting next VEvent
				prevEvent = veto;
				for (int c = 0; c<3; c++) candidate[c].erase(candidate[c].begin(),candidate[c].end());

			} // end EventMatch condition

			// Save qdc into prevqdc before getting next VetoTree entry.
			prevqdc = qdc;
			multiplicity=0;
			qdc_within_Gaussian=0;
			qdc_outside_Gaussian=0;
			underthresh_but_within_Gaussian=0;

		}	// End loop over VetoTree entries.

		// ===== 2nd loop output ====

			sprintf(hname,"lTime_vs_sTime_run%i",run);
			lTime_vs_sTime->Write(hname,TObject::kOverwrite);

			sprintf(hname,"eTime_vs_sTime_run%i",run);
			eTime_vs_sTime->Write(hname,TObject::kOverwrite);

		// ==========================

		delete VetoTree;
		delete MGTree;
		filesScanned++;
	} // End loop over files.

	

	// === END OF SCAN Output & Plotting ===
	printf("\n====================\nFinished loop over files.\n");
	
	TotalMultiplicity->Write("TotalMultiplicity",TObject::kOverwrite);

	TCanvas *c1 = new TCanvas("c1", "Bob Ross's Canvas",600,600);
	c1->SetGrid();
	LEDFreq->SetMarkerColor(4);
	LEDFreq->SetMarkerStyle(21);
	LEDFreq->SetMarkerSize(0.5);
	LEDFreq->SetTitle("LED Frequency (histo method)");
	LEDFreq->GetXaxis()->SetTitle("Frequency (Hz)");
	LEDFreq->Draw("ALP");	
	LEDFreq->Write("LEDFreq",TObject::kOverwrite);


	// Output canvasses of interest.
	Char_t OutputName[200];
	//sprintf(OutputName,"./output/MF_%s_VetoQDC.C",InputName);
	//vcan0->Print(OutputName);

	// ==========================
	VTree->Write();
	RootFile->Close();
	cout << "Wrote ROOT file." << endl;
}

// This should be able to take multiple inputs, but I can't figure out the ROOT syntax:
// .x muFinder.C++ ("M1BG_debug" "M1BG_Range5")  FAILS.
#ifndef __CINT__
int main(int argc, const char* argv[]) {

	//gROOT->ProcessLine(".x /Users/wisecg/Applications/mgsw/MGDO/Majorana/LoadMGDOMJClasses.C");
	cout << "made it to main arg" << endl;
	
	size_t size;
	if (argc>0) {
		for(int i = 0; i < argc; i++) {
			printf( "Input arg %d: %s\n", i, argv[i] );
			size = sizeof argv[i]/sizeof(size_t);
			string Input(argv[i],size);
			//muFinder(Input);
		}
	}
	else {
		string Input = "";
		cout << "Input: " << Input << endl;
		//muFinder(Input);
	}
}
#endif
