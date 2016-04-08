// Finds muons.  Optionally writes some output.
// Clint Wiseman, USC/Majorana
// 3/9/2016

#include "vetoScan.hh"

using namespace std;

void muFinder(string Input, int *thresh, bool root, bool list)
{
	// Custom SW Threshold (obtained from vetoThreshFinder)
	int swThresh[32] = {0};
	if (thresh != NULL) {
		cout << "muFinder is using these SW thresholds: " << endl;
		memcpy(swThresh,thresh,sizeof(swThresh));
		for (int j=0;j<32;j++) cout << j << ":" << swThresh[j] << " ";
		cout << endl;
	}
	else {
		cout << "thresh is NULL!  Using default." << endl;
		for (int j=0;j<32;j++) swThresh[j] = 500;
	}

	// Input a list of run numbers
	ifstream InputList(Input.c_str());
	if(!InputList.good()) {
		cout << "Couldn't open " << Input << endl;
		return;
	}

	// Set up output files
	string Name = Input;
	Name.erase(Name.find_last_of("."),string::npos);
	Name.erase(0,Name.find_last_of("\\/")+1);

	// Output 1: Text file muon list (used in skim files)
	string outName = "./output/MuonList_"+Name+".txt";
	ofstream MuonList(outName.c_str());
	if (!list) MuonList.close();

	// Output 2: ROOT output
	int isGood;
	int run;
	long start;
	long stop;
	long prevStopTime = 0;
	double d;
	int CoinType[32];
	int CutType[32];
	int highestMultip = 0;
	double LEDfreq = 0;
	double LEDrms = 0;
	double xTime = 0;
	double x_deltaT = 0;	
	Char_t OutputFile[200];
	sprintf(OutputFile,"./output/%s.root",Name.c_str());
	TFile *RootFile = new TFile(OutputFile, "RECREATE"); 	
  	TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."
	TTree *vetoEvent = new TTree("vetoEvent","MJD Veto Events");
	MJVetoEvent out;
	if (root) {
		vetoEvent->Branch("events","MJVetoEvent",&out,32000,1);
		vetoEvent->Branch("CoinType[32]",CoinType,"CoinType[32]/I");
		vetoEvent->Branch("CutType[32]",CutType,"CutType[32]/I");
		vetoEvent->Branch("LEDfreq",&LEDfreq);
		vetoEvent->Branch("LEDrms",&LEDrms);
		vetoEvent->Branch("highestMultip",&highestMultip);
		vetoEvent->Branch("start",&start,"start/L");
		vetoEvent->Branch("stop",&stop,"stop/L");
		vetoEvent->Branch("xTime",&xTime);
		vetoEvent->Branch("x_deltaT",&x_deltaT);
	}

	// Loop over files.
	while(!InputList.eof()){

		// initialize
		InputList >> run;
		GATDataSet *ds = new GATDataSet(run);
		TChain *v = ds->GetVetoChain();
		long vEntries = v->GetEntries();
		MJTRun *vRun = new MJTRun();
		MGTBasicEvent *vEvent = new MGTBasicEvent(); 
		unsigned int mVeto = 0;
		uint32_t vBits = 0;
		v->SetBranchAddress("run",&vRun);
		v->SetBranchAddress("mVeto",&mVeto);
		v->SetBranchAddress("vetoEvent",&vEvent);
		v->SetBranchAddress("vetoBits",&vBits);
		v->GetEntry(0);
		start = (long)vRun->GetStartTime();
		stop = (long)vRun->GetStopTime();
		d = (double)(stop-start);

		printf("\n=========== Scanning Run %i: %li entries. ===========\n",run,vEntries);
		printf("start: %li  stop: %li  d: %.0f\n",start,stop,d);

		
		// ========= 1st loop over veto entries - Measure LED frequency. =========
		// 
		// Goal is to measure the LED frequency, to be used in the second loop as a 
		// time cut. This is done by finding the maximum bin of a delta-t histogram.
		// This section of the code uses a weak multiplicity threshold of 20 -- it 
		// doesn't need to be exact, and should also work for runs where there were 
		// only 24 panels installed.
		// 
		bool badLEDFreq = false;
		MJVetoEvent prev;
		char hname[200];
		sprintf(hname,"LEDDeltaT_run%i",run);
		TH1F *LEDDeltaT = new TH1F(hname,hname,100000,0,100); // 0.001 sec/bin
		highestMultip = 0;	// try to predict how many panels there are for this run.
		long skippedEvents = 0;
		long corruptScaler = 0;
		bool foundFirst = false;
		int firstGoodEntry = 0;
		MJVetoEvent first;
		highestMultip=0;
		vector<double> LocalEntryTime;
		vector<double> LocalEntryNum;
		vector<bool> LocalBadScalers;	
		double lastGoodTime = 0;
		int pureLEDcount = 0;
		for (long i = 0; i < vEntries; i++) 
		{
			v->GetEntry(i);
			MJVetoEvent veto;
			veto.SetSWThresh(swThresh);	
	    	isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true);
    		
    		// find event time and fill vectors
			if (!veto.GetBadScaler()) {
				LocalBadScalers.push_back(0);
				xTime = veto.GetTimeSec();
			}
			else {
				LocalBadScalers.push_back(1);
				corruptScaler++;
				xTime = ((double)i / vEntries) * d;
			}
			
	    	// fill vectors
	    	// (the time vectors are revised in the second loop)
			LocalEntryNum.push_back(i);		
			LocalEntryTime.push_back(xTime);

			// skip bad entries (true = print contents of skipped event)
	    	if (CheckForBadErrors(veto,i,isGood,false)) continue;

	    	if (veto.GetMultip() > highestMultip && veto.GetMultip() < 33) {
	    		highestMultip = veto.GetMultip();
	    		cout << "Finding highest multiplicity: " << highestMultip << "  entry: " << i << endl;
	    	}

	    	// Save the first good entry number for the SBC offset
			if (isGood == 1 && !foundFirst) {
				first = veto;
				foundFirst = true;
				firstGoodEntry = i;
			}
			
	    	// Very simple LED tag.
			if (veto.GetMultip() >= 20) {
				LEDDeltaT->Fill(veto.GetTimeSec()-prev.GetTimeSec());
				pureLEDcount++;
			}
			
			// end of loop : save things
			prev = veto;
			lastGoodTime = xTime;
		}

		// make sure the local vectors are all the same size
		if (LocalEntryNum.size() != LocalEntryTime.size())
		printf("Warning! Local vectors are not the same size!\n");

		// if duration is corrupted, use the last good timestamp as the duration.
		if (d == 0) {
			printf("Corrupted duration. Using last good timestamp: %.2f\n",lastGoodTime);
			d = lastGoodTime;
		}

		// Find the SBC offset		
		double SBCOffset = first.GetTimeSBC() - first.GetTimeSec();
		printf("First good entry: %i  SBCOffset: %.2f\n",firstGoodEntry,SBCOffset);

		// Find the LED frequency	
		if (skippedEvents > 0) printf("Skipped %li of %li entries.\n",skippedEvents,vEntries);
		if (corruptScaler > 0) printf("Corrupt scaler: %li of %li entries (%.2f%%) .\n"
			,corruptScaler,vEntries,100*(double)corruptScaler/vEntries);
		
		bool LEDTurnedOff = false;
		if (highestMultip < 20) {
			printf("Warning!  LED's may be off!\n");
			LEDTurnedOff = true;
		}
		LEDrms = 0;
		LEDfreq = 0;
		int dtEntries = LEDDeltaT->GetEntries();
		if (dtEntries > 0) {
			int maxbin = LEDDeltaT->GetMaximumBin();
			LEDDeltaT->GetXaxis()->SetRange(maxbin-100,maxbin+100); // looks at +/- 0.1 seconds of max bin.
			LEDrms = LEDDeltaT->GetRMS();
			if (LEDrms==0) LEDrms = 0.1;
			LEDfreq = 1/LEDDeltaT->GetMean();
		}
		else {
			printf("Warning! No multiplicity > 20 events!!\n");
			LEDrms = 9999;
			LEDfreq = 9999;
			LEDTurnedOff = true;
		}
		double LEDperiod = 1/LEDfreq;
		delete LEDDeltaT;

		// Ralph suggests we vary these parameters and do a study of the accidentals.
		// double RMSTimeWindow = 10 * LEDrms;
		double RMSTimeWindow = 0.1;
		int highMultipThreshold = highestMultip - 10;
		printf("HM: %i LED_f: %.8f LED_t: %.8f RMS: %8f\n",highestMultip,LEDfreq,1/LEDfreq,LEDrms);
		printf("LED window: %.2f  Multip Threshold: %i\n",RMSTimeWindow,highMultipThreshold);
		
		// set a flag for "bad LED" (usually a short run causes it)
		// and replace the period with the "simple" one if possible
		badLEDFreq = false;
		if (LEDperiod > 9 || vEntries < 100) 
		{
			printf("Warning: Short run.\n");
			if (pureLEDcount > 3) {
				printf("   From histo method, LED freq is %.2f.\n   Reverting to the approx rate (%.2fs) ... \n"
					,LEDfreq,(double)pureLEDcount/d);
				LEDperiod = d/pureLEDcount;
			}
			else { 
				printf("   Warning: LED info is corrupted!  Will not use LED period information for this run.\n");
				LEDperiod = 9999;
				badLEDFreq = true;
			}
		}

		// ========= 2nd loop over veto entries - Find muons! =========
		//
		prev.Clear();
		MJVetoEvent prevLED;
		double xTimePrev = 0;
		double x_deltaTPrev = 0;
		double xTimePrevLED = 0;
		bool firstLED = false;
		bool IsLEDPrev = false;
		int almostMissedLED = 0;
		for (long i = 0; i < vEntries; i++) 
		// for (long i = 100; i < 110; i++) 
		{
			v->GetEntry(i);
			MJVetoEvent veto;
			veto.SetSWThresh(swThresh);	
	    	isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true);

	    	// find event time 
			if (!veto.GetBadScaler()) {
				xTime = veto.GetTimeSec();
			}
			else if (run > 8557 && veto.GetTimeSBC() < 2000000000) {
				xTime = veto.GetTimeSBC() - SBCOffset;
				double interpTime = InterpTime(i,LocalEntryTime,LocalEntryNum,LocalBadScalers);
				printf("Entry %li : SBC method: %.2f  Interp method: %.2f  sbc-interp: %.2f\n",i,xTime,interpTime,xTime-interpTime);
			}
			else {
				double eTime = ((double)i / vEntries) * d;
				xTime = InterpTime(i,LocalEntryTime,LocalEntryNum,LocalBadScalers);
				printf("Entry %li : Entry method: %.2f  Interp method: %.2f  eTime-interp: %.2f\n",i,eTime,xTime,eTime-xTime);
			}
			LocalEntryTime[i] = xTime;	// replace entry with the more accurate one
			
			// look at delta-t between events
			double dt = xTime - xTimePrev;
			if (dt > LEDperiod + RMSTimeWindow && i > 0){
				printf("High delta-T event: Entry %li, Prev %li.  dt = %.2f  xTime = %.2f  xTimePrev = %.2f  window: dt > %.2fs\n"
					,i,i-1,dt,xTime,xTimePrev,LEDperiod+RMSTimeWindow);
			}
			
			// skip bad entries (true = print contents of skipped event)
	    	if (CheckForBadErrors(veto,i,isGood,false)) continue;

			//----------------------------------------------------------	    	
	    	// 1: Energy Cut
	    	// The measured muon energy threshold is QDC = 500.  
	    	// Set TRUE if at least TWO panels are over 500.
	    	//
	    	bool EnergyCut = false;
	    	
	    	int over500Count = 0;
	    	for (int q = 0; q < 32; q++) {
	    		if (veto.GetQDC(q) > 500) 
	    			over500Count++;
	    	}
	    	if (over500Count >= 2) EnergyCut = true;

	    	//----------------------------------------------------------
			// 2: Time of event.
			// Employ alternate methods if the scaler is corrupted.
			// Should implement an estimate of the error when alternate methods are used.
			// 
			bool ApproxTime = false;

			xTime = -1; 

			if (!veto.GetBadScaler()) xTime = veto.GetTimeSec();
			else if (run > 8557 && veto.GetTimeSBC() < 2000000000) {
				xTime = veto.GetTimeSBC() - SBCOffset;
				ApproxTime = true;
			}
	    	else {
	    		xTime = ((double)i / vEntries) * d;
	    		ApproxTime = true;
	    	}

			//----------------------------------------------------------
			// 3. Time Cut / LED Tag.
			// 
			// TRUE if an event PASSES (i.e. is physics.)  FALSE if an event is an LED.
			//
			// If LED's are turned off or the frequency measurement is bad, we revert
			// to a simple multiplicity threshold.  
			//
			bool TimeCut = true;
			bool IsLED = false;

			// Set Cut
			x_deltaT = xTime - xTimePrevLED;
			if (!badLEDFreq && fabs(LEDperiod - x_deltaT) < RMSTimeWindow && veto.GetMultip() > highMultipThreshold) {
				TimeCut = false;
				IsLED = true;
			}
			
			// almost missed a high-multiplicity event somehow ...
			else if (!badLEDFreq && fabs(LEDperiod - x_deltaT) >= (LEDperiod - RMSTimeWindow) && veto.GetMultip() > highMultipThreshold)	{
				TimeCut = false;
				IsLED = true;	
				almostMissedLED++;
				cout << "Almost missed LED:\n";
	
				// check this entry
				printf("Current: %-3li  m %-3i LED? %i t %-6.2f LEDP %-5.2f  XDT %-6.2f LEDP-XDT %-6.2f RMSW %-6.2f\n"
					,i,veto.GetMultip(),IsLED,xTime,LEDperiod,x_deltaT,LEDperiod-x_deltaT,RMSTimeWindow);

				// check previous entry
				printf("Previous: %-3li  m %-3i LED? %i t %-6.2f LEDP %-5.2f  XDT %-6.2f LEDP-XDT %-6.2f RMSW %-6.2f\n"
					,i-1,prev.GetMultip(),IsLEDPrev,xTimePrev,LEDperiod,x_deltaTPrev,LEDperiod-x_deltaTPrev,RMSTimeWindow);
			}	
			else TimeCut = true;
			
			// Grab first LED
			if (!firstLED && veto.GetMultip() > highMultipThreshold) {
				printf("Found first LED.  i %-2li m %-2i t %-5.2f\n\n",i,veto.GetMultip(),xTime);
				IsLED=true;
				firstLED=true;
				TimeCut=false;
			}

			// If frequency measurement is bad, revert to standard multiplicity cut
			if (badLEDFreq && veto.GetMultip() >= highestMultip-5){
				IsLED = true;
				TimeCut = false;
			}

			// If LED is off, all events pass time cut.
			if (LEDTurnedOff) {
				IsLED = false;
				TimeCut = true;
			}

			// Check output
			// printf("%-3li  m %-3i LED? %i t %-6.2f LEDP %-5.2f  XDT %-6.2f LEDP-XDT %-6.2f RMSW %-6.2f\n"
				// ,i,veto.GetMultip(),IsLED,xTime,LEDperiod,x_deltaT,LEDperiod-x_deltaT,RMSTimeWindow);

			//----------------------------------------------------------
			// 4: Hit Pattern
			// Map hits above SW threshold to planes and count the hits.
			// 
			bool plane[12] = {false};
			int planeHits[12] = {0};
			for (int k = 0; k < 32; k++) 
			{
				if (veto.GetQDC(k) > veto.GetSWThresh(k))
				{
					if (PanelMap(k)==0) { plane[0]=true; planeHits[0]++; }			// 0: Lower Bottom
					else if (PanelMap(k)==1) { plane[1]=true; planeHits[1]++; }		// 1: Upper Bottom
					else if (PanelMap(k)==2) { plane[2]=true; planeHits[2]++; }		// 3: Inner Top
					else if (PanelMap(k)==3) { plane[3]=true; planeHits[3]++; }		// 4: Outer Top
					else if (PanelMap(k)==4) { plane[4]=true; planeHits[4]++; }		// 5: Inner North
					else if (PanelMap(k)==5) { plane[5]=true; planeHits[5]++; }		// 6: Outer North
					else if (PanelMap(k)==6) { plane[6]=true; planeHits[6]++; }		// 7: Inner South
					else if (PanelMap(k)==7) { plane[7]=true; planeHits[7]++; }		// 8: Outer South
					else if (PanelMap(k)==8) { plane[8]=true; planeHits[8]++; }		// 9: Inner West
					else if (PanelMap(k)==9) { plane[9]=true; planeHits[9]++; }		// 10: Outer West
					else if (PanelMap(k)==10) { plane[10]=true; planeHits[10]++; }	// 11: Inner East
					else if (PanelMap(k)==11) { plane[11]=true; planeHits[11]++; }	// 12: Outer East
				}
			}

			//----------------------------------------------------------
			// 5: Muon Identification
			// Use EnergyCut, TimeCut, and the Hit Pattern to identify them sumbitches.

			// reset
			for (int r = 0; r < 32; r++) {CoinType[r]=0; CutType[r]=0;}

			if (TimeCut && EnergyCut)
			{
				// 0. Everything that passes TimeCut and EnergyCut.
				CoinType[0] = true;

				// 1. Definite Vertical Muons
				if (plane[0] && plane[1] && plane[2] && plane[3]) {
					CoinType[1] = true;
					printf("Found Vertical Muon. QDC:%i  Mult:%i\n",veto.GetTotE(),veto.GetMultip());
					printf("%-3li  m %-3i LED? %i t %-6.2f  XDT %-6.2f  LEDP-XDT %-6.2f  RMSW %-6.2f\n\n"
						,i,veto.GetMultip(),IsLED,xTime,x_deltaT,LEDperiod-x_deltaT,RMSTimeWindow);
				}

				// 2. Both top or side layers + both bottom layers.
				if ((plane[0] && plane[1]) && ((plane[2] && plane[3]) || (plane[4] && plane[5])
					|| (plane[6] && plane[7]) || (plane[8] && plane[9]) || (plane[10] && plane[11]))) {
					CoinType[2] = true;
					
					// show output if we haven't seen it from CT1 already
					if (!CoinType[1]) { 
						printf("Found Vert/Side Muon. QDC:%i  Mult:%i\n",veto.GetTotE(),veto.GetMultip());
						printf("%-3li  m %-3i LED? %i t %-6.2f  XDT %-6.2f  LEDP-XDT %-6.2f  RMSW %-6.2f\n\n"
							,i,veto.GetMultip(),IsLED,xTime,x_deltaT,LEDperiod-x_deltaT,RMSTimeWindow);
					}
				}

				// Other coincidence types can be found by parsing the ROOT output.
			}

			//----------------------------------------------------------
			// 6: Output
			// The skim file used to take a text file of muon candidate events.
			// Additionally, write the ROOT file containing all the real data.
			// 

			// Write a text file
			if (list) {
				char buffer[200];
				if (CoinType[1] || CoinType[0]) {
					int type;
					if (CoinType[0]) type = 1;
					if (CoinType[1]) type = 2;
					sprintf(buffer,"%i %li %.8f %i %i\n",run,start,xTime,type,veto.GetBadScaler());
					MuonList << buffer;
				}
				// This is Jason's TYPE 3: flag runs with gaps since the last stop time.
				if ((start - prevStopTime) > 10 && i == 0) {
					sprintf(buffer,"%i %li 0.0 3 0\n",run,start);
					MuonList << buffer;
				}
			}

			// Assign all bools calculated to the int array CutType[32];
			CutType[0] = LEDTurnedOff;
			CutType[1] = EnergyCut;
			CutType[2] = ApproxTime;
			CutType[3] = TimeCut;
			CutType[4] = IsLED;
			CutType[5] = firstLED;
			CutType[6] = badLEDFreq;
			for (int c = 0; c < 12; c++) CutType[7+c] = plane[c];

			// Write ROOT output
			if (root) {
				out = veto;
				vetoEvent->Fill();
			}

			//----------------------------------------------------------
			if (IsLED) {
				prevLED = veto;
				xTimePrevLED = xTime;
			}
			IsLEDPrev = IsLED;
			prev = veto;
			xTimePrev = xTime;
			x_deltaTPrev = x_deltaT;
	    } 	

	    // End of run summaries.
		if (almostMissedLED > 0) cout << "\nWarning, almost missed " << almostMissedLED << " LED events.\n";

	    // done with this run.
		delete ds;
		prevStopTime = stop;
	}

	printf("\n===================== End of Scan. =====================\n");


	if (list) MuonList.close();
	if (root) vetoEvent->Write();
	RootFile->Close();
}