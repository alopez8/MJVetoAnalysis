#include "vetoScan.hh"

using namespace std;

void vetoPerformance(string Input, int *thresh, bool runBreakdowns) 
{
	// input a list of run numbers
	ifstream InputList(Input.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << Input << endl;
    	return;
    }
    int filesScanned = 0;	// 1-indexed.

    // output a ROOT file
	string Name = Input;
	Name.erase(Name.find_last_of("."),string::npos);
	Name.erase(0,Name.find_last_of("\\/")+1);
	Char_t OutputFile[200];
	sprintf(OutputFile,"./output/vPerf_%s.root",Name.c_str());
	TFile *RootFile = new TFile(OutputFile, "RECREATE"); 	
  	TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."
	RootFile->mkdir("rawQDC");
	if (runBreakdowns) RootFile->mkdir("runPlots");

	// global counters
	const int nErrs = 18;
	int globalErrorCount[nErrs] = {0};
	int globalRunsWithErrors[nErrs] = {0};
	int globalRunsWithErrorsAtBeginning[nErrs] = {0};
	int globalErrorAtBeginningCount[nErrs] = {0};
	int SJCount = 0;
	vector<double> runs;
	vector<double> freqs;
	vector<double> ErrCountEntry;
	vector<double> EntryTime;
	vector<double> EntryNum;
	vector<int> HighDTEvent;
	vector<double> SJTime;
	vector<int> SJIndex;
	long totEntries = 0;
	long totDuration = 0;
	int totHighDT = 0;
	int totHighDTwBTS = 0;
	
	// global histograms and graphs
	TGraph *gRunVsLEDFreq;				// depends on: runs & freqs
	TGraph *gErrorCountEntryVsTime; 	// depends on: ErrorCountEntry & EntryTime
	TGraph *gErrorCountEntryVsEntryNum; // depends on: ErrorCountEntry & EntryNum

	TH1D *TotalMultip = new TH1D("TotalMultip","Events over threshold",33,0,33);
	TotalMultip->GetXaxis()->SetTitle("number of panels hit");
	
	TH1D *TotalEnergy = new TH1D("TotalEnergy","Total QDC from events",100,0,60000);
	TotalEnergy->GetXaxis()->SetTitle("energy (QDC)");

	TH1D *deltaT = new TH1D("deltaT","Time between successive entries",200,0,20);
	deltaT->GetXaxis()->SetTitle("seconds");
	
	TH1D *TotalEnergyNoLED = new TH1D("TotalEnergyNoLED","Total QDC from non-LED events",100,0,60000);
	TotalEnergyNoLED->GetXaxis()->SetTitle("energy (QDC)");

	TH1D *QDC_over_Multip = new TH1D("QDC_over_Multip","Average QDC from events",1000,0,5000);
	QDC_over_Multip->GetXaxis()->SetTitle("Average energy (QDC)");
	
	TH1D *TimestampBadEntry = new TH1D("TimestampBadEntry"," Timestamp of entries with > 2 errors",3650,0,3650);
	TimestampBadEntry->GetXaxis()->SetTitle("seconds");
	
	TH1D *hRawQDC[32];
	char hname[50];
	for (int i=0; i<32; i++)
	{
		sprintf(hname,"hRawQDC%d",i);
		hRawQDC[i] = new TH1D(hname,hname,4200,0,4200);
	}
	
	// ==========================loop over input files==========================
	//
	while(!InputList.eof())
	{
		int run = 0;
		InputList >> run;
		filesScanned++;
		
		// initialize
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
		long start = (long)vRun->GetStartTime();
		long stop = (long)vRun->GetStopTime();
		double duration = (double)(stop - start);
		totEntries += vEntries;
		totDuration += (long)duration;

		// run-by-run variables
		int errorCount[nErrs] = {0};
		vector<double> LocalErrCountEntry;
		vector<double> LocalEntryTime;
		vector<double> LocalEntryNum;
		vector<bool> LocalBadScalers;	

		// run-by-run histos and graphs
		sprintf(hname,"%d_LEDDeltaT",run);
		TH1D *LEDDeltaT = new TH1D(hname,hname,100000,0,100); // 0.001 sec/bin
		TH1D *deltaTRun = NULL;
		TGraph *gMultipVsTimeRun = NULL;
		TGraph *gSTimeVsfIndex = NULL;	
		TGraph *gEventCountScaler = NULL;
		TGraph *gEventCountQDC1 = NULL;
		TGraph *gEventCountQDC2 = NULL;
		if (runBreakdowns)
		{
			sprintf(hname,"%d_deltaT", run);
			deltaTRun = new TH1D(hname,hname,200,0,20);

			sprintf(hname,"%d_MultipVsTime", run);
			gMultipVsTimeRun = new TGraph(vEntries);
			gMultipVsTimeRun->SetName(hname);
			
			sprintf(hname,"%d_STimeVsfIndex", run);
			gSTimeVsfIndex = new TGraph(vEntries);
			gSTimeVsfIndex->SetName(hname);
			
			sprintf(hname,"%d_EventCountScaler", run);
			gEventCountScaler = new TGraph(vEntries);
			gEventCountScaler->SetName(hname);
			
			sprintf(hname,"%d_EventCountQDC1", run);
			gEventCountQDC1 = new TGraph(vEntries);
			gEventCountQDC1->SetName(hname);
			
			sprintf(hname,"%d_EventCountQDC2", run);
			gEventCountQDC2 = new TGraph(vEntries);
			gEventCountQDC2->SetName(hname);
		}

		printf("\n======= Scanning run %i, %li entries, %.0f sec. =======\n",run,vEntries,duration);
		MJVetoEvent prev;
		MJVetoEvent first;
		bool foundFirst = false;
		int firstGoodEntry = 0;
		int pureLEDcount = 0;
		bool errorRunBools[nErrs] = {0};
		bool errorRunBeginningBools[nErrs] = {0};
		int highestMultip = 0;
		double xTime = 0;
		double lastGoodTime = 0;
		bool FirstHighMultip = false;

		// ====================== First loop over entries =========================
		for (int i = 0; i < vEntries; i++)
		{
			v->GetEntry(i);
			MJVetoEvent veto;
			veto.SetSWThresh(thresh);	
	    	int isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true); // true: force-write event with errors.
	    	
	    	// count up error types
	    	int errorsThisEntry = 0; 
	    	if (isGood != 1) 
	    	{	    		
	    		for (int j=0; j<nErrs; j++) if (veto.GetError(j)==1) 
	    		{
	    			errorCount[j]++;
	    			errorsThisEntry++;
	    			errorRunBools[j]=true;
	    			if (i < 10) {
	    				errorRunBeginningBools[j]=true;
	    				globalErrorAtBeginningCount[j]++;
	    			}
	    		}
	    	}
				
	    	// find event time and fill vectors
			if (!veto.GetBadScaler()) {
				LocalBadScalers.push_back(0);
				xTime = veto.GetTimeSec();
			}
			else {
				LocalBadScalers.push_back(1);
				xTime = ((double)i / vEntries) * duration;
			}
			
	    	// fill vectors
	    	// (the time vectors are revised in the second loop)
			EntryNum.push_back(i);
			EntryTime.push_back(xTime);
			ErrCountEntry.push_back(errorsThisEntry);
			LocalEntryNum.push_back(i);		
			LocalEntryTime.push_back(xTime);
			LocalErrCountEntry.push_back(errorsThisEntry);
			
			// skip bad entries (true = print contents of skipped event)
	    	if (CheckForBadErrors(veto,i,isGood,false)) continue;

    		// Save the first good entry number for the SBC offset
			if (isGood == 1 && !foundFirst) {
				first = veto;
				foundFirst = true;
				firstGoodEntry = i;
			}
			
			// find the highest multiplicity in this run (used in 2nd loop)
	    	if (veto.GetMultip() > highestMultip && veto.GetMultip() < 33) {
	    		highestMultip = veto.GetMultip();
	    		cout << "Finding highest multiplicity: " << highestMultip << "  entry: " << i << endl;
	    	}

	    	// very simple LED tag 
			if (veto.GetMultip() >= 20) {
				LEDDeltaT->Fill(veto.GetTimeSec()-prev.GetTimeSec());
				pureLEDcount++;
			}
			
			// end of loop : save things
			prev = veto;
			lastGoodTime = xTime;
			
		}

		// Make sure the local vectors are all the same size
		if ((LocalEntryNum.size() != LocalEntryTime.size()) || (LocalEntryNum.size() != LocalErrCountEntry.size()))
		printf("Warning! Local vectors are not the same size!\n");

		// if duration is corrupted, use the last good timestamp as the duration.
		if (duration == 0) {
			printf("Corrupted duration. Using last good timestamp: %.2f\n",lastGoodTime);
			duration = lastGoodTime;
			totDuration += duration;
		}

		// find the SBC offset		
		double SBCOffset = first.GetTimeSBC() - first.GetTimeSec();
		printf("First good entry: %i  SBCOffset: %.2f\n",firstGoodEntry,SBCOffset);

		// find the LED frequency, set time window, 
		double RMSTimeWindow = 0.1;
		printf("\"Simple\" LED count: %i.  Approx rate: %.3f\n",pureLEDcount,pureLEDcount/duration);
		double LEDrms = 0;
		double LEDfreq = 0;
		int dtEntries = LEDDeltaT->GetEntries();
		if (dtEntries > 0) {
			int maxbin = LEDDeltaT->GetMaximumBin();
			LEDDeltaT->GetXaxis()->SetRange(maxbin-100,maxbin+100); // looks at +/- 0.1 seconds of max bin.
			LEDrms = LEDDeltaT->GetRMS();
			LEDfreq = 1/LEDDeltaT->GetMean();
		}
		else {
			printf("Warning! No multiplicity > 20 events!!\n");
			LEDrms = 9999;
			LEDfreq = 9999;
		}
		double LEDperiod = 1/LEDfreq;
		printf("Histo method: LED_f: %.8f LED_t: %.8f RMS: %8f\n",LEDfreq,LEDperiod,LEDrms);
		delete LEDDeltaT;
		if (LEDfreq != 9999 && vEntries > 100) {
			runs.push_back(run);
			freqs.push_back(LEDfreq);
		}

		// set a flag for "bad LED" (usually a short run causes it)
		// and replace the period with the "simple" one if possible
		bool badLEDFreq = false;
		if (LEDperiod > 9 || vEntries < 100) 
		{
			printf("Warning: Short run.\n");
			if (pureLEDcount > 3) {
				printf("   From histo method, LED freq is %.2f.\n   Reverting to the approx rate (%.2fs) ... \n"
					,LEDfreq,(double)pureLEDcount/duration);
				LEDperiod = duration/pureLEDcount;
			}
			else { 
				printf("   Warning: LED info is corrupted!  Will not use LED period information for this run.\n");
				LEDperiod = 9999;
				badLEDFreq = true;
			}
		}

		// add error counts to global totals
		for (int q = 0; q < nErrs; q++) {
			if (errorRunBools[q]) {
				globalRunsWithErrors[q]++;
				if (q == 1) printf("Missing Channels in run %d\n",run);
				if (q == 6) printf("Duplicate Channels in run %d\n",run);
				if (q == 7) printf("Hardware Count Mismatch in run %d\n",run);
			}	
			if (errorRunBeginningBools[q]) globalRunsWithErrorsAtBeginning[q]++;
		}
		
		// ====================== Second loop over entries =========================
		//
		double xTimePrev = 0;
		int TimeMethod = 0; //1 = scaler, 2 = SBC, 3 = interp
		int TimeMethodPrev = 0;
		double STime = 0;
		double STimePrev = 0;
		int SIndex = 0;
		int SIndexPrev = 0;
		double SSlope = 0;
		double SSlopePrev = 0;
		
		for (int i = 0; i < vEntries; i++)
		{
			// this time we don't skip anything until all the time information is found.
			v->GetEntry(i);
			MJVetoEvent veto;
			veto.SetSWThresh(thresh);	
	    	int isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true);	// true: force-write event with errors.
			bool isSJ = false;
			
	    	// find event time 
			if (!veto.GetBadScaler()) {
				xTime = veto.GetTimeSec();
				STime = veto.GetTimeSec();
				SIndex = veto.GetScalerIndex();
				SSlope = (STime-STimePrev)/((double)SIndex-(double)SIndexPrev);
				TimeMethod = 1;
			}
			else if (run > 8557 && veto.GetTimeSBC() < 2000000000) {
				xTime = veto.GetTimeSBC() - SBCOffset;
				double interpTime = InterpTime(i,LocalEntryTime,LocalEntryNum,LocalBadScalers);
				printf("Entry %i : SBC method: %.2f  Interp method: %.2f  sbc-interp: %.2f\n",i,xTime,interpTime,xTime-interpTime);
				TimeMethod = 2;
			}
			else {
				double eTime = ((double)i / vEntries) * duration;
				xTime = InterpTime(i,LocalEntryTime,LocalEntryNum,LocalBadScalers);
				printf("Entry %i : Entry method: %.2f  Interp method: %.2f  eTime-interp: %.2f\n",i,eTime,xTime,eTime-xTime);
				TimeMethod = 3;
			}
			LocalEntryTime[i] = xTime;	// replace entry with the more accurate one
			
			// look at delta-t between events
			double dt = xTime - xTimePrev;
			deltaT->Fill(dt);
			if (runBreakdowns) { 
				deltaTRun->Fill(dt);
				gMultipVsTimeRun->SetPoint(i,LocalEntryTime[i],veto.GetMultip());
				if (!veto.GetBadScaler()) {
					gSTimeVsfIndex->SetPoint(i,veto.GetScalerIndex(),veto.GetTimeSec());		
				}
				gEventCountScaler->SetPoint(i,xTime,veto.GetSEC());
				gEventCountQDC1->SetPoint(i,xTime,veto.GetQEC());
				gEventCountQDC2->SetPoint(i,xTime,veto.GetQEC2());
			}
			if (dt > LEDperiod + RMSTimeWindow && i > 0){
				//printf("High delta-T event: Entry %i, Prev %i.  dt = %.2f  xTime = %.2f (Method: %d) xTimePrev = %.2f (MethodPrev: %d)  window: dt > %.2fs\n"
				//	,i,i-1,dt,xTime,TimeMethod,xTimePrev,TimeMethodPrev,trueLEDperiod+RMSTimeWindow);
				HighDTEvent.push_back(i-1);
				HighDTEvent.push_back(i);
				totHighDT++;
				if (LocalBadScalers[i-1] == 1 || LocalBadScalers[i] == 1) totHighDTwBTS++;
			}
			
			//look for scaler jumps
			if (SIndex != SIndexPrev && SIndexPrev != 0){
				if( ((STime-STimePrev) > 25 && fabs((SSlope-SSlopePrev)/SSlopePrev) > 0.8) || (STime-STimePrev) > 25 && fabs(SSlope) > 0.19) { //by examining known Scaler Jumps
					SJCount++;
					isSJ = true;
					printf("Scaler Jump found!!! Run: %d  |  Entry: %d  |  ScalerTime: %f  |  ScalerTimePrev: %f  |  ScalerIndex: %d  |  ScalerIndexPrev: %d\n",run,i,STime,STimePrev,SIndex,SIndexPrev); 
				
				}
			}	
			
			// save previous xTime
			xTimePrev = xTime;
			STimePrev = STime;
			SIndexPrev = SIndex;
			if (!isSJ) SSlopePrev = SSlope;
			TimeMethodPrev = TimeMethod;
			
			// skip bad entries (true = print contents of skipped event)
	    	if (CheckForBadErrors(veto,i,isGood,false)) continue;

			// fill energy/multiplicity histos
	    	TotalEnergy->Fill(veto.GetTotE());
	    	TotalMultip->Fill(veto.GetMultip());
			QDC_over_Multip->Fill(veto.GetTotE()/(double)veto.GetMultip());	    	
	    	
	    	for (int j = 0; j < 32; j++) 
	    		hRawQDC[j]->Fill(veto.GetQDC(j));
			
			if (veto.GetMultip() <= 20) 
				TotalEnergyNoLED->Fill(veto.GetTotE());
			
			if (veto.GetMultip() < highestMultip && veto.GetMultip() > 8){
				
				printf("Found event with multiplicity > 8 and < highestMultip ... Multip: %i  Entry: %i\n",veto.GetMultip(),i);
				if (!FirstHighMultip){
					printf("First Strange Multip Event: Entry %d\n",(int)LocalEntryNum[i]);
					veto.Print();
					FirstHighMultip =  true;
				}	
				
			}
			
		}		
		
		cout << "=================== End Run " << run << ". =====================\n";
		for (int i = 0; i < nErrs; i++) {
			if (errorCount[i] > 0) {
				printf("%i: %i errors\t(%.2f%% of total)\n",i,errorCount[i],100*(double)errorCount[i]/vEntries);
				globalErrorCount[i] += errorCount[i];
			}
		}
		// end of run cleanup
		LocalBadScalers.clear();
		LocalEntryNum.clear();
		LocalEntryTime.clear();
		LocalErrCountEntry.clear();
		HighDTEvent.clear();
		if (runBreakdowns) 
		{
			RootFile->cd("runPlots");
			sprintf(hname,"%d_deltaT", run);
			deltaTRun->Write(hname,TObject::kOverwrite); 

			sprintf(hname,"%d_MultipVsTime", run);
			gMultipVsTimeRun->SetMarkerColor(4);
			gMultipVsTimeRun->SetMarkerStyle(21);
			gMultipVsTimeRun->SetMarkerSize(0.5);
			gMultipVsTimeRun->SetLineColorAlpha(kWhite,0);
			gMultipVsTimeRun->Write(hname,TObject::kOverwrite);
			
			sprintf(hname,"%d_STimeVsfIndex", run);
			gSTimeVsfIndex->GetXaxis()->SetTitle("Scaler Index");
			gSTimeVsfIndex->GetYaxis()->SetTitle("Scaler Time (sec)");
			gSTimeVsfIndex->SetMarkerColor(4);
			gSTimeVsfIndex->SetMarkerStyle(21);
			gSTimeVsfIndex->SetMarkerSize(0.5);
			gSTimeVsfIndex->SetLineColorAlpha(kWhite,0);
			gSTimeVsfIndex->Write(hname,TObject::kOverwrite);
			
			sprintf(hname,"%d_EventCountScaler", run);
			gEventCountScaler->SetMarkerStyle(20);
			gEventCountScaler->SetMarkerColor(2);
			gEventCountScaler->SetLineColorAlpha(kWhite,0);
			gEventCountScaler->Write(hname,TObject::kOverwrite);
			
			sprintf(hname,"%d_EventCountQDC1", run);
			gEventCountQDC1->SetMarkerStyle(21);
			gEventCountQDC1->SetMarkerColor(4);
			gEventCountQDC1->SetLineColorAlpha(kWhite,0);
			gEventCountQDC1->Write(hname,TObject::kOverwrite);
			
			sprintf(hname,"%d_EventCountQDC2", run);
			gEventCountQDC2->SetMarkerStyle(22);
			gEventCountQDC2->SetMarkerColor(6);
			gEventCountQDC2->SetLineColorAlpha(kWhite,0);
			gEventCountQDC2->Write(hname,TObject::kOverwrite);	

			delete deltaTRun;
			delete gMultipVsTimeRun;
			delete gSTimeVsfIndex;
			delete gEventCountScaler;
			delete gEventCountQDC1;
			delete gEventCountQDC2;
			RootFile->cd();
		}	
	}
	
	for (int i = 0; i < (int)ErrCountEntry.size(); i++){
		if (ErrCountEntry[i] > 2) TimestampBadEntry->Fill(EntryTime[i]);
	}
	
	cout << "\n\n================= END OF SCAN. =====================\n";
	printf("%i runs, %li total events, total duration: %ld seconds.\n",filesScanned,totEntries,totDuration);

	// check for errors and print summary if we find them.
	bool foundErrors = false;
	for (int i = 0; i < nErrs; i++) { 
		if (globalErrorCount[i] > 0) {
			foundErrors = true; 
			break;
		}
	}

	if (foundErrors) 
	{
		printf("\nError summary:\n");
		for (int i = 0; i < nErrs; i++) 
		{
			if (globalErrorCount[i] > 0) 
			{
				foundErrors = true;
				printf("%i: %i events\t(%.2f%%)\t"
					,i,globalErrorCount[i],100*(double)globalErrorCount[i]/totEntries);
				printf("%i runs\t(%.2f%%)\n"
					,globalRunsWithErrors[i],100*(double)globalRunsWithErrors[i]/filesScanned);
			}
		}
		if (totHighDT>0) printf("High-DeltaT Events: %i  High DT Events with BadScaler: %i\n\n",totHighDT,totHighDTwBTS);
		printf("Number of Scaler Jumps: %d",SJCount);
		
		printf("\nBeginning of runs (i < 10) error summary:\n");
		for (int i = 0; i < nErrs; i++) 
		{
			if (globalErrorAtBeginningCount[i] > 0) 
			{
				printf("%i: %i events\t (%.2f%%)\t",
					i,globalErrorAtBeginningCount[i],100*(double)globalErrorAtBeginningCount[i]/totEntries);
				printf("%i runs\t(%.2f%%)\n",
					globalRunsWithErrorsAtBeginning[i],100*(double)globalRunsWithErrorsAtBeginning[i]/filesScanned);
			}
		}
		printf("\nFor reference, error types are:\n");
		cout << "1. Missing channels (< 32 veto datas in event) " << endl;
		cout << "2. Extra Channels (> 32 veto datas in event) " << endl; 
		cout << "3. Scaler only (no QDC data) " << endl;
		cout << "4. Bad Timestamp: FFFF FFFF FFFF FFFF " << endl;
		cout << "5. QDCIndex - ScalerIndex != 1 or 2 " << endl;
		cout << "6. Duplicate channels (channel shows up multiple times) " << endl;
		cout << "7. HW Count Mismatch (SEC - QEC != 1 or 2) " << endl;
		cout << "8. MJTRun run number doesn't match input file" << endl;
		cout << "9. MJTVetoData cast failed (missing QDC data)" << endl;
		cout << "10. Scaler EventCount doesn't match ROOT entry" << endl;
		cout << "11. Scaler EventCount doesn't match QDC1 EventCount" << endl;
		cout << "12. QDC1 EventCount doesn't match QDC2 EventCount" << endl;
		cout << "13. Indexes of QDC1 and Scaler differ by more than 2" << endl;
		cout << "14. Indexes of QDC2 and Scaler differ by more than 2" << endl;
		cout << "15. Indexes of either QDC1 or QDC2 PRECEDE the scaler index" << endl;
		cout << "16. Indexes of either QDC1 or QDC2 EQUAL the scaler index" << endl;
		cout << "17. Unknown Card is present." << endl;
	}
	
	// write global plots
	gRunVsLEDFreq = new TGraph(runs.size(),&(runs[0]),&(freqs[0]));
	gRunVsLEDFreq->SetTitle("LED Frequency vs Run Number");
	gRunVsLEDFreq->GetXaxis()->SetTitle("Run Number");
	gRunVsLEDFreq->GetYaxis()->SetTitle("LED Freq (Hz)");
	gRunVsLEDFreq->SetMarkerColor(4);
	gRunVsLEDFreq->SetMarkerStyle(21);
	gRunVsLEDFreq->SetMarkerSize(0.5);
	gRunVsLEDFreq->SetLineColorAlpha(kWhite,0);
	gRunVsLEDFreq->Write("RunVsLEDFreq",TObject::kOverwrite);
	
	gErrorCountEntryVsTime = new TGraph(EntryTime.size(),&(EntryTime[0]),&(ErrCountEntry[0]));
	gErrorCountEntryVsTime->SetTitle("Error Count Vs Entry Time");
	gErrorCountEntryVsTime->GetXaxis()->SetTitle("Entry Time (sec)");
	gErrorCountEntryVsTime->GetYaxis()->SetTitle("Error Count");
	gErrorCountEntryVsTime->SetMarkerColor(4);
	gErrorCountEntryVsTime->SetMarkerStyle(21);
	gErrorCountEntryVsTime->SetMarkerSize(0.5);
	gErrorCountEntryVsTime->SetLineColorAlpha(kWhite,0);
	gErrorCountEntryVsTime->Write("ErrorCountEntryVsTime",TObject::kOverwrite);
	
	gErrorCountEntryVsEntryNum = new TGraph(EntryNum.size(),&(EntryNum[0]),&(ErrCountEntry[0]));
	gErrorCountEntryVsEntryNum->SetTitle("Error Count vs Entry Number");
	gErrorCountEntryVsEntryNum->GetXaxis()->SetTitle("Entry Number");
	gErrorCountEntryVsEntryNum->GetYaxis()->SetTitle("Error Count");
	gErrorCountEntryVsEntryNum->SetMarkerColor(4);
	gErrorCountEntryVsEntryNum->SetMarkerStyle(21);
	gErrorCountEntryVsEntryNum->SetMarkerSize(0.5);
	gErrorCountEntryVsEntryNum->SetLineColorAlpha(kWhite,0);
	gErrorCountEntryVsEntryNum->Write("ErrorCountEntryVsEntryNum",TObject::kOverwrite);

	TotalMultip->Write("TotalMultip",TObject::kOverwrite);
	TotalEnergy->Write("TotalEnergy",TObject::kOverwrite);
	TotalEnergyNoLED->Write("TotalEnergyNoLED",TObject::kOverwrite);
	QDC_over_Multip->Write("QDC_over_Multip",TObject::kOverwrite);
	TimestampBadEntry->Write("TimestampBadEntry",TObject::kOverwrite);

	
	deltaT->Write("deltaT",TObject::kOverwrite);

	RootFile->cd("rawQDC");
	for (int i=0;i<32;i++)
	{	
		sprintf(hname,"hRawQDC%d",i);
		hRawQDC[i]->Write(hname,TObject::kOverwrite);
	}
	
	RootFile->Close();
	cout << "\nWrote ROOT file." << endl;
}