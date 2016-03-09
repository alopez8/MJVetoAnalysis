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
	if (runBreakdowns) RootFile->mkdir("runPlots");

	// global counters
	const int nErrs = 18;
	int globalErrorCount[nErrs] = {0};
	int globalRunsWithErrors[nErrs] = {0};
	int globalRunsWithErrorsAtBeginning[nErrs] = {0};
	int globalErrorAtBeginningCount[nErrs] = {0};

	// global histograms and graphs
	TH1D *TotalMultip = new TH1D("TotalMultip","Events over threshold",33,0,33);
	TotalMultip->GetXaxis()->SetTitle("number of panels hit");
	TH1D *TotalEnergy = new TH1D("TotalEnergy","Total QDC from events",100,0,60000);
	TotalEnergy->GetXaxis()->SetTitle("energy (QDC)");
	TH1D *TotalEnergyNoLED = new TH1D("TotalEnergyNoLED","Total QDC from non-LED events",100,0,60000);
	TotalEnergyNoLED->GetXaxis()->SetTitle("energy (QDC)");
	TH1F *hRawQDC[32];
	char hname[50];
	for (int i=0; i<32; i++){
		sprintf(hname,"hRawQDC%d",i);
		hRawQDC[i] = new TH1F(hname,hname,4200,0,4200);
	}
	TGraph *gRunVsLEDFreq;	// initialized at end of scan
	vector<double> runs;
	vector<double> freqs;
	TH1D *deltaT = new TH1D("deltaT","Time between successive entries",200,0,20);
	deltaT->GetXaxis()->SetTitle("seconds");
	
	// ============================================================================================
	// loop over input file
	int run = 0;
	long totEntries = 0;
	long totDuration = 0;
	while(!InputList.eof())
	{
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
		double duration = ds->GetRunTime()/CLHEP::second;	
		if (duration == 0) duration = (double)(stop - start);
		totEntries += vEntries;
		totDuration += (long)duration;

		// run-by-run variables
		int errorCount[nErrs] = {0};

		// run-by-run histos and graphs
		sprintf(hname,"%d_LEDDeltaT",run);
		TH1D *LEDDeltaT = new TH1D(hname,hname,100000,0,100); // 0.001 sec/bin
		
		TH1D *deltaTRun = NULL;

		TGraph *gMultipVsTimeRun = NULL;
		if (runBreakdowns)
		{
			sprintf(hname,"%d_deltaT", run);
			deltaTRun = new TH1D(hname,hname,200,0,20);

			sprintf(hname,"%d_MultipVsTime", run);
			gMultipVsTimeRun = new TGraph(vEntries);
			gMultipVsTimeRun->SetName(hname);
		}

		printf("========== Scanning run %i, %li entries, Duration: %f.=============\n",run,vEntries,duration);
		MJVetoEvent prev;
		MJVetoEvent first;
		bool foundFirst = false;
		int firstGoodEntry = 0;
		int pureLEDcount = 0;
		bool errorRunBools[nErrs] = {0};
		bool errorRunBeginningBools[nErrs] = {0};
		// ====================== First loop over entries =========================
		//
		for (int i = 0; i < vEntries; i++)
		{
			v->GetEntry(i);
			MJVetoEvent veto;
			veto.SetSWThresh(thresh);	
	    	int isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true);
	    	if (isGood != 1) 
	    	{	    		
	    		for (int j=0; j<nErrs; j++) if (veto.GetError(j)==1) {
	    			errorCount[j]++;
	    			errorRunBools[j]=true;
	    			if (i < 10) {
	    				errorRunBeginningBools[j]=true;
	    				globalErrorAtBeginningCount[j]++;
	    			}
	    		}
	    	}
	    	if (CheckForBadErrors(veto,i,isGood,true)) continue;

    		// Save the first good entry number for the SBC offset
			if (isGood == 1 && !foundFirst) {
				first = veto;
				foundFirst = true;
				firstGoodEntry = i;
			}

			// fill QDC histos
	    	for (int j = 0; j < 32; j++) hRawQDC[j]->Fill(veto.GetQDC(j));
	    	TotalEnergy->Fill(veto.GetTotE());
			if (veto.GetMultip() < 20) TotalEnergyNoLED->Fill(veto.GetTotE());

	    	// fill multiplicity histo
	    	TotalMultip->Fill(veto.GetMultip());
			
	    	// very simple LED tag
			if (veto.GetMultip() >= 20) {
				LEDDeltaT->Fill(veto.GetTimeSec()-prev.GetTimeSec());
				pureLEDcount++;
			}

			// save this entry
			prev = veto;
		}

		// Find the SBC offset		
		double SBCOffset = first.GetTimeSBC() - first.GetTimeSec();
		printf("First good entry: %i  SBCOffset: %.2f\n",firstGoodEntry,SBCOffset);

		// Find the LED frequency 
		printf("\"Pure\" LED count: %i.  Approx rate: %.3f\n",pureLEDcount,pureLEDcount/duration);
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
		printf("LED_f: %.8f LED_t: %.8f RMS: %8f\n",LEDfreq,LEDperiod,LEDrms);
		delete LEDDeltaT;
		if (LEDfreq != 9999 && vEntries > 100) {
			runs.push_back(run);
			freqs.push_back(LEDfreq);
		}

		// Count runs with errors
		for (int q = 0; q < nErrs; q++) {
			if (errorRunBools[q]) globalRunsWithErrors[q]++;
			if (errorRunBeginningBools[q]) globalRunsWithErrorsAtBeginning[q]++;
		}

		// ====================== Second loop over entries =========================
		//
		double xTime = 0;
		double xTimePrev = 0;
		for (int i = 0; i < vEntries; i++)
		{
			v->GetEntry(i);
			MJVetoEvent veto;
			veto.SetSWThresh(thresh);	
	    	int isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true);
	    	if (CheckForBadErrors(veto,i,isGood,false)) continue;

	    	if (!veto.GetBadScaler()) xTime = veto.GetTimeSec();
	    	else if (run > 8557) xTime = veto.GetTimeSBC() - SBCOffset;
	    	else xTime = ((double)i / vEntries) * duration;

	    	// fill time histograms.
	    	deltaT->Fill(xTime-xTimePrev);
	    	if (runBreakdowns) { 
	    		deltaTRun->Fill(xTime-xTimePrev);
		    	gMultipVsTimeRun->SetPoint(i,xTime,veto.GetMultip());
		    }

	    	xTimePrev = xTime;
		}
		cout << "================ End Run " << run << ". ===================\n";
		for (int i = 0; i < nErrs; i++) {
			if (errorCount[i] > 0) {
				printf("%i: %i errors\t(%.2f%% of total)\n",i,errorCount[i],100*(double)errorCount[i]/vEntries);
				globalErrorCount[i] += errorCount[i];
			}
		}

		// write run-by-run plots to file and delete
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
			delete deltaTRun;

			delete gMultipVsTimeRun;
			RootFile->cd();
		}
	}

	cout << "\n\n============== END OF SCAN. ==================\n";
	
	// Error count summaries.
	printf("Error summary: %li events total, %i runs.\n",totEntries,filesScanned);

	for (int i = 0; i < nErrs; i++) {
		if (globalErrorCount[i] > 0) 
		{
			printf("%i: %i events\t(%.2f%%)\t"
				,i,globalErrorCount[i],100*(double)globalErrorCount[i]/totEntries);
			printf("%i runs\t(%.2f%%)\n"
				,globalRunsWithErrors[i],100*(double)globalRunsWithErrors[i]/filesScanned);
		}
	}
	printf("\nBeginning of runs (i < 10) error summary: %i runs total.\n",filesScanned);
	for (int i = 0; i < nErrs; i++) {
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

	// write global plots
	gRunVsLEDFreq = new TGraph(runs.size(),&(runs[0]),&(freqs[0]));
	gRunVsLEDFreq->Write("RunVsLEDFreq",TObject::kOverwrite);

	TotalMultip->Write("TotalMultip",TObject::kOverwrite);
	TotalEnergy->Write("TotalEnergy",TObject::kOverwrite);
	TotalEnergyNoLED->Write("TotalEnergyNoLED",TObject::kOverwrite);
	deltaT->Write("deltaT",TObject::kOverwrite);
	for (int i=0;i<32;i++){
		sprintf(hname,"hRawQDC%d",i);
		hRawQDC[i]->Write(hname,TObject::kOverwrite);
	}

	// done!
	RootFile->Close();
	cout << "\nWrote ROOT file." << endl;
}