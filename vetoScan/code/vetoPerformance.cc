#include "vetoScan.hh"

using namespace std;

void vetoPerformance(string Input) 
{
	// Input a list of run numbers
	ifstream InputList(Input.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << Input << endl;
    	return;
    }
	
	const int nErrs = 18;
	int globalErrorCount[nErrs] = {0};

	int run = 0;
	long totEntries = 0;

	while(!InputList.eof())
	{
		InputList >> run;
		
		// standard initialization
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
		// long start = (long)vRun->GetStartTime();
		// long stop = (long)vRun->GetStopTime();
		// duration = ds->GetRunTime()/CLHEP::second;	

		int isGood = 1;
		int errorCount[nErrs] = {0};
		totEntries += vEntries;

		printf("========== Scanning run %i, %li entries.=============\n",run,vEntries);
		for (int i = 0; i < vEntries; i++)
		{
			v->GetEntry(i);
			MJVetoEvent veto;
			veto.SetQDCThreshold();	
	    	isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run,true);
	    	if (isGood != 1) 
	    	{	    		
	    		for (int j=0; j<nErrs; j++) if (veto.error[j]==1) errorCount[j]++;
	    		cout << endl;
	    		cout << "Entry: " << i << endl;
	    		veto.Print();
	    		continue;
	    	}

	    	// DON'T FORGET to code bad timestamps back in as "OK!"
	    	// They are being skipped right now.

	    	// Check SBC Time and adjust for offset.
			//veto.timeSBC = veto.timeSBC - SBCOffset - start;
			//double BareSBC = (long)vData[0]->GetSBCTSsec()-start;
			//printf("\nstart: %li  SBC offset: %li  Bare SBC - start: %.2f  Scaler: %.8f  SBC: %.8f  Diff: %.8f\n"
			//	,start,SBCOffset,BareSBC,veto.timeSec,veto.timeSBC,veto.timeSec-veto.timeSBC);

		}
		cout << "================ End Run " << run << ". ===================\n";
		for (int i = 0; i < nErrs; i++) {
			if (errorCount[i] > 0) {
				printf("%i: %i errors\t(%.2f%% of total)\n",i,errorCount[i],100*(double)errorCount[i]/vEntries);
				globalErrorCount[i] += errorCount[i];
			}
		}
	}

	cout << "\n\n============== END OF SCAN. ==================\n";
	for (int i = 0; i < nErrs; i++) {
		if (globalErrorCount[i] > 0) {
			printf("%i: %i errors\t(%.2f%% of total)\n"
				,i,globalErrorCount[i],100*(double)globalErrorCount[i]/totEntries);
		}
	}
}