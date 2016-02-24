#include "vetoScan.hh"

// I'm sick of programming in QDC thresholds by hand.
// Figure them out for me, computer!
//
// Also, try to catch when a QDC pedestal moves from run to run.
//
void vetoThreshFinder(string Input, bool runHistos){

	// Input a list of run numbers
	ifstream InputList(Input.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << Input << endl;
    	return;
    }

	// Strip off path and extension: use for output files.
	string Name = Input;
	Name.erase(Name.find_last_of("."),string::npos);
	Name.erase(0,Name.find_last_of("\\/")+1);

	// Return ONE picture of 32 panels' raw spectrum, 
	// with 32 big red vertical lines at the location
	// the program decided to place the threshold.
	//
	TH1F *hRawQDC[32];  
	int bins = 500;
	int lower = 0;
	int upper = 500;
	char hname[50];
	for (int i = 0; i < 32; i++) {
		sprintf(hname,"hRawQDC%d",i);
		hRawQDC[i] = new TH1F(hname,hname,bins,lower,upper);
	}
	bool pedestalShift = false;
	int runThresh[32] = {0};	// run-by-run threshold
	int prevThresh[32] = {0};	
		
	int run = 0;
	int filesScanned = 0;
	while(!InputList.eof()){

		// initialize 
		InputList >> run;
		GATDataSet *ds = new GATDataSet(run);

		// standard veto initialization block
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

		printf("\n========= Scanning Run %i: %li entries. =========\n",run,vEntries);

		// Use super-low QDC threshold for this.
		// This should cause all entries to have a multiplicity of 32
		int def[32] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

		// Run-by-run histograms, trying to catch a changing QDC pedestal.
		// These are not drawn.
		TH1F *hRunQDC[32];  
		for (int i = 0; i < 32; i++) {
			sprintf(hname,"hRunQDC%d",i);
			hRunQDC[i] = new TH1F(hname,hname,bins,lower,upper);
		}
		
		int highestMultip = 0;	// try to predict how many panels there are for this run.
		long skippedEvents = 0;
		long corruptScaler = 0;
		int isGood = 0;
		for (long i = 0; i < vEntries; i++) 
		{
			v->GetEntry(i);

			// Fill MJVetoEvent 
			MJVetoEvent veto;
			veto.SetQDCThreshold(def);	
	    	isGood = veto.WriteEvent(i,vRun,vEvent,vBits,run);
	    	if (CheckForBadErrors(veto,i,isGood,false)) continue;
	    	
	    	if (veto.badScaler) corruptScaler++;
	    	if (veto.multip > highestMultip) highestMultip = veto.multip;

	    	// Fill raw histogram under 500
	    	for (int q = 0; q < 32; q++) {
	    		hRawQDC[q]->Fill(veto.QDC[q]);
	    		hRunQDC[q]->Fill(veto.QDC[q]);
			}
		}
		if (skippedEvents > 0) printf("Skipped %li of %li entries.\n",skippedEvents,vEntries);
		if (corruptScaler > 0) printf("Corrupt scaler: %li of %li entries (%.2f%%) .\n"
			,corruptScaler,vEntries,100*(double)corruptScaler/vEntries);

		// Set up a 32-panel plot for each run if runHistos = true.
		TCanvas *runHist = new TCanvas("run","veto low QDC",800,600);
		runHist->Divide(8,4);

		// Calculate the run-by-run threshold location.
		// Throw a warning if a pedestal shifts by more than 5%.
		for (int c = 0; c < 32; c++) 
		{
			runThresh[c] = FindQDCThreshold(hRunQDC[c],c,false);
			double ratio = (double)runThresh[c]/prevThresh[c];
			if (filesScanned !=0 && (ratio > 1.1 || ratio < 0.9)) 
			{
				printf("Warning! Found pedestal shift! Panel: %i  Previous: %i  This run: %i \n"
					,c,prevThresh[c],runThresh[c]);
				pedestalShift = true;
			}

			// save threshold for next scan
			prevThresh[c] = runThresh[c];

			if (runHistos) {
				runHist->cd(c+1);
				hRunQDC[c]->Draw();
			}

		}
		if (runHistos) {
			char runName[200];
			sprintf(runName,"./output/QDCLow_%s_%i.C",Name.c_str(),run);
			runHist->Print(runName);
		}

		// clear memory
		delete runHist;
		for (int c=0;c<32;c++) delete hRunQDC[c];
		
		// done with this run
		filesScanned++;
	}
	cout << "\n==================== End of Scan. ====================\n\n";

	// Output: Find the QDC Pedestal location in each channel.
	// Give a threshold that is 20 QDC above this location, and output a plot
	// that confirms this choice.

	int thresh[32] = {9999};
	TCanvas *vcan0 = new TCanvas("vcan0","veto QDC thresholds, panels 1-32",0,0,800,600);
	vcan0->Divide(8,4,0,0);
	for (int i=0; i<32; i++)
	{
		vcan0->cd(i+1);
		TVirtualPad *vpad0 = vcan0->cd(i+1); vpad0->SetLogy();

		// find overall threshold for this panel
		thresh[i] = FindQDCThreshold(hRawQDC[i],i,true);

		// reset histo range and draw
		hRawQDC[i]->GetXaxis()->SetRange(lower,upper);
		hRawQDC[i]->Draw();

		double ymax = hRawQDC[i]->GetMaximum();
		TLine *line = new TLine(thresh[i],0,thresh[i],ymax+10);
		line->SetLineColor(kRed);
		line->SetLineWidth(2.0);
		line->Draw();		
	}

	// End of Scan Output

	if (pedestalShift) {
		printf("Warning: Found a pedestal shift by more than 5%% of its previous value.\n");
		printf("         You may want to go back and examine the output.\n");
		printf("         It can be caused by including a run with very few entries.\n");
	}

	cout << "\nMeasured QDC thresholds:\n\n";
	for (int r = 0; r < 32; r++) printf("[%i] %i  ",r,thresh[r]);
	cout << "\n\n";

	cout << "int thresh[32] = {";
	for (int r = 0; r < 31; r++) cout << thresh[r] << ",";
	cout << thresh[31] << "};" << "\n\n";
    	
   	// Write canvas
	Char_t OutputName[200];	
	sprintf(OutputName,"./output/SWThresh_%s.C",Name.c_str());
	vcan0->Print(OutputName);
	
}