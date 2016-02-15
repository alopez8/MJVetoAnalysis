#include "vetoScan.hh"

// I'm sick of programming in QDC thresholds by hand.
// Figure them out for me, computer!
//
// Also, try to catch when a QDC pedestal moves from run to run.
//
void vetoThreshFinder(string Input, int forceThresh){

	// Return ONE picture of 32 panels' raw spectrum, 
	// with 32 big red vertical lines at the location
	// the program decided to place the threshold.
	//
	TH1F *hRawQDC[32];  
	int bins = 250;
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
		
	// Input a list of run numbers
	ifstream InputList(Input.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << Input << endl;
    	return;
    }
	int run = 0;
	int filesScanned = 0;

	// Strip off path and extension: use for output files.
	string Name = Input;
	Name.erase(Name.find_last_of("."),string::npos);
	Name.erase(0,Name.find_last_of("\\/")+1);


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
		// This will cause all entries to have a multiplicity of 32
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
	    	
	    	// Simple error checking
	    	// (Remember, WriteEvent returns easily.)
	    	if (veto.badScaler) corruptScaler++;
	    	if (veto.multip > highestMultip) highestMultip = veto.multip;
	    	if (isGood != 1) {
	    		//cout << "Skipped event " << i << ". Checking bits:" << endl;
	    		//veto.CheckVetoBits();
	    		skippedEvents++;
	    		continue;
	    	}

	    	// Fill raw histogram under 500
	    	for (int q = 0; q < 32; q++) {
	    		hRawQDC[q]->Fill(veto.QDC[q]);
	    		hRunQDC[q]->Fill(veto.QDC[q]);
			}
		}
		if (skippedEvents > 0) printf("Skipped %li of %li entries.\n",skippedEvents,vEntries);
		if (corruptScaler > 0) printf("Corrupt scaler: %li of %li entries (%.2f%%) .\n"
			,corruptScaler,vEntries,100*(double)corruptScaler/vEntries);

		// Calculate the run-by-run threshold location.
		// Throw a warning if a pedestal shifts by more than 5%.
		for (int c = 0; c < 32; c++) {
			int bin = hRunQDC[c]->GetMaximumBin();
			double xval = hRunQDC[c]->GetXaxis()->GetBinCenter(bin);
			runThresh[c] = xval+20;

			double ratio = (double)runThresh[c]/prevThresh[c];
			if (filesScanned !=0 && (ratio > 1.05 || ratio < 0.95)) {
				printf("Warning! Found pedestal shift!\nPanel: %i  Previous: %i  This run: %i \n"
					,c,prevThresh[c],runThresh[c]);
				pedestalShift = true;
			}

			// save threshold for next scan
			prevThresh[c] = runThresh[c];

			// clear the histo from memory
			delete hRunQDC[c];
		}
		
		// done with this run
		filesScanned++;
	}
	cout << "\n==================== End of Scan. ====================\n\n";

	// Output: Find the QDC Pedestal location in each channel.
	// Give a threshold that is 20 QDC above this location, and output a plot
	// that confirms this choice.

	int thresh[32] = {9999};
	//TF1 *fits[32];	// try to fit the QDC pedestal
	TCanvas *vcan0 = new TCanvas("vcan0","veto QDC thresholds, panels 1-32",0,0,800,600);
	vcan0->Divide(8,4,0,0);
	for (int i=0; i<32; i++)
	{
		vcan0->cd(i+1);
		TVirtualPad *vpad0 = vcan0->cd(i+1); vpad0->SetLogy();

		// Method 1: Find QDC value of maximum bin.
		//
		int binmax = hRawQDC[i]->GetMaximumBin();
		double x = hRawQDC[i]->GetXaxis()->GetBinCenter(binmax);
		double ymax = hRawQDC[i]->GetMaximum();
		hRawQDC[i]->Draw();
		thresh[i] = x+20;
		TLine *line = new TLine(thresh[i],0,thresh[i],ymax+10);
		line->SetLineColor(kRed);
		line->SetLineWidth(2.0);
		line->Draw();		
		
		// Method 2: Try to fit the pedestal to a gaussian.
		// This does not work very well, the pedestal is too narrow.
		//
		// fits[i] = new TF1("fits","gaus",0,500);
		// fits[i]->SetParameter(1,x);	
		// fits[i]->SetParameter(2,5);	// sigma of 5qdc is probably quite generous
		// hRawQDC[i]->Fit(fits[i],"q");
		// hRawQDC[i]->Draw();
		// char buffer[2000];
		// if (hRawQDC[i]->GetEntries() > 0) 
		// {
		// 	double NDF = (double)fits[i]->GetNDF();
		// 	if (fits[i]->GetNDF() == 0) NDF = 0.0000001;
		// 	sprintf(buffer,"%i  %.1f  %.1f  %.1f  %.1f  %.1f",i,
		// 		fits[i]->GetParameter(1),fits[i]->GetParError(1),		// mean (LED center)
		// 		fits[i]->GetParameter(2),fits[i]->GetParError(2),		// sigma
		// 		fits[i]->GetChisquare()/NDF); // X^2/NDF (closer to 1, better fit.)
		// 	cout << "  " << buffer << endl;
		// }
	}

	// End of Scan Output

	if (pedestalShift) {
		printf("Warning: Found a pedestal shift by more than 5%% of its previous value.\n");
		printf("         You may want to go back and examine the output.\n");
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