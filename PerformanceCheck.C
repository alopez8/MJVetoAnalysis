#ifndef __CINT__
#include <iostream>
#include <fstream>
#include <sstream>
#include <fstream>

#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"

#include "GATDataSet.hh"
#include "MJVetoEvent.hh"
#endif

using namespace std;

// ==================================================
// Processing Functions
// ==================================================
long GetStartUnixTime(GATDataSet ds)
{
	TChain *c = ds.GetVetoChain();
	MJTRun *runInfo = new MJTRun();
	c->SetBranchAddress("run",&runInfo);
	c->GetEntry(0);
	return (long)runInfo->GetStartTime();
}

long GetStopUnixTime(GATDataSet ds)
{
	TChain *c = ds.GetVetoChain();
	MJTRun *runInfo = new MJTRun();
	c->SetBranchAddress("run",&runInfo);
	c->GetEntry(0);
	return (long)runInfo->GetStopTime();
}

int GetNumFiles(string arg)
{
	int run = 0; 

	ifstream InputList;
	InputList.open(arg.c_str());
	if(!InputList.good()) {
    	cout << "Couldn't open " << arg << " !" << endl;
    	return 0;
    }

	int filesToScan = 0;
  	while(!InputList.eof()) { 
  		InputList >> run; 
  		filesToScan++; 
  	}
  	cout << "Scanning " << filesToScan << " files." << endl;	  	
  	InputList.close();
  	return filesToScan;
}
// ROOT color wheel: https://root.cern.ch/root/html/TColor.html
int color(int i)
{
	if (i == 0) return kRed-2;  
	if (i == 1) return kBlue;  
	if (i == 2) return kBlue+6;
	if (i == 3) return kAzure; 
	if (i == 4) return kMagenta;
	if (i == 5) return kGreen;
	if (i == 6) return kSpring;
	return 0;
}

// ==================================================
// Processing Functions
// ==================================================
// Might move this to a separate .cc file.
// Intended to emulate builtVetoCal with extra error checking.
// Check as much as possible without muon identification.
void PerformanceCheck(string file)
{	
	// Set up output file
	string filename = file;
	filename.erase(filename.find_last_of("."), string::npos); 
	char OutputFile[200];
	sprintf(OutputFile,"./output/VPC_%s.root",filename.c_str());
	TFile *RootFile = new TFile(OutputFile, "RECREATE");	
  	TH1::AddDirectory(kFALSE); // Global flag: "When a (root) file is closed, all histograms in memory associated with this file are automatically deleted."
	int filesToScan = GetNumFiles(file);

	// Plots
 	bool PlotMultiplicity = true;	// flag to plot multiplicity for EACH RUN
	TH1D *TotalMultiplicity = new TH1D("TotalMultiplicity","Events over threshold",32,0,32);
 	TotalMultiplicity->GetXaxis()->SetTitle("number of panels hit");
	TH1F *hRawQDC[32];  
	TH1F *hCutQDC[32];
	TH1F *hThreshQDC[32];
	TH1F *hLEDCutQDC[32];
 	const int nqdc_bins=4200;
	const int ll_qdc=0;
	const int ul_qdc=4200;
	Char_t hname[50];
	for (Int_t i=0; i<32; i++){
		sprintf(hname,"hRawQDC%d",i);
		hRawQDC[i] = new TH1F(hname,hname,nqdc_bins,ll_qdc,ul_qdc);
		sprintf(hname,"hCutQDC%d",i);
		hCutQDC[i] = new TH1F(hname,hname,nqdc_bins,ll_qdc,ul_qdc);
		sprintf(hname,"hThreshQDC%d",i);
		hThreshQDC[i] = new TH1F(hname,hname,500,ll_qdc,500);
		sprintf(hname,"hLEDCutQDC%d",i);
		hLEDCutQDC[i] = new TH1F(hname,hname,500,ll_qdc,500);
	}

  	// Loop over files in dataset
  	int run = 0;			// run number
  	long start = 0;			// unix timestamps
  	long stop = 0;
	int filesScanned = 0;	// counter
	vector<int> xaxis[7];	// save run numbers for TGraph
	vector<int> yaxis[7];	// save veto bit error counts 

	int vErrorCountTotal[11] = {0};
	long TotalEntries = 0;
	long totalDuration = 0;

	// run error statistics
	int errorCount0 = 0;
	int errorCount1 = 0;
	int errorCount2 = 0;
	int errorCount3 = 0;
	int errorCount4 = 0;
	int errorCount5 = 0;
	int errorCount6 = 0;
	int errorCount7 = 0;
	int errorCount8 = 0;
	int errorCount9 = 0;
	int errorCount10 = 0;

	// event error statistics
	int eventError0 = 0;
	int eventError1 = 0;
	int eventError2 = 0;
	int eventError3 = 0;
	int eventError4 = 0;
	int eventError5 = 0;
	int eventError6 = 0;
	int eventError7 = 0;
	int eventError8 = 0;
	int eventError9 = 0;
	int eventError10 = 0;

	// Keep track of runs with bad errors
	ofstream badFiles;
	badFiles.open("./output/P3KJR_BadRuns.txt");
	
  	
  	ifstream InputList;
  	InputList.open(file.c_str());
	while(!InputList.eof())
	{
		InputList >> run;
		for(int i = 0; i < 7; i++) xaxis[i].push_back(run);
		filesScanned++;
		GATDataSet ds(run);
		start = GetStartUnixTime(ds);
		stop = GetStopUnixTime(ds);
			
		// standard veto initialization block
		// with some error checking
		TChain *v = ds.GetVetoChain();
		long vEntries = v->GetEntries();
		MJTRun *vRun = new MJTRun();
		MGTBasicEvent *vEvent = new MGTBasicEvent(); 
		unsigned int mVeto = 0;
		uint32_t vBits = 0;
		if (v->GetListOfBranches()->FindObject("run")) 
			v->SetBranchAddress("run",&vRun);
		else { cout << "Couldn't find branch: \"run\"" << endl; break; }
		if (v->GetListOfBranches()->FindObject("mVeto")) 
			v->SetBranchAddress("mVeto",&mVeto);
		else { cout << "Couldn't find branch: \"mVeto\"" << endl; break; }
		if (v->GetListOfBranches()->FindObject("vetoEvent")) 
			v->SetBranchAddress("vetoEvent",&vEvent);
		else { cout << "Couldn't find branch: \"vetoEvent\"" << endl; break; }
		if (v->GetListOfBranches()->FindObject("vetoBits")) 
			v->SetBranchAddress("vetoBits",&vBits);
		else { cout << "Couldn't find branch: \"vetoBits\"" << endl; break; }
		v->GetEntry(0);
		MJTVetoData *vData[32];
		for (int i=0; i<32; i++) { vData[i] = dynamic_cast<MJTVetoData*>(vEvent->GetDetectorData()->At(i)); }
		long SBCOffset = vData[0]->GetSBCTSsec() - start;
		//long SBCOffset = -15049;
		//vRun->ListRunBits();	// Dump run bits
		TotalEntries += vEntries;

		printf("\n ========= Scanning Run %i: %li entries. =========\n",run,vEntries);

		// only write these if their bools are set = true
		//TH1D *CorruptionInTime = new TH1D("CorruptionInTime","corrupted entries during run",(int)duration/5,0,(int)duration);
		//CorruptionInTime->GetXaxis()->SetTitle("time (5 sec / bin)");
		TH1D *OneRunMultiplicity = new TH1D("multiplicity","multiplicity of veto entries",32,0,32);
		OneRunMultiplicity->GetXaxis()->SetTitle("number of panels hit");		
		
		// Run-level checks
		bool highCorruption = false;
		
		// veto error variables
		bool vBitTripped = false;
		int eventsSkipped = 0;
		int vPrevErrorCount[11] = {0};
		int vErrorCount[11] = {0};
		bool LEDsOff = true;

		bool error0 = false;
		bool error1 = false;
		bool error2 = false;
		bool error3 = false;
		bool error4 = false;
		bool error5 = false;
		bool error6 = false;
		bool error7 = false;
		bool error8 = false;
		bool error9 = false;
		bool error10 = false;

		MJVetoEvent prev;
		
		// Loop over entries
		//for (long i = 0; i < 10; i++) 
		//for (long i = vEntries-10; i<vEntries; i++)
		for (int i = 0; i < vEntries; i++)
		{
			v->GetEntry(i);

			// Fill MJVetoEvent object
			MJVetoEvent veto;
			veto.SetQDCThreshold(run,1);	// kept separate to allow the user to specify their own.
    		int isGood = veto.WriteEvent(vRun,vEvent,vData,vBits,run);

    		// Check SBC Time and adjust for offset.
			veto.timeSBC = veto.timeSBC - SBCOffset - start;
			double BareSBC = (long)vData[0]->GetSBCTSsec()-start;
			//printf("\nstart: %li  SBC offset: %li  Bare SBC - start: %.2f  Scaler: %.8f  SBC: %.8f  Diff: %.8f\n"
			//	,start,SBCOffset,BareSBC,veto.timeSec,veto.timeSBC,veto.timeSec-veto.timeSBC);

    		//veto.Print();
    		//cout << "Wrote event " << i << ".  isGood: " << isGood << endl;
			
			// Check veto bits for errors (defined in MGDO/Majorana/MJTypes)
			if (MJBits::GetBit(vBits, MJVetoBits::kMissingChannels)) 	{vErrorCount[0]++; 	error0=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kExtraChannels)) 		{vErrorCount[1]++; 	error1=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kDuplicateChannel)) 	{vErrorCount[2]++; 	error2=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kBadTimeStamp)) 		{vErrorCount[3]++; 	error3=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kQDCOutOfSequence)) 	{vErrorCount[4]++; 	error4=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kScalerOnly)) 		{vErrorCount[5]++; 	error5=true;}
			if (MJBits::GetBit(vBits, MJVetoBits::kHWCountMismatch)) 	{vErrorCount[6]++; 	error6=true;}
			if (isGood == -1) {vErrorCount[8]++;   	error8=true;}// QDC entries don't exist but kScalerOnly was not tripped
    		if (isGood == -2) {vErrorCount[9]++;   	error9=true;}// Unknown card in crate
    		if (isGood == -3) {vErrorCount[10]++;  	error10=true;}// vEvent run number doesn't match the current run

			for (int j = 0; j < 10; j++) 
			{ 
				// Show first instance of a bit being thrown.
				if (vErrorCount[j] == 1 && vPrevErrorCount[j] == 0) 
				
				// Show all instances of bits being thrown.
				//if (vErrorCount[j] > 0) 

				// Show the first and second instance of bits being thrown.
				if ((vErrorCount[j] == 1 && vPrevErrorCount[j] == 0) || (vErrorCount[j] == 2 && vPrevErrorCount[j] == 1))
				{
					vBitTripped = true;
					printf(" New Error: Bit %i  Count: %i  Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
						,j,vErrorCount[j],veto.run,i,veto.timeSec,((double)i/vEntries)*100); 
					
					// Info for fIndex and EventCount mismatches.
					printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
						,veto.QDCchans,veto.ScalerIndex,veto.QDC1Index,veto.QDC2Index
						,veto.ScalerIndex - veto.QDC1Index, veto.ScalerIndex - veto.QDC2Index
						,veto.SEC,veto.QEC);

					// Show previous event.
					if (i>0) {
						printf("   Last event: Run: %i  Entry: %li  Time: %.5f  %.2f%%.\n"
							,prev.run,i-1,prev.timeSec,((double)(i-1)/vEntries)*100);
						
						printf("   Channels: %lu  scaler_Ind: %li  qdc1_Ind: %li  qdc2_Ind: %li  s-QDC1: %li  s-qdc2: %li  SEC: %li  QEC: %li\n\n"
							,prev.QDCchans,prev.ScalerIndex,prev.QDC1Index,prev.QDC2Index
							,prev.ScalerIndex - prev.QDC1Index, prev.ScalerIndex - prev.QDC2Index
							,prev.SEC,prev.QEC);
					}
				}
			}
			memcpy(vPrevErrorCount, vErrorCount, sizeof(vPrevErrorCount));
			if (i > 0) prev = veto;
			else prev.Clear();

			if (isGood == 0 || isGood == -1 || isGood == -2 || isGood == -3) {
				eventsSkipped++;
				veto.Clear();
				prev = veto;
				continue;
			}
 
			// Check for LED entries in a simple way
			if (veto.multip > 25) LEDsOff = false;

			// Compare multiplicities
			//printf("mVeto: %i  multip: %i\n",mVeto,veto.multip);

			// Multiplicity and QDC plots
			TotalMultiplicity->Fill(veto.multip);			
	        if (PlotMultiplicity) OneRunMultiplicity->Fill(veto.multip);
			for (int k = 0; k < 32; k++) {
				hRawQDC[k]->Fill(veto.QDC[k]);
				if (veto.QDC[k] < 500) hThreshQDC[k]->Fill(veto.QDC[k]);
				if (veto.QDC[k] > veto.QDC_SWThresh[k]) hCutQDC[k]->Fill(veto.QDC[k]);
			}
	
		} // End loop over entries.
		
		printf(" ================= End of Run %i. =================\n",run);

		printf(" Duration: %i seconds.\n",stop-start);

		if (error0) { errorCount0++; badFiles << run << endl; }
		if (error1) errorCount1++;
		if (error2) errorCount2++;
		if (error3) errorCount3++;
		if (error4) errorCount4++;
		if (error5) errorCount5++;
		if (error6) errorCount6++;
		if (error8) errorCount8++;
		if (error9) errorCount9++;
		if (error10) errorCount10++;

		if (!LEDsOff) { vErrorCount[7] = 1; erroCount7++;}
		if (PlotMultiplicity) {
			char outfile2[200];
			sprintf(outfile2,"Multiplicity_Run%i",run);
			OneRunMultiplicity->Write(outfile2,TObject::kOverwrite);
		}

		if (vBitTripped) {
			double scalerCorruption = (((double)vErrorCount[3])/vEntries)*100;
			printf(" Skipped %i events. (%.2f%% of total). \n Veto Errors: \n",eventsSkipped,((double)eventsSkipped/vEntries)*100);
			if (vErrorCount[0] > 0) printf(" 0: Missing Channels (< 32 veto datas in event) : %i\n",vErrorCount[0]);
			if (vErrorCount[1] > 0) printf(" 1: Extra Channels (> 32 veto datas in event) : %i\n",vErrorCount[1]);
			if (vErrorCount[2] > 0) printf(" 2: Duplicate Channels (any channel shows up multiple times) : %i\n",vErrorCount[2]);
			if (vErrorCount[3] > 0) printf(" 3: Bad Timestamps: %i of %li entries: %.2f%%.\n",vErrorCount[3],vEntries,scalerCorruption);
			if (vErrorCount[4] > 0) printf(" 4: Scaler found w/ no QDC data: %i\n",vErrorCount[4]);
			if (vErrorCount[5] > 0) printf(" 5: fIndex Errors: (QDCIndex - ScalerIndex != 1 or 2) : %i\n",vErrorCount[5]);
			if (vErrorCount[6] > 0) printf(" 6: EventCount Errors: (SEC - QEC != 1 or 2) : %i\n",vErrorCount[6]);
			if (vErrorCount[7] > 0) printf(" 7: LED's not activated this run : %i\n",vErrorCount[7]);
			if (vErrorCount[8] > 0) printf(" 8: No QDC entries but bit 4 not thrown : %i\n",vErrorCount[8]);
			if (vErrorCount[9] > 0) printf(" 9: Unknown card in crate : %i \n",vErrorCount[9]);
			if (vErrorCount[10]> 0) printf(" 10: Incorrect run number in MJTRun vRun : %i \n",vErrorCount[10]);
			cout << " ====================================================\n";
		}

		eventError0 += vErrorCount[0];
		eventError1 += vErrorCount[1];
		eventError2 += vErrorCount[2];
		eventError3 += vErrorCount[3];
		eventError4 += vErrorCount[4];
		eventError5 += vErrorCount[5];
		eventError6 += vErrorCount[6];
		eventError7 += vErrorCount[7];
		eventError8 += vErrorCount[8];
		eventError9 += vErrorCount[9];
		eventError10 += vErrorCount[10];

		// Push veto bits into yaxis for TGraph
		for (int q = 0; q < 7; q++) {
			yaxis[q].push_back(vErrorCount[q]);
		}

	} // End loop over file list

	printf ("\n \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ END OF SCAN //////////////////////\n");

	printf("Error Summary: \n");
	printf("Runs scanned: %i\n",filesToScan);
	printf(" 0: Missing Channels (< 32 veto datas in event) : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%) \n"
		,errorCount0,filesToScan,100*(double)errorCount0/filesToScan,eventError0,TotalEntries,100*(double)eventError0/TotalEntries);
	printf(" 1: Extra Channels (> 32 veto datas in event) : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%) \n"
		,errorCount1,filesToScan,100*(double)errorCount1/filesToScan,eventError1,TotalEntries,100*(double)eventError1/TotalEntries);
	printf(" 2: Duplicate Channels (any channel shows up multiple times) : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,errorCount2,filesToScan,100*(double)errorCount2/filesToScan,eventError2,TotalEntries,100*(double)eventError2/TotalEntries);
	printf(" 3: Bad Timestamps: %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,errorCount3,filesToScan,100*(double)errorCount3/filesToScan,eventError3,TotalEntries,100*(double)eventError3/TotalEntries);
	printf(" 4: Scaler found w/ no QDC data: %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,errorCount4,filesToScan,100*(double)errorCount4/filesToScan,eventError4,TotalEntries,100*(double)eventError4/TotalEntries);
	printf(" 5: fIndex Errors: (QDCIndex - ScalerIndex != 1 or 2) : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,errorCount5,filesToScan,100*(double)errorCount5/filesToScan,eventError5,TotalEntries,100*(double)eventError5/TotalEntries);
	printf(" 6: EventCount Errors: (SEC - QEC != 1 or 2) : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,errorCount6,filesToScan,100*(double)errorCount6/filesToScan,eventError6,TotalEntries,100*(double)eventError6/TotalEntries);
	printf(" 7: LED's not activated this run : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,errorCount7,filesToScan,100*(double)errorCount7/filesToScan,eventError7,TotalEntries,100*(double)eventError7/TotalEntries);
	printf(" 8: No QDC entries but bit 4 not thrown : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,errorCount8,filesToScan,100*(double)errorCount8/filesToScan,eventError8,TotalEntries,100*(double)eventError8/TotalEntries);
	printf(" 9: Unknown card in crate : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,errorCount9,filesToScan,100*(double)errorCount9/filesToScan,eventError9,TotalEntries,100*(double)eventError9/TotalEntries);
	printf(" 10: Incorrect run number in MJTRun vRun : %i of %i runs (%.2f %%) and %i of %i events (%.2f %%)\n"
		,errorCount10,filesToScan,100*(double)errorCount10/filesToScan,eventError10,TotalEntries,100*(double)eventError10/TotalEntries);


	badFiles.close();



	//TotalCorruptionInTime->Write("TotalCorruptionInTime",TObject::kOverwrite);
	TotalMultiplicity->Write("TotalMultiplicity",TObject::kOverwrite);

	// Plot errors vs. run number
	// Find graph parameters
	int lastxaxis = 0;
	for (int i = 0; i < xaxis[0].size(); i ++ ) 
		lastxaxis = xaxis[0][i];
	int lastyaxis = 0;
	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < yaxis[i].size(); j++) {
			if (yaxis[i][j] > lastyaxis) lastyaxis = yaxis[i][j];
		}
	}
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas");
	c1->SetGrid();
	TGraph *g[7];
	TH1F *hr = c1->DrawFrame(xaxis[0][0],0,lastxaxis,lastyaxis);
	hr->SetXTitle("Run");
	hr->SetYTitle("Errors");
	hr->SetTitle("Veto Errors");
	TLegend *legend = new TLegend(0.8,0.58,0.9,0.9);
	for (int i = 0; i < 7; i++){ 
		g[i] = new TGraph(filesScanned,&(xaxis[i][0]),&(yaxis[i][0]));
		g[i]->SetLineColor(color(i));
		g[i]->SetMarkerColor(color(i));
		g[i]->SetMarkerStyle(20);
		g[i]->SetLineWidth(2);	// want to make error bars thicker
		g[i]->Draw("LP");
	}
	legend->AddEntry(g[0],"Missing Channels","lep");
	legend->AddEntry(g[1],"Extra Channels","lep");
	legend->AddEntry(g[2],"Duplicate Channels","lep");
	legend->AddEntry(g[3],"Bad Timestamps","lep");
	legend->AddEntry(g[4],"Scaler No QDC","lep");
	legend->AddEntry(g[5],"fIndex Mismatch","lep");
	legend->AddEntry(g[6],"EventCount Mismatch","lep");
	legend->Draw();
	c1->Write("ErrorPlot",TObject::kOverwrite);
	delete c1;

	// QDC plots & calibration table:
	// gaus: A gaussian with 3 parameters: f(x) = p0*exp(-0.5*((x-p1)/p2)^2)).
	TF1 *fits[32];
	TCanvas *vcan0 = new TCanvas("FittedCutQDC","cut & fitted QDC, panels 1-32",0,0,800,600);
	vcan0->Divide(8,4,0,0);
	TCanvas *vcan1 = new TCanvas("QDCThresholds","veto QDC thresholds, panels 1-32",0,0,800,600);
	vcan1->Divide(8,4,0,0);	
	Char_t buffer[2000];
	printf("\n  LED Peak Information:\n  Panel / Mean,error / Sigma,error / Chi-square/NDF (~1?)\n");
	for (Int_t i=0; i<32; i++)
	{	
		vcan0->cd(i+1);
		TVirtualPad *vpad0 = vcan0->cd(i+1); vpad0->SetLogy();
		hThreshQDC[i]->Write(); // write the low-QDC part of the spectrum separately
		hCutQDC[i]->Write();    // write the cut QDC

		// if we have entries above the cut threshold, fit them
		if (hCutQDC[i]->GetEntries() > 0) {
			hRawQDC[i]->Write();
			hCutQDC[i]->Fit("gaus","q");
			hCutQDC[i]->Draw();
			fits[i] = hCutQDC[i]->GetFunction("gaus");
			float NDF = (float)fits[i]->GetNDF();
			if (fits[i]->GetNDF() == 0) NDF = 0.0000001;
			sprintf(buffer,"%i  %.1f  %.1f  %.1f  %.1f  %.1f",i,
				fits[i]->GetParameter(1),fits[i]->GetParError(1),		// mean (LED center)
				fits[i]->GetParameter(2),fits[i]->GetParError(2),		// sigma
				fits[i]->GetChisquare()/NDF); // X^2/NDF (closer to 1, better fit.)
			cout << "  " << buffer << endl;
		}
		else {
			hRawQDC[i]->Write();	// just write the raw QDC without fitting
		}	    	
   		// plot low-QDC range separately
  		vcan1->cd(i+1);
  		TVirtualPad *vpad1 = vcan1->cd(i+1); vpad1->SetLogy();
  		hThreshQDC[i]->Draw();
	}
	vcan0->Write((filename+"_QDC").c_str());
	vcan1->Write((filename+"_ThreshQDC").c_str());
	

	// All done!
	RootFile->Close();
}