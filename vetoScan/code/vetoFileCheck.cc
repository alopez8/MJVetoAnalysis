/*
	vetoFileCheck.C
	Clint Wiseman, USC/Majorana
	August 2015.

	Macro to check a list of run numbers and determine if the built files exist.
	Will also check if files have been blinded and aren't readable.

	Usage: 
	root[0] .X vetoFileCheck.C 
*/

#include "vetoScan.hh"

using namespace std;

void vetoFileCheck(string Input)
{
	int mode = 1; // switch: 0 for local files, 1 for pdsf files

	// Input a list of run numbers
	if (Input == "") Input = "Default.txt";
	
	ifstream InputList;
	InputList.open(Input.c_str());
	Char_t GATFile[200];
	Char_t BuiltFile[200];

	double durationTotal = 0;
	Int_t run;
	while(!InputList.eof()){

		// initialize 
		InputList >> run;
		if (mode==0) sprintf(BuiltFile,"~/dev/datasets/builtVeto/OR_run%i.root",run);
		//else if (mode==1) sprintf(BuiltFile,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/P3KJR/OR_run%u.root",run); 
		else if (mode==1) sprintf(BuiltFile,"/global/project/projectdirs/majorana/data/mjd/surfprot/data/built/P3END/OR_run%u.root",run); 

		if (mode==0) sprintf(GATFile,"~/dev/datasets/builtVeto/mjd_run%i.root",run);
		//else if (mode==1) sprintf(GATFile,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/gatified/P3KJR/mjd_run%u.root",run); 
		else if (mode==1) sprintf(GATFile,"/global/project/projectdirs/majorana/data/mjd/surfprot/data/gatified/P3END/mjd_run%u.root",run); 



		//cout << "Checking for built & gatified runs, run " << run << endl;

		// if file doesn't exist, ROOT will fail to open it.
		TFile *f1 = new TFile(BuiltFile);
    	f1->Close();

		// Check also that the duration is not corrupted!
		Float_t duration = 0;
		TChain *MGTree = new TChain("MGTree");
		MGTree->AddFile(BuiltFile);
		MJTRun *MyRun = new MJTRun();
		MGTree->SetBranchAddress("run",&MyRun);
        MGTree->GetEntry(0);
        duration = MyRun->GetStopTime() - MyRun->GetStartTime();
        //if (duration >= 3595 && duration <= 3605) cout << run << endl;
		if (duration <= 0 || duration > 4000 ) {
        	printf("\nRun %i has duration %.0f, skipping file!\n\n",run,duration);
        	continue;
        }
	durationTotal+=duration;

    	TFile *f2 = new TFile(GATFile);
    	f2->Close();
    }
	cout << "Total duration: " << durationTotal << " seconds" << endl; 
}
