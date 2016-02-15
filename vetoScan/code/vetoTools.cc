// ==================================================
// Processing Functions
// ==================================================

#include "vetoScan.hh"

using namespace std;

void Test()
{
	string s1("../somepath/somemorepath/somefile.ext");
	string s2("..\\somepath\\somemorepath\\somefile.ext");
	cout << s1.substr(0, s1.find_last_of("\\/")) << endl;
	cout << s2.substr(0, s2.find_last_of("\\/")) << endl;

	string s3("./runs/P3K93_VetoOnlyRuns.txt");
	cout << s3.substr(s3.find_last_of("\\/")+1,string::npos) << endl;
}

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

// For tagging plane-based coincidences.
// This uses a zero-indexed map.
int PanelMap(int i){
	// 0: Lower Bottom
	// 1: Upper Bottom
	// 2: Inner Top
	// 3: Outer Top
	// 4: Inner North
	// 5: Outer North
	// 6: Inner South
	// 7: Outer South
	// 8: Inner West
	// 9: Outer West
	// 10: Inner East
	// 11: Outer East 
			
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

	else if (i == 17) return 3; // Top outer
	else if (i == 18) return 3; // Top outer
	else if (i == 20) return 2; // Top inner
	else if (i == 21) return 2; // Top inner

	else if (i == 15) return 5; // North outer
	else if (i == 16) return 5; // North outer
	else if (i == 19) return 4; // North inner
	else if (i == 23) return 4; // North inner

	else if (i == 24) return 6; // South inner
	else if (i == 25) return 5; // South outer
	else if (i == 26) return 6; // South inner
	else if (i == 27) return 5; // South outer
	
	else if (i == 12) return 8; // West inner
	else if (i == 13) return 8; // West inner
	else if (i == 14) return 9; // West outer
	else if (i == 22) return 9; // West outer
	
	else if (i == 28) return 10; // East inner
	else if (i == 29) return 11; // East outer
	else if (i == 30) return 10; // East inner
	else if (i == 31) return 11; // East outer

	else return -1;
}

/*
int GetQDCThreshold(int* thresh){
	
	cout << 1 << endl;
	
	int result[32] = {0};

	// Default QDC threshold 
	int def[32] = {500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,
		500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,};

	// Before 32 veto panels (JDY_01.txt)
	// After 32 veto panels (JDY_02.txt)
	// P3K93 (veto only, 24 panels)
	// P3KJR (24 panels)
	// P3KJR (32 panels)
	
	if (!thresh) {
		cout << 2 << endl;
		memcpy(result, def, sizeof(def));
	}
	// Custom QDC threshold 
	else if (thresh) {
		cout << 3 << endl;
		memcpy(result, thresh, sizeof(def));
	}
	else std::cout << "Couldn't set threshold!" << std::endl;

	cout << 4 << endl;
	return result;
}
*/