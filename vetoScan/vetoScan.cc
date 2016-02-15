// MJD veto analysis suite.
// Uses the November '15 MJD built data format.
// 
// Run with: "make" and ./vetoScan InputFile.txt
// 
// Clint Wiseman, University of South Carolina
// 1/23/2016

#include "vetoScan.hh"

using namespace std;

int main(int argc, char** argv) 
{
	// Input a file from the shell.
	if (argc < 2) {
		cout << "Usage: ./vetoScan <input_runs_file>" << endl;
		return 1;
	}
	// Might want to adjust output verbosity from shell args
	string arg = argv[1];

	// Test();
	// vetoFileCheck(arg);
	// vetoPerformance(arg);
	// vetoThreshFinder(arg,300);	// 300: force threshold

	cout << "\nCletus codes good." << endl;
}