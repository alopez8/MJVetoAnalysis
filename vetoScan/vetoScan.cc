// MJD veto analysis suite.
// Uses the November '15 MJD built data format.
// 
// Clint Wiseman, University of South Carolina
// 1/23/2016

#include "vetoScan.hh"

using namespace std;

static const char Usage[] =
"\n"
"Usage: vetoScan [options].\n"
"     REQUIRED for most routines : \n"
"     -F (--file) ./path/to/your/runList.txt\n\n"
"Additional options:\n"
"     -h (--help) : Print usage info\n"
"     -S (--serial) : Set the part number (P3JDY, etc.)  REQUIRED to use checkFiles.\n"
"     -f (--checkFiles) : Check that files exist.\n"
"                       : Options: `checkBuilt`, `checkGAT`, `checkBoth`, `checkGDS`, `checkAll`\n" 
"                       : (checkBuilt and checkGAT both require PDSF)\n"
"     -T (--swThresh) : Set QDC software threshold using `vetoSWThresholds.txt`\n"
"     -p (--perfCheck) : Veto performance check (data quality).\n"
"                      : Option: `runs`, `totals`\n"
"                      : If -T is specified, user picks which SW thresholds to use.\n"
"\n";


int main(int argc, char** argv) 
{
	if (argc == 1) {
		cout << Usage;
	}

	// Parse command line arguments with getopt_long:
	// http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html
	//
	string file = "", partNum = "", threshName = "";
	bool findThresh=0;
	bool perfCheck=0, fileCheck=0;
	bool checkBuilt=0,checkGAT=0,checkGDS=0;
	bool runBreakdowns=0;
	//
	int c;
	int option_index = 0;
	while (1)
	{
		static struct option long_options[] = {
			{"help", no_argument, 0, 'h'},
			{"file", required_argument, 0, 'F'},
			{"serial",required_argument, 0, 'S'},
			{"checkFiles", required_argument, 0, 'f'},
			{"swThresh", required_argument, 0, 'T'},
			{"perfCheck", required_argument, 0, 'p'},
		};

		c = getopt_long (argc, argv, "hF:S:f:T:p:",long_options,&option_index);
		if (c == -1) break;

		switch (c)
		{
		case 'h':
			cout << Usage;
			break;
		case 'F':
			file = string(optarg);
			cout << "Using input file: " << file << endl;
			break;
		case 'S':
			partNum = string(optarg);
			cout << "Using part number: " << partNum << endl;
			break;
		case 'f':
			fileCheck=1;
			if (string(optarg) == "checkBuilt") 		checkBuilt=1;
			else if (string(optarg) == "checkGAT")		checkGAT=1;
			else if (string(optarg) == "checkGDS")		checkGDS=1;
			else if (string(optarg) == "checkBoth") {	checkBuilt=1; checkGAT=1; }
			else if (string(optarg) == "checkAll")  {	checkBuilt=1; checkGAT=1; checkGDS=1; }
			printf("Checking for ... built? %i  gatified? %i  GATDataSet? %i\n",checkBuilt,checkGAT,checkGDS);
			break;
		case 'T':
			threshName = string(optarg);
			cout << "Using thresholds from " << threshName << endl;
			break;
	    case 'p': 
	    	perfCheck=1; 
	    	if (string(optarg) == "runs") runBreakdowns=1;
	    	else if (string(optarg) == "totals") runBreakdowns=0;
	    	break;
		case '?':
		    if (isprint (optopt))  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
		    else fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
		    cout << "Call vetoScan -h for usage options.\n";
		    return 1;
		default:
		    abort ();
		}
	}

	// =======================================================
	// Run selected routines
	//
	int thresh[32] = {0};

	if (fileCheck) 	vetoFileCheck(file,partNum,checkBuilt,checkGAT,checkGDS);
	if (findThresh)	vetoThreshFinder(file);
	if (perfCheck)	
	{	
		if (threshName != "") GetQDCThreshold(file,thresh,threshName);
		else GetQDCThreshold(file,thresh);
		vetoPerformance(file,thresh,runBreakdowns);
	}
	// =======================================================

	cout << "\nCletus codes good." << endl;

}