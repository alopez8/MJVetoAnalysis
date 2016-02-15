// MJD veto analysis suite.
// Uses the November '15 MJD built data format.
// 
// Run with: "make" and ./vetoScan InputFile.txt
// 
// Clint Wiseman, University of South Carolina
// 1/23/2016


#ifndef VETOSCAN_H_GUARD
#define VETOSCAN_H_GUARD

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"
#include "TLine.h"

#include "MJVetoEvent.hh"
#include "GATDataSet.hh"

using namespace std;

// Processing (defined in vetoTools.cc)
void Test();
long GetStartUnixTime(GATDataSet ds);
long GetStopUnixTime(GATDataSet ds);
int GetNumFiles(string arg);
int color(int i);
int PanelMap(int i);
//int* GetQDCThreshold(int* thresh);

// Analysis
void vetoFileCheck(string Input = "");
void vetoPerformance(string file);
void vetoThreshFinder(string arg, int forceThresh = 500);

#endif