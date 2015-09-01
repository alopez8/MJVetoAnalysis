/*
	CheckFiles.C
	Clint Wiseman, USC/Majorana
	August 2015.

	Macro to check a list of run numbers and determine if the built files
	exist.
	Will also check if files have been blinded and aren't readable.

	Usage: 
	root[0] .X CheckFiles.C 
	root[0] .X CheckFiles.C ("InputFile_WithNoExtension")
*/

void CheckFiles(string Input = ""){

	int mode = 1; // switch: 0 for local files, 1 for pdsf files

	// Input a list of run numbers
	if (Input == "") Char_t InputName[200] = "builtVeto_DebugList";
	else Char_t InputName[200] = Input.c_str();
	Char_t InputFile[200];
	sprintf(InputFile,"%s.txt",InputName);
	ifstream InputList;
	InputList.open(InputFile);
	Char_t TheFile[200];

	Int_t run;
	while(!InputList.eof()){

		// initialize 
		InputList >> run;
		if (mode==0) sprintf(TheFile,"~/dev/datasets/builtVeto/OR_run%i.root",run);
		else if (mode==1) sprintf(TheFile,"/global/project/projectdirs/majorana/data/mjd/surfmjd/data/built/P3JDY/OR_run%u.root",run); 

		// if file doesn't exist, ROOT will quit.
	 	TFile f(TheFile);
	 	if (!f || f->IsZombie()) {
       		printf("Cannot open %s!", TheFile);
       		break;
    	}
    	else f.Close();
    }

}
