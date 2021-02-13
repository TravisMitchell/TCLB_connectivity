#include <stdio.h>
#include "Preprocessor.h"
#include <iostream>
#include <cstdlib>
#include "string.h"




int main(int argc, char** argv) {
    // check we have enough parameters
    if(argc != 10) {
		printf("Format: connectivity {STL_file} {nx} {ny} {nz} {x0} {y0} {z0} {d} {Q}\n");
		return 1;
	}

    Preprocessor preproc;
    
    // parse parameters
    char* stl_filename = argv[1];
    preproc.nx = atoi(argv[2]);
    preproc.ny = atoi(argv[3]);
    preproc.nz = atoi(argv[4]);
    preproc.x0 = atoi(argv[5]);
    preproc.y0 = atoi(argv[6]);
    preproc.z0 = atoi(argv[7]);
    preproc.d_ = atoi(argv[8]);
    preproc.Q_ = atoi(argv[9]);
    
    char* stlfn = (char*) malloc(strlen(stl_filename) + 4 + 1);
    strcpy(stlfn, stl_filename);
    strcat(stlfn, ".stl");

    preproc.loadSTL(stlfn);
    preproc.generateKDTree();
    preproc.generateDomain();
    preproc.writeToConnectivity(stl_filename);
}