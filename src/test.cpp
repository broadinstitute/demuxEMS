#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "VCFLoader.hpp"

using namespace std;

VCFLoader vcf_loader;

int main(int argc, char* argv[]) {
	if (argc != 2) {
		printf("Usage: test input.vcf\n");
		exit(-1);
	}

	vcf_loader.loadVCF(argv[1]);

	SNPMapType& snp_map = vcf_loader.getMap();

	// cout<< snp_map["1"][0].getDonorGenotype(4)<< endl;
	cout<< snp_map["X"][0].getPos()<< endl;
	
	return 0;
}
