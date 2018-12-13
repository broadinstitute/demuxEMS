#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "VCFLoader.hpp"
#include "SamParser.hpp"

using namespace std;

int num_threads;
htsThreadPool p = {NULL, 0};

VCFLoader vcf_loader;
SamParser *parser;

int main(int argc, char* argv[]) {
	if (argc != 3) {
		printf("Usage: scan_snp input.vcf.gz input.bam\n");
		exit(-1);
	}

	num_threads = 4;
	if (num_threads > 1) assert(p.pool = hts_tpool_init(num_threads));
	parser = new SamParser(argv[2], num_threads > 1 ? &p : NULL);

	const bam_hdr_t* header = parser->getHeader();
	for (int i = 0; i < header->n_targets; ++i) {
		string chr_name = string(header->target_name[i]);
		if (chr_name.substr(0, 3) == "chr" && chr_name != "chrM") cout<< chr_name<< endl;

	}
	// vcf_loader.loadVCF(argv[1]);

	// SNPMapType& snp_map = vcf_loader.getMap();

	delete parser;
	
	return 0;
}
