#include <cstdio>
#include <cstdint>
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
#include "BamAlignment.hpp"
#include "DataCollector.hpp"

using namespace std;

int num_threads;
htsThreadPool p = {NULL, 0};

VCFLoader vcf_loader;
SamParser *parser;
BamAlignment ba(true); // important to set true here, otherwise will read paired-end reads
DataCollector *data_collector;


int main(int argc, char* argv[]) {
	if (argc != 4 && argc != 6) {
		printf("Usage: scan_snp input.vcf.gz input.bam output_name [-p number_of_threads]\n");
		exit(-1);
	}

	if (argc == 6) {
		assert(!strcmp(argv[4], "-p"));
		num_threads = atoi(argv[5]);
	}

	vcf_loader.loadVCF(argv[1]);

	// num_threads = 4;
	if (num_threads > 1) assert(p.pool = hts_tpool_init(num_threads));
	parser = new SamParser(argv[2], num_threads > 1 ? &p : NULL);

	const bam_hdr_t* header = parser->getHeader();
	vector<string> chr_names;
	for (int i = 0; i < header->n_targets; ++i) {
		string chr_name = string(header->target_name[i]);
		if (chr_name.substr(0, 3) == "chr" && chr_name != "chrM") chr_names.push_back(chr_name.substr(3));
		else chr_names.push_back(std::move(chr_name));
	}

	vcf_loader.reOrderSNP(chr_names);

	int cnt = 0, nalign = 0;
	data_collector = new DataCollector(vcf_loader);
	while (ba.read(parser)) {
		if (ba.isAligned() & 1) {
			int ret = data_collector->collectData(ba);
			if (ret == 0) break;
			nalign += ret == 1;
		}
		++cnt;
		if (cnt % 1000000 == 0) printf("PROCESSED %d lines, %d aligned.\n", cnt, nalign);
	}
	printf("BAM file is parsed!\n");

	data_collector->outputDataMatrix(argv[3]);
	printf("Output files are generated!\n");
	
	delete data_collector;
	delete parser;
	if (num_threads > 1) hts_tpool_destroy(p.pool);

	return 0;
}
