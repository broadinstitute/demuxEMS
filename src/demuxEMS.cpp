#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <string>
#include <fstream>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "VCFLoader.hpp"
#include "SufficientStatistics.hpp"
#include "SamParser.hpp"
#include "BamAlignment.hpp"
#include "NucleotideDist.hpp"
#include "DemuxAlgo.hpp"

using namespace std;



int num_threads;
htsThreadPool p = {NULL, 0};

VCFLoader vcf_loader;
SamParser *parser;
BamAlignment ba(true); // important to set true here, otherwise will read paired-end reads

SufficientStatistics ss;
NucleotideDist *nd;
DemuxAlgo *algo;


int main(int argc, char* argv[]) {
	if (argc != 5) {
		printf("Usage: demuxEMS input.vcf.gz input.bam barcodes.tsv output.txt\n");
		exit(-1);
	}

	vcf_loader.loadVCF(argv[1]);

	num_threads = 4;
	if (num_threads > 1) assert(p.pool = hts_tpool_init(num_threads));
	parser = new SamParser(argv[2], num_threads > 1 ? &p : NULL);

	const bam_hdr_t* header = parser->getHeader();
	vector<string> chr_names;
	for (int i = 0; i < header->n_targets; ++i)
		chr_names.emplace_back(header->target_name[i]);
	vcf_loader.reOrderSNP(chr_names);

	int cnt = 0, nalign = 0;
	nd = new NucleotideDist(vcf_loader);
	while (ba.read(parser)) {
		if (ba.isAligned() & 1) {
			int ret = nd->collectData(ba);
			if (ret == 0) break;
			nalign += ret == 1;
		}
		++cnt;
		if (cnt % 1000000 == 0) printf("PROCESSED %d lines, %d aligned.\n", cnt, nalign);
	}

	delete parser;
	if (num_threads > 1) hts_tpool_destroy(p.pool);

	nd->collectEmptyBarcodes();

	nd->loadBarcodeList(argv[3]);
	nd->parseData(ss);
	delete nd;

	algo = new DemuxAlgo(vcf_loader.getNumDonor(), ss);
	algo->estimate_background();
	delete algo;
	// ofstream fout(argv[4]);
	// fout<< ss.cellxgeno[0].size()<< endl;
	// for (int i = 0; i < (int)ss.cellxgeno[0].size(); ++i) {
	// 	fout<< ss.cellxgeno[0][i];
	// 	for (int j = 0; j < NucDist::size; ++j) fout<< '\t'<< ss.cellxnd[0][i].dist[j];
	// 	fout<< endl;
	// }
	// fout.close();

	return 0;
}
