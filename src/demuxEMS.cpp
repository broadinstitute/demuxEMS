#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <string>
#include <vector>
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
double prior_noise, prior_donor, tol;

htsThreadPool p = {NULL, 0};

VCFLoader vcf_loader;
SamParser *parser;
BamAlignment ba(true); // important to set true here, otherwise will read paired-end reads

SufficientStatistics ss;
NucleotideDist *nd;
DemuxAlgo *algo;


int main(int argc, char* argv[]) {
	if (argc < 5) {
		printf("Usage: demuxEMS input.vcf.gz input.bam barcodes.tsv output_name [-p num_threads] [--prior-noise prior_noise] [--prior-donor prior_donor] [--tol tol]\n");
		exit(-1);
	}

	num_threads = 1;
	prior_noise = 1.0;
	prior_donor = 0.0;
	tol = 1e-6;
	for (int i = 5; i < argc; ++i) {
		if (!strcmp(argv[i], "-p")) num_threads = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "--prior-noise")) prior_noise = atof(argv[i + 1]);
		if (!strcmp(argv[i], "--prior-donor")) prior_donor = atof(argv[i + 1]);
		if (!strcmp(argv[i], "--tol")) tol = atof(argv[i + 1]);
	}

	vcf_loader.loadVCF(argv[1]);

	if (num_threads > 1) assert(p.pool = hts_tpool_init(num_threads));
	parser = new SamParser(argv[2], num_threads > 1 ? &p : NULL);

	const bam_hdr_t* header = parser->getHeader();
	vector<string> chr_names;
	for (int i = 0; i < header->n_targets; ++i)
		chr_names.emplace_back(header->target_name[i]);
	vcf_loader.reOrderSNP(chr_names);

	while (vcf_loader.isValid()) {
		const SNPType& snp = vcf_loader.getSNP();
		for (int i = 0; i < vcf_loader.getNumDonor(); ++i) {
			int geno = snp.getDonorGenotype(i);
			if (geno < 0 || geno > 2) {
				printf("snp_id = %d, tid = %d, pos = %d, donor_id = %d, geno = %d\n", vcf_loader.getSnpID(), vcf_loader.getTid(), vcf_loader.getPos(), i, geno);
				exit(-1);
			}
		}
		vcf_loader.next();
	}
	vcf_loader.reset();

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
	algo->demuxEMS(num_threads);
	algo->writeOutputs(argv[4], vcf_loader.getDonorNames());
	delete algo;

	return 0;
}
