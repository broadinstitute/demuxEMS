#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "VCFLoader.hpp"
#include "SamParser.hpp"
#include "CIGARstring.hpp"
#include "SEQstring.hpp"
#include "QUALstring.hpp"
#include "BamAlignment.hpp"


using namespace std;

int num_threads;
htsThreadPool p = {NULL, 0};

VCFLoader vcf_loader;
SamParser *parser;
BamAlignment ba(true); // important to set true here, otherwise will read paired-end reads

int tid, pos, rpos; // rpos, read pos
int optype, oplen;
int ub, snp_rpos; // ub, upper bound

CIGARstring ci;
SEQstring si;
QUALstring qi;

uint8_t* tag_p;
char tag_type;

unordered_map<string, unordered_set<int> > cells;
pair<string, unordered_set<int> > insert_pair;
unordered_set<int>* snp_set;

int main(int argc, char* argv[]) {
	if (argc != 4) {
		printf("Usage: scan_snp input.vcf.gz input.bam output.txt\n");
		exit(-1);
	}

	vcf_loader.loadVCF(argv[1]);

	num_threads = 4;
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

	cells.clear();
	insert_pair.second.clear();

	int cnt = 0;
	while (ba.read(parser)) {
		if (ba.isAligned() & 1) {
			tid = ba.get_tid();
			pos = ba.getLeftMostPos();

			while (vcf_loader.isValid() && ((vcf_loader.getTid() < tid) || ((vcf_loader.getTid() == tid) && (vcf_loader.getPos() < pos)))) vcf_loader.next();
			if (!vcf_loader.isValid()) break;

			if (vcf_loader.getTid() == tid) {
				vcf_loader.save_pointer();

				ba.getCIGAR(ci); ci.setCurrent();
				ba.getSEQ(si); si.setCurrent();
				ba.getQUAL(qi); qi.setCurrent();
				tag_p = nullptr;

				rpos = 0;
				for (int i = 0; i < ci.getLen(); ++i) {
					oplen = ci.oplenAt(i);
					optype = ci.optypeAt(i);

					if (ci.opchrAt(i) == 'M') {
						while (vcf_loader.getPos() < pos) {
							vcf_loader.next();
							if (!vcf_loader.isValid() || vcf_loader.getTid() > tid) break;
						}
						ub = pos + oplen;
						while (vcf_loader.isValid() && vcf_loader.getTid() == tid && vcf_loader.getPos() < ub) {
							snp_rpos = rpos + (vcf_loader.getPos() - pos);
							const SNPType& snp = vcf_loader.getSNP();
							if (si.baseAt(snp_rpos) == snp.getAlt() && qi.qualAt(snp_rpos) >= 30) {
								if (tag_p == nullptr) {
									ba.findTag("RG", tag_p, tag_type);
									assert(tag_type == 'Z');
									insert_pair.first = string(ba.tag2Z(tag_p));
									auto result = cells.insert(insert_pair);
									snp_set = &(result.first->second);
								}
								snp_set->insert(vcf_loader.getSnpID());
							}
							vcf_loader.next();
						}
					}

					rpos += (optype & 1) ? oplen : 0;
					pos += (optype & 2) ? oplen : 0;
				}
				
				vcf_loader.load_pointer();
			}
		}

		++cnt;
		if (cnt % 1000000 == 0) printf("PROCESSED %d lines, seen %d cells.\n", cnt, (int)cells.size());
	}

	ofstream fout(argv[3]);
	fout<< "BARCODE\tNSNP"<< endl;
	for (auto&& apair : cells) {
		fout<< apair.first << '\t'<< apair.second.size() << endl;
	}
	fout.close();

	delete parser;
	if (num_threads > 1) hts_tpool_destroy(p.pool);

	return 0;
}
