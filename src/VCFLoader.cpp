/* Copyright (c) 2018
   Bo Li (The Broad Institute of MIT and Harvard)
   libo@broadinstitute.org

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.   

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA
*/

#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "VCFLoader.hpp"



const int SNPType::SHIFT2[16] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};

int SNPType::nDonor;

int SNPType::getDonorGenotype(int donor) const {
	assert(nDonor > donor && genotypes != nullptr);
	return (genotypes[donor >> SHIFT1] >> SHIFT2[donor & BASE1]) & BASE2;
}

void SNPType::setDonorGenotype(int donor, std::string genotype) {
	if (genotypes == nullptr) {
		int s = (nDonor >> SHIFT1) + ((nDonor & BASE1) > 0);
		genotypes = new uint64_t[s];
		for (int i = 0; i < s; ++i) genotypes[i] = -1;
	}
   
	uint64_t gvalue = BASE2;

	if (genotype == "0/0") gvalue = 0;
	else if (genotype == "0/1") gvalue = 1;
	else if (genotype == "1/1") gvalue = 2;

	int pos1 = donor >> SHIFT1;
	int pos2 = donor & BASE1;

	genotypes[pos1] = genotypes[pos1] - (((genotypes[pos1] >> SHIFT2[pos2]) & BASE2) << SHIFT2[pos2]) + (gvalue << SHIFT2[pos2]);
}



void VCFLoader::loadVCF(std::string input_vcf_file) {
	int nDonor;
	std::ifstream fin;
	boost::iostreams::filtering_istream gin;
	std::string line, field;
	std::istringstream strin;

	bool is_gzip = input_vcf_file.substr(input_vcf_file.length() - 3, 3) == ".gz";

	if (is_gzip) {
		fin.open(input_vcf_file, std::ios_base::in | std::ios_base::binary);
		gin.push(boost::iostreams::gzip_decompressor());
		gin.push(fin);
	} 
	else {
		fin.open(input_vcf_file);
		gin.push(fin);
	}

	while (std::getline(gin, line) && line[0] == '#' && line[1] == '#');

	assert(line[0] == '#');
	strin.clear();
	strin.str(line.substr(1));

	std::string fixed[] = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"};
	for (int i = 0; i < 8; ++i)
		assert(std::getline(strin, field, '\t') && field == fixed[i]);
	if (std::getline(strin, field, '\t')) {
		assert(field == "FORMAT");
		donor_names.clear();
		while (std::getline(strin, field, '\t')) {
			donor_names.push_back(field);
		}
	}

	nDonor = donor_names.size();
	SNPType::setNumDonor(nDonor);
	snpMap.clear();

	std::string curr_chrom = "", chrom;
	int pos;
	char ref, alt;
	SNPVecType *snp_vec = nullptr, new_vec;

	new_vec.clear();
	nsnp = 0;
	int cnt = 0;
	while (std::getline(gin, line)) {
		++cnt;

		strin.clear();
		strin.str(line);
		assert(std::getline(strin, chrom, '\t')); // CHROM

		assert(std::getline(strin, field, '\t')); // POS
		pos = std::stoi(field) - 1; // convert to 0-based position

		assert(std::getline(strin, field, '\t')); // ID

		assert(std::getline(strin, field, '\t')); // REF
		if (field.length() != 1 || !(field[0] == 'A' || field[0] == 'C' || field[0] == 'G' || field[0] == 'T')) continue;
		ref = field[0];

		assert(std::getline(strin, field, '\t')); // ALT
		if (field.length() != 1 || !(field[0] == 'A' || field[0] == 'C' || field[0] == 'G' || field[0] == 'T')) continue;
		alt = field[0];

		assert(std::getline(strin, field, '\t')); // QUAL
		
		assert(std::getline(strin, field, '\t')); // FILTER
		if (field != "PASS") continue;

		assert(std::getline(strin, field, '\t')); // INFO

		if (chrom != curr_chrom) {
			auto iter = snpMap.find(chrom);
			if (iter == snpMap.end()) {
				snpMap[chrom] = new_vec;
			}
			curr_chrom = chrom;
			snp_vec = &snpMap[chrom];
		}

		snp_vec->emplace_back(pos, ref, alt);
		++nsnp;

		if (nDonor > 0) {
			assert(std::getline(strin, field, '\t') && field.substr(0, 2) == "GT"); // FORMAT
			for (int i = 0; i < nDonor; ++i) {
				assert(std::getline(strin, field, '\t')); // Donor
				snp_vec->back().setDonorGenotype(i, field.substr(0, field.find_first_of(':')));
			}
		}
				
		if (cnt % 1000000 == 0) printf("Processed %d lines, loaded %d SNPs.\n", cnt, nsnp);
	}

	printf("Loaded %d SNPs in total.\n", nsnp);

	fin.close();
	if (is_gzip) gin.reset();
}

void VCFLoader::reOrderSNP(std::vector<std::string>& chrom_names) {
	chrom_ids.clear();
	chrom_snps.clear();

	int tid = 0;
	for (auto&& chrom_name : chrom_names) {
		auto iter = snpMap.find(chrom_name);
		if (iter != snpMap.end()) {
			chrom_ids.push_back(tid);
			chrom_snps.push_back(&(iter->second));
		}
		++tid;
	}

	nchr = chrom_ids.size();
	reset();
}
