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

#ifndef DATACOLLECTOR_H_
#define DATACOLLECTOR_H_

#include <cstdint>
#include <cstring>
#include <string>
#include <functional>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <boost/functional/hash.hpp>

#include "VCFLoader.hpp"
#include "SufficientStatistics.hpp"
#include "BamAlignment.hpp"



struct BarcodeUMI {
	uint64_t barcode, umi;

	bool operator==(const BarcodeUMI &o) const {
		return barcode == o.barcode && umi == o.umi;
	}
};

namespace std
{
	template<> struct hash<BarcodeUMI>
	{
		typedef BarcodeUMI argument_type;
		typedef std::size_t result_type;
		result_type operator()(argument_type const& k) const noexcept
		{
			result_type seed = 0;
			boost::hash_combine(seed, k.barcode);
			boost::hash_combine(seed, k.umi);
			return seed;
		}
	};
};


typedef std::unordered_map<BarcodeUMI, NucDist> NucDistMap;

typedef std::vector<NucDist> NucDistVec;
typedef std::unordered_map<int, NucDistVec> SNP2NucDist;
typedef std::unordered_map<uint64_t, SNP2NucDist> Barcode2SNP;



// struct GenoKeyType {
// 	const uint64_t* genotypes;

// 	GenoKeyType(const uint64_t* genotypes) : genotypes(genotypes) {}

// 	bool operator< (const GenoKeyType& o) const {
// 		for (int i = 0; i < SNPType::getGenoSize(); ++i) 
// 			if (genotypes[i] != o.genotypes[i])
// 				return genotypes[i] < o.genotypes[i];
// 		return false;
// 	}
// };

// typedef std::unordered_map<int, NucDist*> GenoID2NucDist;

class DataCollector {
public:
	DataCollector(VCFLoader& vcf_loader, std::string barcode_tag = "CB", std::string umi_tag = "UB") : vcf_loader(vcf_loader), barcode_tag(barcode_tag), umi_tag(umi_tag) {
		barcode_len = -1;

		max_snpid = 0;
		snp_nucdist_vec = new NucDistMap*[vcf_loader.getTotalSNP()];
		memset(snp_nucdist_vec, 0, sizeof(NucDistMap*));
		
		nobs = 0;
		snpid_vec.clear();
		data_matrix.clear();

		insert_pair_bu2n.second.clear();
		insert_pair_b2s.second.clear();
		insert_pair_s2n.second.clear();

		// barcodexumi.clear();
		// empty_barcodes.clear();
		// insert_pair_b2u.second.clear();
	}

	~DataCollector() {
		delete[] snp_nucdist_vec;
	}

	int collectData(BamAlignment& ba); // 0, run out of SNPs; 1, valid; 2, skip

	// process SNPs that are loaded into snp_nucdist_vec but not processed into data_matrix
	void wrapUp() {
		while (vcf_loader.getSnpID() < max_snpid) {
			transferToDataMatrix(vcf_loader.getSnpID());
			vcf_loader.next();
		}
	}

	void outputDataMatrix(const std::string& output_name);

	// void collectEmptyBarcodes();

	// void loadBarcodeList(std::string barcode_list_file);

	// void parseData(SufficientStatistics& ss); // parse data and put sufficient statistics into ss

private:
	VCFLoader& vcf_loader;
	std::string barcode_tag, umi_tag;
	int barcode_len;

	int max_snpid; // max_snpid, max_snpid - 1 is the maximum snpID ever explored.
	NucDistMap** snp_nucdist_vec;
	std::pair<BarcodeUMI, NucDist> insert_pair_bu2n;

	int nobs;
	std::vector<int> snpid_vec; // observed snpids;

	Barcode2SNP data_matrix;
	std::pair<uint64_t, SNP2NucDist> insert_pair_b2s;
	std::pair<int, NucDistVec> insert_pair_s2n;


	// int empty_upper_umi;
	// std::unordered_map<uint64_t, std::unordered_set<uint64_t> > barcodexumi; // this method may underestimate # of UMIs 
	// std::pair<uint64_t, std::unordered_set<uint64_t> > insert_pair_b2u;
	// std::unordered_set<uint64_t> empty_barcodes;


	// std::unordered_map<std::string, int> barcode2bid; // bid barcode_id
	// std::vector<std::string> barcode_vec;

	// std::map<GenoKeyType, int> geno2gid; // genotype to genotype id

	// std::vector<GenoID2NucDist> cells;

	void transferToDataMatrix(int snp_id) {
		if (snp_nucdist_vec[snp_id] != nullptr) {
			for (auto& kv : *snp_nucdist_vec[snp_id]) {
				kv.second.normalize();
				insert_pair_b2s.first = kv.first.barcode;
				auto b2s_res = data_matrix.insert(insert_pair_b2s);
				SNP2NucDist& s2n = b2s_res.first->second;
				insert_pair_s2n.first = nobs;
				auto s2n_res = s2n.insert(insert_pair_s2n);
				s2n_res.first->second.push_back(std::move(kv.second));
			}
			++nobs;
			snpid_vec.push_back(snp_id);
			delete snp_nucdist_vec[snp_id];
			snp_nucdist_vec[snp_id] = nullptr;
		}
	}
};

#endif
