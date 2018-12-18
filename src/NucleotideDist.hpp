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

#ifndef NUCLEOTIDEDIST
#define NUCLEOTIDEDIST

#include <cstring>
#include <string>
#include <functional>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include "VCFLoader.hpp"
#include "SufficientStatistics.hpp"
#include "BamAlignment.hpp"



struct BarcodeUMI {
	std::string barcode, umi;

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
			return std::hash<std::string>{}(k.barcode + k.umi);
		}
	};
};

typedef std::unordered_map<BarcodeUMI, NucDist> NucDistMap;

struct GenoKeyType {
	const uint64_t* genotypes;

	GenoKeyType(const uint64_t* genotypes) : genotypes(genotypes) {}

	bool operator< (const GenoKeyType& o) const {
		for (int i = 0; i < SNPType::getGenoSize(); ++i) 
			if (genotypes[i] != o.genotypes[i])
				return genotypes[i] < o.genotypes[i];
		return false;
	}
};

typedef std::unordered_map<int, NucDist*> GenoID2NucDist;

class NucleotideDist {
public:
	NucleotideDist(VCFLoader& vcf_loader, std::string barcode_tag = "CB", std::string umi_tag = "UB", int empty_upper_umi = 50) : vcf_loader(vcf_loader), barcode_tag(barcode_tag), umi_tag(umi_tag), empty_upper_umi(empty_upper_umi) {
		nsnp = vcf_loader.getTotalSNP();
		snp_nucdist_vec = new NucDistMap*[nsnp];
		memset(snp_nucdist_vec, 0, sizeof(NucDistMap*));
		barcodexumi.clear();
		empty_barcodes.clear();
	}

	~NucleotideDist() {
		for (int i = 0; i < nsnp; ++i)
			if (snp_nucdist_vec[i] != nullptr) delete snp_nucdist_vec[i];
		delete[] snp_nucdist_vec;
	}

	int collectData(BamAlignment& ba); // 0, run out of SNPs; 1, valid; 2, skip

	void collectEmptyBarcodes();

	void loadBarcodeList(std::string barcode_list_file);

	void parseData(SufficientStatistics& ss); // parse data and put sufficient statistics into ss

private:
	VCFLoader& vcf_loader;
	std::string barcode_tag, umi_tag;

	int empty_upper_umi;
	std::unordered_map<std::string, std::unordered_set<std::string> > barcodexumi;
	std::unordered_set<std::string> empty_barcodes;

	int nsnp;
	NucDistMap** snp_nucdist_vec;

	std::pair<BarcodeUMI, NucDist> insert_pair;

	std::unordered_map<std::string, int> barcode2bid; // bid barcode_id
	std::vector<std::string> barcode_vec;

	std::map<GenoKeyType, int> geno2gid; // genotype to genotype id

	std::vector<GenoID2NucDist> cells;
};

#endif
