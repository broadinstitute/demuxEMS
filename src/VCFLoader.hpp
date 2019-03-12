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

#ifndef VCFLOADER_H_
#define VCFLOADER_H_

#include <cstdint>
#include <cassert>
#include <string>
#include <vector>
#include <unordered_map>

class SNPType {
public:
	SNPType(int pos, char ref, char alt) : pos(pos), ref(ref), alt(alt), genotypes(nullptr) {
	}

	SNPType(SNPType&& o) noexcept {
		if (this != &o) moveFrom(std::move(o));
	}

	SNPType(const SNPType& o) { // if noexcept copy constructor + noexcept destructor, will use move semantics , but still need copy constructors...
		assert(false);
	}

	SNPType& operator=(SNPType&& o) noexcept {
		if (this != &o) moveFrom(std::move(o));
		return *this;
	}

	SNPType& operator=(const SNPType& o) {
		assert(false);
	}

	~SNPType() noexcept {
		if (genotypes != nullptr) delete[] genotypes;
	}

	void moveFrom(SNPType&& o) {
		pos = o.pos;
		ref = o.ref;
		alt = o.alt;
		genotypes = o.genotypes;
		o.genotypes = nullptr;
	}

	int getPos() const { return pos; }
	char getRef() const { return ref; }
	char getAlt() const { return alt; }
	int getDonorGenotype(int donor) const;

	const uint64_t* getGenotypeVec() const { return genotypes; }

	void setPos(int pos) { this->pos = pos; }
	void setRef(char ref) { this->ref = ref; }
	void setAlt(char alt) { this->alt = alt; }
	void setDonorGenotype(int donor, std::string genotype);


	static void setNumDonor(int num) { 
		nDonor = num;
		geno_size = (nDonor >> SHIFT1) + ((nDonor & BASE1) > 0);
	}

	static int getGenoSize() {
		return geno_size;
	}

private:
	int pos; // chromosome position, 0-based
	char ref, alt; // reference allele and alternative allele
	uint64_t* genotypes; // genotypes for each donor, 0: R/R; 1: R/A; 2: A/A; 3: unknown

	static const int SHIFT1 = 4; // one uint64_t stores genotype info for 16 donors
	static const int BASE1 = 15; // % 16 == & 15
	static const int SHIFT2[16]; // shifts within in one uint64_t
	static const uint64_t BASE2 = 3; // within one uint64_t, & 3 to get the genotype

	static int nDonor, geno_size;
};

typedef std::vector<SNPType> SNPVecType;
typedef std::unordered_map<std::string, SNPVecType> SNPMapType;

class VCFLoader {
public:
	VCFLoader() {
		nDonor = nsnp = nchr = cur_chr = cur_vecp = snp_id = cur_chr_copy = cur_vecp_copy = snp_id_copy = 0;
	}

	void loadVCF(std::string input_vcf_file);

	int getNumDonor() const { return nDonor; }

	const std::vector<std::string>& getDonorNames() const { return donor_names; }
	
	int getTotalSNP() const { return nsnp; }

	void reOrderSNP(const std::vector<std::string>& reorder_names);

	void reset() {
		cur_chr = cur_vecp = snp_id = cur_chr_copy = cur_vecp_copy = snp_id_copy = 0;
	}

	void save_pointer() {
		cur_chr_copy = cur_chr;
		cur_vecp_copy = cur_vecp;
		snp_id_copy = snp_id;
	}

	void load_pointer() {
		cur_chr = cur_chr_copy;
		cur_vecp = cur_vecp_copy;
		snp_id = snp_id_copy;
	}

	bool isValid() const { 
		return snp_id < nsnp;
	}

	int getTid() const {
		return chrom_ids[cur_chr];
	}

	std::string getChromName() const {
		return chrom_names[cur_chr];
	}

	int getPos() const {
		return chrom_snps[cur_chr]->at(cur_vecp).getPos();
	}

	char getRef() const {
		return chrom_snps[cur_chr]->at(cur_vecp).getRef();	
	}

	char getAlt() const {
		return chrom_snps[cur_chr]->at(cur_vecp).getAlt();	
	}

	std::string getDonorGenotype(int donor) const {
	switch(chrom_snps[cur_chr]->at(cur_vecp).getDonorGenotype(donor)) {
		case 0: return "0/0";
		case 1: return "0/1";
		case 2: return "1/1";
		default: assert(false);
	}

	const SNPType& getSNP() const {
		return chrom_snps[cur_chr]->at(cur_vecp);
	}

	int getSnpID() const {
		return snp_id;
	}

	void next() {
		++cur_vecp;
		if (cur_vecp >= (int)chrom_snps[cur_chr]->size()) {
			++cur_chr;
			cur_vecp = 0;
		}
		++snp_id;
	}

	void locate(int target_snp_id) {
		assert(snp_id <= target_snp_id && target_snp_id < nsnp);
		
		int snp_size;

		snp_size = chrom_snps[cur_chr]->size();
		while (snp_id + (snp_size - cur_vecp) <= target_snp_id) {
			snp_id += snp_size - cur_vecp;
			++cur_chr;
			cur_vecp = 0;
			snp_size = chrom_snps[cur_chr]->size();
		}

		cur_vecp += target_snp_id - snp_id;
		snp_id = target_snp_id;
	}

private:
	int nDonor;
	int nsnp;
	std::vector<std::string> donor_names;
	SNPMapType snpMap;

	int nchr; // number of chromosomes that are in the BAM header and also have SNPs.
	std::vector<int> chrom_ids;
	std::vector<std::string> chrom_names;
	std::vector<SNPVecType*> chrom_snps;

	int cur_chr, cur_vecp, snp_id; // current position in chrom_snps (cur_chr) and current position in chrom_snps[cur_chr]
	int cur_chr_copy, cur_vecp_copy, snp_id_copy; // a copy of the cur_chr & cur_vecp;
};

#endif
