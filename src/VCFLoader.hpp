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

#ifndef VCFLOADER
#define VCFLOADER

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

	void setPos(int pos) { this->pos = pos; }
	void setRef(char ref) { this->ref = ref; }
	void setAlt(char alt) { this->alt = alt; }
	void setDonorGenotype(int donor, std::string genotype);


	static void setNumDonor(int num) { nDonor = num; }

private:
	int pos; // chromosome position
	char ref, alt; // reference allele and alternative allele
	uint64_t* genotypes; // genotypes for each donor, 0: R/R; 1: R/A; 2: A/A; 3: unknown

	static const int SHIFT1 = 4; // one uint64_t stores genotype info for 16 donors
	static const int BASE1 = 15; // % 16 == & 15
	static const int SHIFT2[16]; // shifts within in one uint64_t
	static const uint64_t BASE2 = 3; // within one uint64_t, & 3 to get the genotype

	static int nDonor;
};

typedef std::vector<SNPType> SNPVecType;
typedef std::unordered_map<std::string, SNPVecType> SNPMapType;

class VCFLoader {
public:
	VCFLoader() {}

	void loadVCF(std::string input_vcf_file);

	SNPMapType& getMap() { return snpMap; }

private:
	std::vector<std::string> donor_names;
	SNPMapType snpMap;
};


#endif
