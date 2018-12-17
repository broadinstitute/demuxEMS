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

#include <cmath>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <fstream>

#include "VCFLoader.hpp"
#include "SufficientStatistics.hpp"
#include "CIGARstring.hpp"
#include "SEQstring.hpp"
#include "QUALstring.hpp"
#include "BamAlignment.hpp"
#include "NucleotideDist.hpp"



bool NucleotideDist::collectData(BamAlignment& ba) {
	int tid, pos, rpos; // rpos, read pos
	int optype, oplen;
	int ub, snp_rpos; // ub, upper bound

	CIGARstring ci;
	SEQstring si;
	QUALstring qi;

	uint8_t* tag_p;
	char tag_type;

	char base;
	float p_err, p_crt;

	int snp_id;



	tid = ba.get_tid();
	pos = ba.getLeftMostPos();

	assert(ba.findTag(barcode_tag.c_str(), tag_p, tag_type) && tag_type == 'Z');
	insert_pair.first.barcode = std::move(std::string(ba.tag2Z(tag_p)));
	assert(ba.findTag(umi_tag.c_str(), tag_p, tag_type) && tag_type == 'Z');
	insert_pair.first.umi = std::move(std::string(ba.tag2Z(tag_p)));

	while (vcf_loader.isValid() && ((vcf_loader.getTid() < tid) || ((vcf_loader.getTid() == tid) && (vcf_loader.getPos() < pos)))) vcf_loader.next();
	if (!vcf_loader.isValid()) return false;

	if (vcf_loader.getTid() == tid) {
		vcf_loader.save_pointer();

		ba.getCIGAR(ci); ci.setCurrent();
		ba.getSEQ(si); si.setCurrent();
		ba.getQUAL(qi); qi.setCurrent();

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

					snp_id = vcf_loader.getSnpID();
					if (snp_nucdist_vec[snp_id] == nullptr) {
						snp_nucdist_vec[snp_id] = new NucDistMap();
					}
					auto result = snp_nucdist_vec[snp_id].insert(insert_pair);

					p_err = pow(10.0, -qi.qualAt(snp_rpos) / 10.0);
					p_crt = 1.0 - p_err;
					p_err /= 3.0;

					base = si.baseAt(snp_rpos);

					if (base == snp.getRef()) {
						result.first->second[NucDist::R] += p_crt;
						result.first->second[NucDist::A] += p_err;
						result.first->second[NucDist::O] += p_err * 2;
					}
					else if (base == snp.getAlt()) {
						result.first->second[NucDist::R] += p_err;
						result.first->second[NucDist::A] += p_crt;
						result.first->second[NucDist::O] += p_err * 2;
					}
					else {
						result.first->second[NucDist::R] += p_err;
						result.first->second[NucDist::A] += p_err;
						result.first->second[NucDist::O] += p_crt + p_err;
					}

					vcf_loader.next();
				}
			}

			rpos += (optype & 1) ? oplen : 0;
			pos += (optype & 2) ? oplen : 0;
		}
		
		vcf_loader.load_pointer();
	}

	return true;
}

void NucleotideDist::loadBarcodeList(std::string barcode_list_file) {
	int bid;
	std::ifstream fin(barcode_list_file);
	std::string line;

	barcode2bid.clear();
	barcode_vec.resize(1, ""); // 0 is "" for ambient RNA

	bid = 0;
	while (std::getline(fin, line)) {
		barcode2bid[line] = ++bid;
		barcode_vec.push_pack(line);
	}

	fin.close();
	printf("Barcode list is loaded!\n");
}

void NucleotideDist::parseData(SufficientStatistics& ss) {
	int gid = 0;
	int cur_gid, barcode_id;
	GenoID2NucDist example;

	geno2gid.clear();
	example.clear();
	cells.resize(barcode_vec.size());


	ss.cell_barcodes = barcode_vec;
	ss.diff_genos.clear();


	vcf_loader.reset();
	for (int i = 0; i < nsnp; ++i) {
		if (snp_nucdist_vec[i] != nullptr) {
			const SNPType& snp = vcf_loader.getSNP();
			auto geno2gid_result = geno2gid.insert(std::make_pair(GenoKeyType(snp.getGenotypeVec()), gid));
			cur_gid = geno2gid_result.first->second;

			if (geno2gid_result.second) {
				ss.diff_genos.push_back(&snp);
				++gid;
			}

			for (auto&& kv : *snp_nucdist_vec[i]) {
				kv.second.normalize();
				auto iter = barcode2bid.find(kv.first.barcode);
				barcode_id = iter == barcode2bid.end() ? 0 : iter->second;
				auto result = cells[barcode_id].insert(std::make_pair(cur_gid, &kv.second));
				if (!result.second) {
					result.first->second->add(kv.second);
				}
			}
		}
		vcf_loader.next();
	}

	ss.ncells = ss.cell_barcodes.size() - 1;
	ss.ngenos = ss.diff_genos.size();

	ss.cellxgeno.resize(ss.ncells + 1);
	ss.cellxnd.resize(ss.ncells + 1);

	for (int i = 0; i <= ss.ncells; ++i) {
		ss.cellxgeno[i].resize(cells[i].size());
		ss.cellxnd[i].resize(cells[i].size());
		int j = 0;
		for (auto&& kv : cells[i]) {
			ss.cellxgeno[i][j] = kv.first;
			ss.cellxnd[i][j].copy(*kv.second);
		}
	}
}
