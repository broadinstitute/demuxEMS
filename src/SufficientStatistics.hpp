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

#ifndef SUFFICIENTSTATISTICS_H_
#define SUFFICIENTSTATISTICS_H_

#include <cstring>
#include <string>
#include <vector>

#include "VCFLoader.hpp"

struct NucDist {
	static const int size = 3;
	static const int R = 0; // Ref
	static const int A = 1; // Alt
	static const int O = 2; // Other

	double dist[size]; // 0, Ref; 1, Alt; 2, Other.

	NucDist() { memset(dist, 0, sizeof(dist)); }

	void normalize() {
		double denom = 0.0;
		for (int i = 0; i < size; ++i) denom += dist[i];
		for (int i = 0; i < size; ++i) dist[i] /= denom;
	}

	void add(const NucDist& o) {
		for (int i = 0; i < size; ++i) dist[i] += o.dist[i];
	}

	void copy(const NucDist& o) {
		for (int i = 0; i < size; ++i) dist[i] = o.dist[i];
	}
};

struct SufficientStatistics {
	int ncells; // number of cell barcodes
	std::vector<std::string> cell_barcodes; // cell_barcodes[0] is ambient RNA
	std::vector<const SNPType*> diff_genos; // distinct genotypes

	std::vector<std::vector<int> > cellxgeno; // cell x geno matrix
	std::vector<std::vector<NucDist> > cellxnd; // cell x nuc dist matrix 

	SufficientStatistics() {
		ncells = 0;
		cell_barcodes.clear();
		diff_genos.clear();
		cellxgeno.clear();
		cellxnd.clear();
	}
};

#endif
