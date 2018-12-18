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
#include <cstring>

#include "VCFLoader.hpp"
#include "SufficientStatistics.hpp"
#include "DemuxAlgo.hpp"

const int DemuxAlgo::geno_allele_to_prob[3][3] = {{0, 2, 2}, {1, 1, 2}, {2, 0, 2}}

void DemuxAlgo::estimate_background() {
	double *P_ec, *alpha_ec; // expected count arrays
	double *cond_prob;
	int *obs_types;

	P_ec = new double[nDonor];
	alpha_ec = new double[3];
	cond_prob = new double[nDonor];
	obs_types = new int[nDonor];

	std::vector<int>& genovec = ss.cellxgeno[0];
	std::vector<NucDist>& ndvec = ss.cellxnd[0];
	int s = ndvec.size();
	double denom, frac;
	double a, b, c, delta;

	// E step
	memset(P_ec, 0, sizeof(double) * nDonor);
	memset(alpha_ec, 0, sizeof(double) * 3);
	for (int i = 0; i < s; ++i) {
		const SNPType* snp = ss.diff_genos[genovec[i]];
		for (int j = 0; j < NucDist::size; ++j) {
			denom = 0.0;
			for (int k = 0; k < nDonor; ++k) {
				obs_types[k] = geno_allele_to_prob[snp->getDonorGenotype(k)][j];
				cond_prob[k] = P[k] * getObsProb(obs_types[k]);
				denom += cond_prob[k];
			}
			for (int k = 0; k < nDonor; ++k) {
				cond_prob[k] /= denom;
				frac = ndvec[i].dist[j] * cond_prob[k];
				P_ec[k] += frac;
				alpha_ec[obs_types[k]] += frac;
			}
		}
	}

	// M step
	// Estimate P
	denom = 0.0;
	for (int k = 0; k < nDonor; ++k) denom += P_ec[k];
	for (int k = 0; k < nDonor; ++k) P[k] = P_ec[k] / denom;
	for (int k = 0; k < nDonor; ++k) printf("%.5g\t", P[k]); printf("\n");
	// Estimate alpha
	a = 0.375 * (alpha_ec[0] + alpha_ec[1] + alpha_ec[2]);
	b = -0.25 * (3 * alpha_ec[0] + 2 * alpha_ec[1] + 5 * alpha_ec[2]);
	c = alpha_ec[2];
	delta = sqrt(b * b - 4.0 * a * c);
	printf("%.5g %.5g\n", (-b + delta) / a / 2.0, (-b - delta) / a / 2.0);

	delete[] P_ec;
	delete[] alpha_ec;
	delete[] cond_prob;
	delete[] obs_types;
}
