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

#ifndef DEMUXALGO_H_
#define DEMUXALGO_H_

#include <string>
#include <vector>

#include "SufficientStatistics.hpp"


class DemuxAlgo {
public:
	DemuxAlgo(int nDonor, SufficientStatistics& ss) : nDonor(nDonor), ss(ss), alpha(0.05) { 
		double init_prob;

		P = new double[nDonor];
		init_prob = 1.0 / nDonor;
		for (int i = 0; i < nDonor; ++i) P[i] = init_prob;

		theta = new double*[ss.ncells + 1];
		init_prob = 1.0 / (nDonor + 1);
		for (int i = 1; i <= ss.ncells; ++i) {
			theta[i] = new double[nDonor + 1];
			for (int j = 0; j <= nDonor; ++j) theta[i][j] = init_prob;
		}
	}

	~DemuxAlgo() {
		delete[] P;
		for (int i = 1; i <= ss.ncells; ++i) delete[] theta[i];
		delete[] theta;
	}

	void estimate_background(double tol = 1e-6);

	void demultiplex(int thread_id, int fr, int to, double prior_noise, double prior_donor, double tol);

	void demuxEMS(int num_threads, double prior_noise = 1.0, double prior_donor = 0.0, double tol = 1e-6);

	void writeOutputs(std::string output_name, const std::vector<std::string>& donor_names);

private:
	int nDonor;
	SufficientStatistics& ss;

	double alpha, *P, **theta;

	static const int geno_allele_to_prob[3][3]; // geno: 0, 0/0; 1, 0/1; 2, 1/1; allele: 0, R; 1, A; 2, O; prob: 0, 1 - 0.75 alpha; 1, 0.5 - 0.25 alpha; 2, 0.25 alpha

	double getObsProb(int obs_type) {
		switch(obs_type) {
			case 0: return 1.0 - 0.75 * alpha;
			case 1: return 0.5 - 0.25 * alpha;
			case 2: return 0.25 * alpha;
			default: assert(false);
		}
	}
};

#endif
