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
#include <thread>
#include <vector>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "VCFLoader.hpp"
#include "SufficientStatistics.hpp"
#include "DemuxAlgo.hpp"

const int DemuxAlgo::geno_allele_to_prob[3][3] = {{0, 2, 2}, {1, 1, 2}, {2, 0, 2}};

void DemuxAlgo::estimate_background(double tol) {
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

	double alpha_old, epsilon;
	double *P_old = new double[nDonor];

	int cnt = 0;

	do {
		alpha_old = alpha;
		memcpy(P_old, P, sizeof(double) * nDonor);

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
		// Estimate alpha
		a = 0.375 * (alpha_ec[0] + alpha_ec[1] + alpha_ec[2]);
		b = -0.25 * (3 * alpha_ec[0] + 2 * alpha_ec[1] + 5 * alpha_ec[2]);
		c = alpha_ec[2];
		delta = sqrt(b * b - 4.0 * a * c);
		alpha = (-b - delta) / a / 2.0;

		// Calculate epsilon
		epsilon = fabs(alpha - alpha_old);
		for (int k = 0; k < nDonor; ++k) epsilon += fabs(P[k] - P_old[k]);

		++cnt;
	} while (epsilon > tol);

	delete[] P_ec;
	delete[] alpha_ec;
	delete[] cond_prob;
	delete[] obs_types;

	printf("Background estimated. Total iteration %d, alpha = %.5g, P =", cnt, alpha);
	for (int k = 0; k < nDonor; ++k) printf(" %.5g", P[k]);
	printf("\n");
}

void DemuxAlgo::demultiplex(int thread_id, int fr, int to, double prior_noise, double prior_donor, double tol) {
	int s;
	double epsilon, denom, obs_prob;
	double *cond_prob, *theta_vec, *theta_vec_ec, *theta_vec_old;

	cond_prob = new double[nDonor + 1];
	theta_vec_ec = new double[nDonor + 1];
	theta_vec_old = new double[nDonor + 1];

	for (cell_id = fr; cell_id < to; ++cell_id) {
		std::vector<int>& genovec = ss.cellxgeno[cell_id];
		std::vector<NucDist>& ndvec = ss.cellxnd[cell_id];
		s = ndvec.size();
		theta_vec = theta[cell_id];

		do {
			memcpy(theta_vec_old, theta_vec, sizeof(double) * (nDonor + 1));

			// E step
			memset(theta_vec_ec, 0, sizeof(double) * (nDonor + 1));
			for (int i = 0; i < s; ++i) {
				const SNPType* snp = ss.diff_genos[genovec[i]];
				for (int j = 0; j < NucDist::size; ++j) {
					denom = 0.0;
					cond_prob[nDonor] = 0.0;
					for (int k = 0; k < nDonor; ++k) {
						obs_prob = getObsProb(geno_allele_to_prob[snp->getDonorGenotype(k)][j]);
						cond_prob[k] = theta_vec[k] * obs_prob;
						cond_prob[nDonor] += theta_vec[nDonor] * P[k] * obs_prob;
						denom += cond_prob[k];
					}
					denom += cond_prob[nDonor];
					for (int k = 0; k <= nDonor; ++k) {
						cond_prob[k] /= denom;
						theta_vec_ec[k] += ndvec[i].dist[j] * cond_prob[k];
					}
				}
			}

			// M step
			denom = 0.0;
			for (int k = 0; k < nDonor; ++k) denom += std::max(theta_vec_ec[k] + (prior_donor - 1.0), 0.0);
			denom += std::max(theta_vec_ec[nDonor] + (prior_noise - 1.0), 0.0);
			for (int k = 0; k <= nDonor; ++k) theta_vec[k] = theta_vec_ec[k] / denom;

			// Calculate epsilon
			epsilon = 0.0;
			for (int k = 0; k <= nDonor; ++k) epsilon += fabs(theta_vec[k] - theta_vec_old[k]);

		} while (epsilon > tol);
	}

	delete[] cond_prob;
	delete[] theta_vec_ec;
	delete[] theta_vec_old;

	printf("DemuxAlgo::demultiplex thread %d is finished!\n", thread_id);
}

void DemuxAlgo::demuxEMS(int num_threads, double prior_noise, double prior_donor, double tol) {
	int fr, to;
	int base = ss.ncells / num_threads;
	int left = ss.ncells % num_threads;
	std::vector<std::thread> threads;

	fr = 1;
	threads.clear();
	for (int i = 0; i < num_threads; ++i) {
		to = fr + (base + (i < left));
		threads.push_back(std::thread(&DemuxAlgo::demultiplex, this, i, fr, to, prior_noise, prior_donor, tol));
		fr = to;
	}
	for (int i = 0; i < (int)threads.size(); ++i)
		threads[i].join();
	printf("DemuxAlgo::demuxEMS is finished.\n");
}

void DemuxAlgo::writeOutputs(std::string output_name, const std::vector<std::string>& donor_names) {
	double nsnp;
	std::ofstream fout;

	fout.open(output_name + ".params");
	fout.setf(std::ios::scientific, std::ios::floatfield);
	fout.precision(5);
	fout<< alpha<< std::endl;
	fout<< P[0];
	for (int i = 1; i < nDonor; ++i) fout<< '\t'<< P[i];
	fout<< std::endl;
	fout.close();

	fout.open(output_name + ".results");
	fout.setf(std::ios::scientific, std::ios::floatfield);
	fout.precision(5);
	fout<< "BARCODE\tNSNP\tAMBIENT";
	for (int i = 0; i < nDonor; ++i) fout<< '\t'<< donor_names[i];
	fout<< std::endl;
	for (int i = 1; i <= ss.ncells; ++i) {
		fout<< ss.cell_barcodes[i];
		nsnp = 0.0;
		for (int j = 0; j < (int)cellxnd[i].size(); ++j)
			for (int k = 0; k < NucDist::size; ++k)
				nsnp += cellxnd[i][j].dist[k];
		fout<< '\t'<< int(nsnp + 1e-2);
		fout<< '\t'<< theta[i][nDonor];
		for (int j = 0; j < nDonor; ++j)
			fout<< '\t'<< theta[i][j];
		fout<< endl;
	}	

	printf("DemuxAlgo::writeOutputs is finished.\n");
}
