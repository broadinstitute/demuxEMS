#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "VCFLoader.hpp"
#include "SufficientStatistics.hpp"
#include "DemuxAlgo.hpp"

using namespace std;

SufficientStatistics ss;
DemuxAlgo *algo;

int main(int argc, char* argv[]) {
	if (argc != 3) {
		printf("Usage: debug input.params input.data\n");
		exit(-1);
	}

	SNPType::setNumDonor(8);
	algo = new DemuxAlgo(argv[1], ss);

	ifstream fin(argv[2]);
	string barcode, line;
	int ngeno, geno_id;
	char c;

	ss.ncells = 1;
	fin>>barcode>> ngeno;
	getline(fin, line);

	ss.cell_barcodes.push_back(std::move(""));
	ss.cell_barcodes.push_back(barcode);
	ss.cellxgeno.resize(2);
	ss.cellxnd.resize(2);
	ss.cellxgeno[1].resize(ngeno);
	ss.cellxnd[1].resize(ngeno);

	for (int i = 0; i < ngeno; ++i) {
		fin>> geno_id>> c;
		for (int j = 0; j < 3; ++j) fin>> ss.cellxnd[1][i].dist[j];
		fin>> c>> c;
		SNPType *snp = new SNPType(0, 'A', 'C');
		int geno;
		for (int j = 0; j < 8; ++j) {
			fin>> geno;
			string geno_str = (geno == 0 ? "0/0" : (geno == 1 ? "0/1" : "1/1"));
			snp->setDonorGenotype(j, geno_str);
		}
		getline(fin, line);
		ss.cellxgeno[1][i] = geno_id;
		ss.diff_genos.push_back(snp);
	}
	fin.close();

	// cout<< ss.cell_barcodes[1]<< endl;
	// for (auto&& val : ss.diff_genos) {
	// 	for (int j = 0; j < 8; ++j) cout<< val->getDonorGenotype(j)<< "\t";
	// 	cout<< endl;
	// }
	// cout<< endl;
	// for (auto&& val : ss.cellxgeno[1]) cout<< val<< "\t"; cout<<endl;
	// for (auto&& val : ss.cellxnd[1]) cout<< val.dist[0]<< "\t"<< val.dist[1]<< "\t"<< val.dist[2]<< endl;

	algo->demultiplex(0, 1, 2, 1.0, 1.0, 1e-6);

	delete algo;

	return 0;
}
