#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "VCFLoader.hpp"

using namespace std;

VCFLoader vcf_loader;

int main(int argc, char* argv[]) {
	if (argc != 2) {
		printf("Usage: test input.vcf");
	}

	string line;
	ifstream fin(argv[1]);
	cout<< getline(fin, line).good()<< endl;

	return 0;
}
