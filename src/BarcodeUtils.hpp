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

#ifndef BARCODEUTILS_H_
#define BARCODEUTILS_H_

#include <cassert>
#include <cstdint>
#include <string>
#include <vector>

const int STEP = 3;
const int BASE = 7;
const int UPPER = 21;
const int NNUC = 5; // ACGTN

const char id2base[NNUC] = {'A', 'C', 'G', 'T', 'N'};

// Assume char's range is -128..127
const int CHAR_RANGE = 128;
static std::vector<int> init_base2id() {
	std::vector<int> vec(CHAR_RANGE, -1);
	vec['a'] = vec['A'] = 0;
	vec['c'] = vec['C'] = 1;
	vec['g'] = vec['G'] = 2;
	vec['t'] = vec['T'] = 3;
	vec['n'] = vec['N'] = 4;

	return vec;
}

static const std::vector<int> base2id = init_base2id();



static std::vector<std::vector<uint64_t> > init_aux_arr() {
	std::vector<std::vector<uint64_t> > aux_arr;
	for (int i = 0; i < UPPER; ++i) {
		std::vector<uint64_t> arr(NNUC + 1, 0);
		for (uint64_t j = 0; j < NNUC; ++j) arr[j] = j << (STEP * i);
		arr[NNUC] = uint64_t(BASE) << (STEP * i);
		aux_arr.push_back(arr);
	}
	return aux_arr;
}

static const std::vector<std::vector<uint64_t> > aux_arr = init_aux_arr();

uint64_t barcode_to_binary(const std::string& barcode) {
	uint64_t binary_id = 0;
	char c;

	assert(barcode.length() <= UPPER);
	for (auto&& it = barcode.begin(); it != barcode.end(); ++it) {
		c = *it;
		assert(base2id[c] >= 0);
		binary_id <<= STEP;
		binary_id += base2id[c];
	}
	return binary_id;
}

std::string binary_to_barcode(uint64_t binary_id, int len) {
	std::string barcode(len, 0);
	for (int i = len - 1; i >= 0; --i) {
		barcode[i] = id2base[binary_id & BASE];
		binary_id >>= STEP;
	}
	return barcode;
}

#endif 
