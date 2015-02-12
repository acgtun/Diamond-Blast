/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef TRANSLATE_H_
#define TRANSLATE_H_

struct Translator
{

public:

	static const Nucleotide reverseNucleotide[5];
	static const Amino_acid lookup[5][5][5];
	static const Amino_acid lookupReverse[5][5][5];
	static const Amino_acid STOP;

	static Nucleotide getReverseComplement(Nucleotide nucleotide)
	{ return reverseNucleotide[(int) nucleotide]; }

	static Amino_acid getAminoAcid(vector<Nucleotide> const &dnaSequence, size_t pos)
	{ return lookup[(int) dnaSequence[pos]][(int)dnaSequence[pos+1]][(int)dnaSequence[pos+2]]; }

	static Amino_acid getAminoAcidReverse(vector<Nucleotide> const &dnaSequence, size_t pos)
	{ return lookupReverse[(int) dnaSequence[pos + 2]][(int)dnaSequence[pos + 1]][(int)dnaSequence[pos]]; }

	static int translate(vector<Nucleotide> const &dnaSequence, vector<Amino_acid> *proteinSequences)
	{
		size_t length_ = dnaSequence.size();
		size_t r = length_ - 2;
		unsigned pos = 0;
		unsigned i = 0;
		while(r > 2) {
			proteinSequences[0][i] = getAminoAcid(dnaSequence, pos++);
			proteinSequences[3][i] = getAminoAcidReverse(dnaSequence, --r);
			proteinSequences[1][i] = getAminoAcid(dnaSequence, pos++);
			proteinSequences[4][i] = getAminoAcidReverse(dnaSequence, --r);
			proteinSequences[2][i] = getAminoAcid(dnaSequence, pos++);
			proteinSequences[5][i] = getAminoAcidReverse(dnaSequence, --r);
			++i;
		}
		if(r) {
			proteinSequences[0][i] = getAminoAcid(dnaSequence, pos++);
			proteinSequences[3][i] = getAminoAcidReverse(dnaSequence, --r);
		}
		if(r) {
			proteinSequences[1][i] = getAminoAcid(dnaSequence, pos);
			proteinSequences[4][i] = getAminoAcidReverse(dnaSequence, r);
		}
		return 6;
	}

	static Amino_acid const* nextChar(Amino_acid const*p, Amino_acid const*end)
	{
		while(*(p) != STOP && p < end)
			++p;
		return p;
	}

	static void mask_runs(vector<Amino_acid> &query, unsigned run_len)
	{
		Amino_acid *last = &query[0]-1, *i = &query[0], *end = &query.back();
		while (i <= end) {
			if(*i == STOP) {
				if(last != 0 && i - last - 1 < run_len) {
					for(Amino_acid *j = last+1; j < i; ++j)
						*j = Value_traits<Amino_acid>::MASK_CHAR;
				}
				last = i;
			}
			++i;
		}
		if(last != 0 && i - last - 1 < run_len) {
			for(Amino_acid *j = last+1; j < i; ++j)
				*j = Value_traits<Amino_acid>::MASK_CHAR;
		}
	}

	static unsigned computeGoodFrames(vector<Amino_acid>  const *queries, unsigned runLen)
	{
		unsigned set = 0;

		for (unsigned i = 0; i < 6; ++i) {
			if (queries[i].size() > 0) {
				unsigned run = 0;
				Amino_acid const*p =  &(queries[i][0]);
				Amino_acid const*q;
				Amino_acid const*end = p + queries[i].size();
				while((q = nextChar(p, end)) < end) {
					run = q-p;
					if (run >= runLen)
						set |= 1 << i;
					p=q+1;
				}
				run = q-p;
				if (run >= runLen)
					set |= 1 << i;
			}
		}
		return set;
	}

	static void mask_runs(vector<Amino_acid> *queries, unsigned run_len)
	{
		for (unsigned i = 0; i < 6; ++i)
			mask_runs(queries[i], run_len);
	}

};

const Nucleotide Translator::reverseNucleotide[5] = { 3, 2, 1, 0, 4 };

const Amino_acid Translator::lookup[5][5][5] = {
{ { 11,2,11,2,Value_traits<Amino_acid>::MASK_CHAR },
{ 16,16,16,16,16 },
{ 1,15,1,15,Value_traits<Amino_acid>::MASK_CHAR },
{ 9,9,12,9,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
 },
{ { 5,8,5,8,Value_traits<Amino_acid>::MASK_CHAR },
{ 14,14,14,14,14 },
{ 1,1,1,1,1 },
{ 10,10,10,10,10 },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
 },
{ { 6,3,6,3,Value_traits<Amino_acid>::MASK_CHAR },
{ 0,0,0,0,0 },
{ 7,7,7,7,7 },
{ 19,19,19,19,19 },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
 },
{ { 23,18,23,18,Value_traits<Amino_acid>::MASK_CHAR },
{ 15,15,15,15,15 },
{ 23,4,17,4,Value_traits<Amino_acid>::MASK_CHAR },
{ 10,13,10,13,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
 },
{ { Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
} };

const Amino_acid Translator::lookupReverse[5][5][5] = {
{ { 13,10,13,10,Value_traits<Amino_acid>::MASK_CHAR },
{ 4,17,4,23,Value_traits<Amino_acid>::MASK_CHAR },
{ 15,15,15,15,15 },
{ 18,23,18,23,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
 },
{ { 19,19,19,19,19 },
{ 7,7,7,7,7 },
{ 0,0,0,0,0 },
{ 3,6,3,6,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
 },
{ { 10,10,10,10,10 },
{ 1,1,1,1,1 },
{ 14,14,14,14,14 },
{ 8,5,8,5,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
 },
{ { 9,12,9,9,Value_traits<Amino_acid>::MASK_CHAR },
{ 15,1,15,1,Value_traits<Amino_acid>::MASK_CHAR },
{ 16,16,16,16,16 },
{ 2,11,2,11,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
 },
{ { Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
{ Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR,Value_traits<Amino_acid>::MASK_CHAR },
}};

const Amino_acid Translator::STOP (Value_traits<Amino_acid>::from_char('*'));

#endif /* TRANSLATE_H_ */
