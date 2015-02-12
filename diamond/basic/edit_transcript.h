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

#ifndef EDIT_TRANSCRIPT_H_
#define EDIT_TRANSCRIPT_H_

#include <vector>

using std::vector;
using std::endl;

struct Edit_transcript : public vector<char>
{

	typedef enum { match, insertion, deletion } Operation;

	static const char *op_char;
	static const int op_gaps[];

	void push(Operation op, unsigned n)
	{
		for(unsigned i=0;i<n;++i)
			this->push_back(op);
	}

	size_t print_cigar(char *ptr) const
	{
		char *ptr2 = ptr;
		for(const_iterator i = this->end()-1; i >= this->begin();)
			ptr2 += print(ptr2, (Operation)*i, i);
		return ptr2 - ptr;
	}

	unsigned gap_positions() const
	{
		unsigned n = 0;
		for(const_iterator i = this->end()-1; i >= this->begin(); --i)
			n += op_gaps[(long)*i];
		return n;
	}

	template<typename _val>
	size_t print_MD(char *ptr, const sequence<const _val> &query, const sequence<const _val> &subject, unsigned qpos, unsigned spos)
	{
		char *ptr2 = ptr;
		unsigned matches = 0;
		for(const_iterator i = this->end()-1; i >= this->begin();)
			switch(*i) {
			case match:
				if(query[qpos] == mask_critical(subject[spos]))
					++matches;
				else {
					print_matches(ptr2, matches);
					*(ptr2++) = Value_traits<_val>::ALPHABET[mask_critical(subject[spos])];
				}
				++qpos;
				++spos;
				--i;
				break;
			case insertion:
				//*(ptr2++) = Value_traits<_val>::ALPHABET[(query[qpos++])];
				++qpos;
				--i;
				break;
			case deletion:
				print_matches(ptr2, matches);
				*(ptr2++) = '^';
				while(i >= this->begin() && *i == deletion) {
					--i;
					*(ptr2++) = Value_traits<_val>::ALPHABET[mask_critical(subject[spos++])];
				}
				if(query[qpos] != mask_critical(subject[spos]))
					*(ptr2++) = '0';
			}
		print_matches(ptr2, matches);
		return ptr2 - ptr;
	}

	template<typename _val>
	void print(std::ostream &os, const _val *query, const _val *subject)
	{
		print(os, query, deletion);
		os << endl;
		print(os, subject, insertion);
		os << endl;
	}

	template<typename _val>
	void print(std::ostream &os, const _val *s, Operation gap_op)
	{
		for(const_iterator i = this->end()-1; i >= this->begin(); --i)
			if(*i == gap_op)
				os << '-';
			else
				os << Value_traits<_val>::ALPHABET[*(s++)];
	}

private:

	void print_matches(char *&ptr, unsigned &n)
	{
		if(n > 0) {
			ptr += sprintf(ptr, "%u", n);
			n = 0;
		}
	}

	size_t print(char *ptr, Operation op, const_iterator &i) const
	{
		unsigned n = 0;
		while(i >= this->begin() && *i == op) {
			++n;
			--i;
		}
		return sprintf(ptr, "%u%c", n, op_char[op]);
	}

};

const char* Edit_transcript::op_char = "MID";
const int Edit_transcript::op_gaps[] = { 0, 1, 1 };

#endif /* EDIT_TRANSCRIPT_H_ */
