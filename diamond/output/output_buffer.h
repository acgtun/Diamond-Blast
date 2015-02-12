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

#ifndef OUTPUT_BUFFER_H_
#define OUTPUT_BUFFER_H_

#include "../util/text_buffer.h"
#include "output_format.h"

struct Segment_meta_info
{
	Segment_meta_info()
	{ }
	Segment_meta_info(unsigned query_id, unsigned subject_id, int score):
		query_id (query_id),
		subject_id (subject_id),
		len (0),
		score (score)
	{ }
	unsigned query_id, subject_id, len;
	int score;
};

template<typename _val>
struct Output_buffer : public Text_buffer
{

	virtual void print_match(const Output_format<_val> &format,
			const match<_val> &match,
			size_t query_source_len,
			const sequence<const char> &query_name,
			const sequence<const _val> &query,
			unsigned query_id,
			size_t res)
	{ this->reserve(res); format.print_match(*this, match, query_source_len, query_name, query); }

	virtual ~Output_buffer()
	{ }

};

template<typename _val>
struct Temp_output_buffer : public Output_buffer<_val>
{

	virtual void print_match(const Output_format<_val> &format,
			const match<_val> &match,
			size_t query_source_len,
			const sequence<const char> &query_name,
			const sequence<const _val> &query,
			unsigned query_id,
			size_t res)
	{
		Segment_meta_info i (query_id, match.subject_id_, match.score_);
		this->write(i);
		this->reserve(res);
		size_t n = format.print_match(*this, match, query_source_len, query_name, query);
		reinterpret_cast<Segment_meta_info*>(this->operator char *() - n - sizeof(Segment_meta_info))->len = n;
	}

	virtual ~Temp_output_buffer()
	{ }

};

#endif /* OUTPUT_BUFFER_H_ */
