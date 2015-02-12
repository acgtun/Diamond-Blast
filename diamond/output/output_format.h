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

#ifndef OUTPUT_FORMAT_H_
#define OUTPUT_FORMAT_H_

#include "../basic/match.h"
#include "../align/match_func.h"

template<typename _val>
struct Output_format
{
	virtual size_t print_match(Text_buffer &buf,
			const match<_val> &match,
			size_t query_source_len,
			const sequence<const char> &query_name,
			const sequence<const _val> &query) const = 0;
	virtual void print_header(Output_stream &file) const
	{ }
	virtual unsigned max_line_length(size_t max_query_id, size_t max_ref_id, size_t max_query_len) const = 0;
	virtual ~Output_format()
	{ }
};

template<typename _val>
struct Blast_tab_format : public Output_format<_val>
{

	virtual size_t print_match(Text_buffer &buf,
			const match<_val> &match,
			size_t query_source_len,
			const sequence<const char> &query_name,
			const sequence<const _val> &query) const
	{
		assert(match.frame_ < 6);
		//assert(match.traceback_ != 0);

		int raw_score = match.traceback_ != 0 ? match.traceback_->score_ : match.score_;
		local_match<_val> l;

		if(program_options::alignment_traceback)
			l = *match.traceback_;

		size_t n, t = 0;

		buf += (n = print_str(buf, query_name.c_str(), Const::id_delimiters));
		t += n;
		buf += (n = sprintf(buf, "\t"));
		t += n;
		buf += (n = print_str(buf, ref_ids::get()[match.subject_id_].c_str(), Const::id_delimiters));
		t += n;
		buf += (n = sprintf(buf, "\t%.2f\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%.0e\t%.1f",
				static_cast<float>(l.identities_)*100/l.len_,
				l.len_,
				l.mismatches_,
				l.gap_openings_,
				l.query_begin_+1,
				l.query_begin_+1 + (l.query_len_ > 0 ? -1 : 1) + l.query_len_,
				l.subject_begin_+1,
				l.subject_begin_+1 + l.subject_len_ - 1,
				match.evalue_,
				score_matrix::get().bitscore(raw_score)));
		t += n;
		if(program_options::salltitles) {
			buf << '\t';
			t += print_salltitles(buf, ref_ids::get()[match.subject_id_].c_str()) + 1;
		}
		buf << '\n';
		return t+1;
	}

	virtual unsigned max_line_length(size_t max_query_id, size_t max_ref_id, size_t max_query_len) const
	{ return max_query_id + max_ref_id + 128; }

	virtual ~Blast_tab_format()
	{ }

	static size_t print_salltitles(Text_buffer &buf, const char *id)
	{
		size_t n = 0;
		const vector<string> t (tokenize(id, "\1"));
		vector<string>::const_iterator i=t.begin();
		for(;i<t.end()-1;++i) {
			buf << *i << "<>";
			n += i->length() + 2;
		}
		buf << *i;
		n += i->length();
		return n;
	}

	static const Blast_tab_format instance;

};

template<typename _val>
const Blast_tab_format<_val> Blast_tab_format<_val>::instance;

template<typename _val>
struct Sam_format : public Output_format<_val>
{

	virtual size_t print_match(Text_buffer &buf,
			const match<_val> &match,
			size_t query_source_len,
			const sequence<const char> &query_name,
			const sequence<const _val> &query) const
	{
		int raw_score = match.traceback_ != 0 ? match.traceback_->score_ : match.score_;
		local_match<_val> l;
		const unsigned query_aligned_len = abs(match.traceback_->query_len_ / query_len_factor());

		if(program_options::alignment_traceback)
			l = *match.traceback_;

		size_t n, t = 0;

		buf += (n = print_str(buf, query_name.c_str(), Const::id_delimiters));
		t += n;
		buf += (n = sprintf(buf, "\t%i\t", 0));
		t += n;
		buf += (n = print_str(buf, ref_ids::get()[match.subject_id_].c_str(), Const::id_delimiters));
		t += n;
		buf += (n = sprintf(buf, "\t%u\t%u\t",
				l.subject_begin_+1,
				255));
		t += n;

		if(program_options::alignment_traceback) {
			if(l.transcript_)
				buf += (n = l.transcript_->print_cigar(buf));
			else
				buf += (n = sprintf(buf, "%uM", query_aligned_len));
		} else
			buf += (n = sprintf(buf, "*"));
		t += n;

		const unsigned qbegin = query_translated() ? query_translated_begin(*match.traceback_, match.frame_, query_source_len) : l.query_begin_;

		buf += (n = sprintf(buf, "\t*\t0\t0\t"));
		t += n;
		buf += (n = query.print(buf, qbegin, query_aligned_len));
		t += n;
		buf += (n = sprintf(buf, "\t*\tAS:i:%u\tNM:i:%u\tZL:i:%lu\tZR:i:%i\tZE:f:%.1e\tZI:i:%u\tZF:i:%i\tZS:i:%u\tMD:Z:",
				(unsigned)score_matrix::get().bitscore(raw_score),
				l.len_ - l.identities_,
				ref_seqs<_val>::data_->length(match.subject_id_),
				raw_score,
				match.evalue_,
				l.identities_*100/l.len_,
				blast_frame(match.frame_),
				l.query_begin_+1));
		t += n;
		if(l.transcript_)
			buf += (n = l.transcript_->print_MD(buf, query, ref_seqs<_val>::get()[match.subject_id_], qbegin, l.subject_begin_));
		else
			buf += (n = print_MD(buf, query, ref_seqs<_val>::get()[match.subject_id_], qbegin, l.subject_begin_, l.len_));
		t += n;
		buf += sprintf(buf, "\n");
		return t+1;
	}

	virtual void print_header(Output_stream &file) const
	{
		static const char* line = "@HD\tVN:1.5\tSO:query\n\
@PG\tPN:DIAMOND\n\
@mm\tBlastX\n\
@CO\tBlastX-like alignments\n\
@CO\tReporting AS: bitScore, ZR: rawScore, ZE: expected, ZI: percent identity, ZL: reference length, ZF: frame, ZS: query start DNA coordinate\n";
		file.write(line, strlen(line));
	}

	virtual unsigned max_line_length(size_t max_query_id, size_t max_ref_id, size_t max_query_len) const
	{ return max_query_id + max_ref_id + 2*max_query_len + 128; }

	virtual ~Sam_format()
	{ }

	static const Sam_format instance;

};

template<typename _val>
struct Bam_format : public Output_format<_val>
{

	virtual size_t print_match(Text_buffer &buf,
			const match<_val> &match,
			size_t query_source_len,
			const sequence<const char> &query_name,
			const sequence<const _val> &query) const
	{
		int raw_score = match.traceback_ != 0 ? match.traceback_->score_ : match.score_;
		local_match<_val> l;
		const unsigned query_aligned_len = abs(match.traceback_->query_len_ / query_len_factor());

		if(program_options::alignment_traceback)
			l = *match.traceback_;

		//buf << 0;
	}

	virtual ~Bam_format()
	{ }

};

template<typename _val> const Sam_format<_val> Sam_format<_val>::instance;

template<typename _val>
size_t print_MD(char *ptr, const sequence<const _val> &query, const sequence<const _val> &subject, unsigned qpos, unsigned spos, unsigned len)
{
	unsigned n = 0;
	char *ptr2 = ptr;
	for(unsigned i=0;i<len;++i)
		if(query[qpos+i] == mask_critical(subject[spos+i]))
			++n;
		else {
			if(n) {
				ptr2 += sprintf(ptr2, "%u", n);
				n = 0;
			}
			*(ptr2++) = Value_traits<_val>::ALPHABET[mask_critical(subject[spos+i])];
		}
	if(n)
		ptr2 += sprintf(ptr2, "%u", n);
	return ptr2 - ptr;
}

#endif /* OUTPUT_FORMAT_H_ */
