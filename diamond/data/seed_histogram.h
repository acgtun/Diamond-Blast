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

#ifndef SEED_HISTOGRAM_H_
#define SEED_HISTOGRAM_H_

#include <boost/thread.hpp>
#include <boost/atomic.hpp>
#include "../basic/seed.h"
#include "sequence_set.h"
#include "../basic/shape_config.h"

using std::vector;
using boost::thread;
using boost::atomic;

typedef size_t shape_histogram[Const::seqp][Const::seedp];

struct seedp_range
{
	seedp_range():
		begin_ (0),
		end_ (0)
	{ }
	seedp_range(unsigned begin, unsigned end):
		begin_ (begin),
		end_ (end)
	{ }
	bool contains(unsigned i) const
	{ return i >= begin_ && i < end_; }
	unsigned begin() const
	{ return begin_; }
	unsigned end() const
	{ return end_; }
	bool lower(unsigned i) const
	{ return i < begin_; }
	bool lower_or_equal(unsigned i) const
	{ return i < end_; }
private:
	unsigned begin_, end_;
} current_range;

size_t partition_size(const shape_histogram &hst, unsigned p)
{
	size_t s (0);
	for(unsigned i=0;i<Const::seqp;++i)
		s += hst[i][p];
	return s;
}

size_t hst_size(const shape_histogram &hst, const seedp_range &range)
{
	size_t s (0);
	for(unsigned i=range.begin();i<range.end();++i)
		s += partition_size(hst, i);
	return s;
}

struct seed_histogram
{

	seed_histogram()
	{ }

	template<typename _val>
	seed_histogram(const Sequence_set<_val> &seqs, const _val&)
	{
		memset(data_, 0, sizeof(data_));
		const vector<shape_config> cfgs (shape_configs<_val>());
		const vector<size_t> seq_partition (seqs.partition());
#pragma omp parallel for schedule(dynamic)
		for(unsigned seqp=0;seqp<Const::seqp;++seqp)
			build_seq_partition(seqs, seqp, seq_partition[seqp], seq_partition[seqp+1], cfgs);
	}

	const shape_histogram& get(unsigned index_mode, unsigned sid) const
	{ return data_[index_mode-1][sid]; }

	size_t max_chunk_size() const
	{
		size_t max (0);
		::partition p (Const::seedp, program_options::lowmem);
		for(unsigned shape=0;shape < shape_config::get().count();++shape)
			for(unsigned chunk=0;chunk < p.parts; ++chunk)
				max = std::max(max, hst_size(data_[program_options::index_mode-1][shape], seedp_range(p.getMin(chunk), p.getMax(chunk))));
		return max;
	}

private:

	template<typename _val>
	void build_seq_partition(const Sequence_set<_val> &seqs,
			const unsigned seqp,
			const size_t begin,
			const size_t end,
			const vector<shape_config> &cfgs)
	{
		assert(seqp < Const::seqp);
		uint64_t key;
		for(size_t i=begin;i<end;++i) {

			assert(i < seqs.get_length());
			const sequence<const _val> seq = seqs[i];
			if(seq.length() < Const::min_shape_len) continue;
			for(unsigned j=0;j<seq.length()+1-Const::min_shape_len; ++j)
				for(vector<shape_config>::const_iterator cfg = cfgs.begin(); cfg != cfgs.end(); ++cfg) {
					assert(cfg->mode() < Const::index_modes);
					assert(cfg->count() <= Const::max_shapes);
					for(unsigned k=0;k<cfg->count(); ++k)
						if(j+cfg->get_shape(k).length_ < seq.length()+1 && cfg->get_shape(k).set_seed(key, &seq[j]))
							++data_[cfg->mode()][k][seqp][seed_partition(key)];
				}

		}
	}

	template<typename _val>
	static vector<shape_config> shape_configs()
	{
		vector<shape_config> v;
		if(program_options::command == program_options::makedb) {
			for(unsigned i=1;i<=Const::index_modes;++i)
				v.push_back(shape_config (i, _val()));
		} else
			v.push_back(shape_config (program_options::index_mode, _val()));
		return v;
	}

	shape_histogram data_[Const::index_modes][Const::max_shapes];

};

#endif /* SEED_HISTOGRAM_H_ */
