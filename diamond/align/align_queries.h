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

#ifndef ALIGN_QUERIES_H_
#define ALIGN_QUERIES_H_

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3)
#include <parallel/algorithm>
#else
#include "../util/merge_sort.h"
#endif
#include "../search/trace_pt_buffer.h"
#include "../util/map.h"
#include "align_read.h"
#include "../util/task_queue.h"

using std::vector;
using boost::thread;

struct Output_piece
{
	Output_piece():
		text (),
		size ()
	{ }
	Output_piece(char *text, size_t size):
		text (text),
		size (size)
	{ }
	char *text;
	size_t size;
};

struct Output_pieces
{
	Output_piece data[Const::output_files];
};

struct Output_writer
{
	Output_writer(const vector<Output_stream*> &f):
		f_ (f)
	{ }
	void operator()(const Output_pieces &p)
	{
		for(unsigned i=0;i<f_.size();++i) {
			if(p.data[i].text != 0)
				f_[i]->write(p.data[i].text, p.data[i].size);
			free(p.data[i].text);
		}
	}
private:
	const vector<Output_stream*> &f_;
};

template<typename _val, typename _locr, typename _locl, unsigned _d>
void align_queries(typename Trace_pt_buffer<_locr,_locl>::Vector &trace_pts,
		size_t begin,
		size_t end,
		ptr_vector<Output_buffer<_val> > &buffers,
		Statistics &st)
{
	typedef Map<typename vector<hit<_locr,_locl> >::iterator,typename hit<_locr,_locl>::template Query_id<_d> > Map_t;
	Map_t hits (trace_pts.begin() + begin, trace_pts.begin() + end);
	typename Map_t::Iterator i = hits.begin();
	while(i.valid() && !exception_state()) {
		align_read<_val,_locr,_locl>(buffers, st, i.begin(), i.end());
		++i;
	}
}

template<typename _val, typename _locr, typename _locl>
void align_queries(typename Trace_pt_buffer<_locr,_locl>::Vector &trace_pts, const vector<Output_stream*> &output_files)
{
	const size_t max_size=65536,max_segments=4096,min_segments=program_options::threads()*4;
	if(trace_pts.size() == 0)
		return;
	vector<size_t> p;
	if(query_contexts() == 6)
		p = map_partition(trace_pts.begin(), trace_pts.end(), hit<_locr,_locl>::template query_id<6>, max_size, max_segments, min_segments);
	else
		p = map_partition(trace_pts.begin(), trace_pts.end(), hit<_locr,_locl>::template query_id<1>, max_size, max_segments, min_segments);
	const size_t n = p.size() - 1;

	Output_writer writer (output_files);
	Task_queue<Output_pieces,Output_writer> queue (n, program_options::threads()*2, writer);

#pragma omp parallel
	{
		Statistics st;
		size_t i;
		while(queue.get(i) && !exception_state()) {
			try {
				size_t begin = p[i], end = p[i+1];
				ptr_vector<Output_buffer<_val> > buffers = Output_stack<_val>::get().get_buffers();
				if(query_contexts() == 6)
					align_queries<_val,_locr,_locl,6>(trace_pts, begin, end, buffers, st);
				else
					align_queries<_val,_locr,_locl,1>(trace_pts, begin, end, buffers, st);
				Output_pieces p;
				for(unsigned j=0;j<buffers.size();++j)
					p.data[j] = Output_piece(buffers[j].get_begin(), buffers[j].size());
				queue.push(i, p);
			} catch(std::exception &e) {
				exception_state.set(e);
			}
		}
#pragma omp critical
		statistics += st;
	}
	exception_state.sync();
}

template<typename _val, typename _locr, typename _locl>
void align_queries(const Trace_pt_buffer<_locr,_locl> &trace_pts, const vector<Output_stream*> &output_files)
{
	typename Trace_pt_buffer<_locr,_locl>::Vector v;
	size_t q_len = max_id_len(query_ids::get()), ref_len = max_id_len(ref_ids::get());
	log_stream << "ID len " << ref_len << ' ' << q_len << endl;
	Output_stack<_val>::get().set_line_sizes(q_len, ref_len, query_seqs<_val>::data_->len_bounds(0).second);
	//ref_bam_map.init(ref_seqs<_val>::get().get_length());
	for(unsigned bin=0;bin<trace_pts.bins();++bin) {
		log_stream << "Processing query bin " << bin+1 << '/' << trace_pts.bins() << '\n';
		task_timer timer ("Loading trace points", false);
		trace_pts.load(v, bin);
		timer.go("Sorting trace points");
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3)
		__gnu_parallel::sort(v.begin(), v.end());
#else
		merge_sort(v.begin(), v.end(), program_options::threads());
#endif
		timer.go("Computing alignments");
		align_queries<_val,_locr,_locl>(v, output_files);
	}
	//ref_bam_map.finish();
}

#endif /* ALIGN_QUERIES_H_ */
