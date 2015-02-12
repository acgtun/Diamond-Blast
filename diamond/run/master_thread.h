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

#ifndef MASTER_THREAD_H_
#define MASTER_THREAD_H_

#include <iostream>
#include <boost/timer/timer.hpp>
#include "../data/reference.h"
#include "../data/queries.h"
#include "../basic/statistics.h"
#include "../basic/shape_config.h"
#include "../output/join_blocks.h"
#include "../align/align_queries.h"
#include "../search/align_range.h"
#include "../basic/setup.h"

using std::endl;
using std::cout;
using boost::thread;
using boost::timer::cpu_timer;
using boost::ptr_vector;

template<typename _val, typename _locr, typename _locq, typename _locl>
void process_shape(unsigned sid,
		cpu_timer &timer_mapping,
		unsigned query_chunk,
		unsigned ref_chunk,
		char *query_buffer,
		char *ref_buffer)
{
	using std::vector;
	using boost::atomic;

	::partition p (Const::seedp, program_options::lowmem);
	for(unsigned chunk=0;chunk < p.parts; ++chunk) {

		verbose_stream << "Processing query chunk " << query_chunk << ", reference chunk " << ref_chunk << ", shape " << sid << ", index chunk " << chunk << '.' << endl;
		const seedp_range range (p.getMin(chunk), p.getMax(chunk));
		current_range = range;

		task_timer timer ("Building reference index", true);
		typename sorted_list<_locr>::Type ref_idx (ref_buffer,
				*ref_seqs<_val>::data_,
				shape_config::instance.get_shape(sid),
				ref_hst.get(program_options::index_mode, sid),
				range);
		ref_masking.build<_val,_locr>(sid, range, ref_idx);

		timer.go("Building query index");
		timer_mapping.resume();
		typename sorted_list<_locq>::Type query_idx (query_buffer,
				*query_seqs<_val>::data_,
				shape_config::instance.get_shape(sid),
				query_hst->get(program_options::index_mode, sid),
				range);
		timer.finish();

		timer.go("Searching alignments");
#pragma omp parallel
		{
			Statistics stat;
#pragma omp for schedule(dynamic) nowait
			for(unsigned seedp=0;seedp<Const::seedp;++seedp) {
				try {
					align_partition<_val,_locr,_locq,_locl>(seedp,
							stat,
							sid,
							ref_idx.get_partition_cbegin(seedp),
							query_idx.get_partition_cbegin(seedp));
				} catch (std::exception &e) {
					exception_state.set(e);
				}
			}
#pragma omp critical
			statistics += stat;
		}

	}
	timer_mapping.stop();
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void run_ref_chunk(Database_file &db_file,
		cpu_timer &timer_mapping,
		cpu_timer &total_timer,
		unsigned query_chunk,
		unsigned ref_chunk,
		pair<size_t,size_t> query_len_bounds,
		char *query_buffer)
{
	task_timer timer ("Loading reference sequences", true);
	ref_seqs<_val>::data_ = new Sequence_set<_val> (db_file);
	ref_ids::data_ = new String_set<char,0> (db_file);
	db_file.read(&ref_hst, 1);
	setup_search_params(query_len_bounds, ref_seqs<_val>::data_->letters());

	timer.go("Allocating buffers");
	char *ref_buffer = sorted_list<_locr>::Type::alloc_buffer(ref_hst);

	timer.go("Initializing temporary storage");
	timer_mapping.resume();
	Trace_pt_buffer<_locr,_locl>::instance = new Trace_pt_buffer<_locr,_locl> (query_seqs<_val>::data_->get_length()/query_contexts(),
			program_options::tmpdir,
			program_options::mem_buffered());
	timer.finish();
	timer_mapping.stop();

	for(unsigned i=0;i<shape_config::instance.count();++i)
		process_shape<_val,_locr,_locq,_locl>(i, timer_mapping, query_chunk, ref_chunk, query_buffer, ref_buffer);

	timer.go("Closing temporary storage");
	Trace_pt_buffer<_locr,_locl>::instance->close();
	exception_state.sync();

	timer.go("Deallocating buffers");
	delete[] ref_buffer;

	timer_mapping.resume();
	vector<Output_stream*> out;
	if(ref_header.n_blocks > 1) {
		timer.go ("Opening temporary output files");
		for(unsigned i=0;i<Output_stack<_val>::get().size();++i)
			out.push_back(new Temp_output_file (i, ref_chunk, program_options::tmpdir));
	} else {
		for(typename Output_stack<_val>::iterator i = Output_stack<_val>::get().begin(); i != Output_stack<_val>::get().end(); ++i)
			out.push_back(&i->master_file);
	}

	timer.go("Computing alignments");
	align_queries<_val,_locr,_locl>(*Trace_pt_buffer<_locr,_locl>::instance, out);
	delete Trace_pt_buffer<_locr,_locl>::instance;

	if(ref_header.n_blocks > 1) {
		timer.go("Closing the output files");
		for(vector<Output_stream*>::iterator i = out.begin(); i != out.end(); ++i) {
			(*i)->close();
			delete *i;
		}
	}
	timer_mapping.stop();

	timer.go("Deallocating reference");
	delete ref_seqs<_val>::data_;
	delete ref_ids::data_;
	timer.finish();
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void run_query_chunk(Database_file &db_file,
		cpu_timer &timer_mapping,
		cpu_timer &total_timer,
		unsigned query_chunk,
		pair<size_t,size_t> query_len_bounds)
{
	task_timer timer ("Allocating buffers", true);
	char *query_buffer = sorted_list<_locq>::Type::alloc_buffer(*query_hst);
	timer.finish();

	db_file.rewind();
	for(unsigned ref_chunk=0;ref_chunk<ref_header.n_blocks;++ref_chunk)
		run_ref_chunk<_val,_locr,_locq,_locl>(db_file, timer_mapping, total_timer, query_chunk, ref_chunk, query_len_bounds, query_buffer);

	timer.go("Deallocating buffers");
	timer_mapping.resume();
	delete[] query_buffer;

	timer.go("Deallocating queries");
	delete query_seqs<_val>::data_;
	delete query_ids::data_;

	if(ref_header.n_blocks > 1) {
		timer.go("Joining output blocks");
		join_blocks(ref_header.n_blocks, Output_stack<_val>::get());
	}
	timer_mapping.stop();
}

template<typename _val, typename _locr>
void master_thread(Database_file &db_file, cpu_timer &timer_mapping, cpu_timer &total_timer)
{
	shape_config::instance = shape_config (program_options::index_mode, Amino_acid());

	task_timer timer ("Opening the input file", true);
	timer_mapping.resume();
	const Sequence_file_format<Nucleotide> *format_n (guess_format<Nucleotide>(program_options::query_file));
	const Sequence_file_format<Amino_acid> *format_a (guess_format<Amino_acid>(program_options::query_file));
	Input_stream query_file (program_options::query_file, true);
	current_query_chunk=0;

	timer.go("Opening the output files");
	Output_stack<_val>::get().template add<Blast_tab_format<_val> >(program_options::output_file);
	Output_stack<_val>::get().template add<Sam_format<_val> >(program_options::sam_output);
	timer_mapping.stop();
	timer.finish();

	for(;;++current_query_chunk) {
		task_timer timer ("Loading query sequences", true);
		timer_mapping.resume();
		size_t n_query_seqs;
		if(input_sequence_type() == nucleotide)
			n_query_seqs = load_seqs<Nucleotide,_val>(query_file, *format_n, query_seqs<_val>::data_, query_ids::data_, (size_t)(program_options::chunk_size * 1e9));
		else
			n_query_seqs = load_seqs<Amino_acid,_val>(query_file, *format_a, query_seqs<_val>::data_, query_ids::data_, (size_t)(program_options::chunk_size * 1e9));
		if(n_query_seqs == 0)
			break;
		timer.finish();
		query_seqs<_val>::data_->print_stats();

		if(program_options::seg == "yes") {
			timer.go("Running complexity filter");
			Complexity_filter<_val>::get().run(*query_seqs<_val>::data_);
		}

		timer.go("Building query histograms");
		query_hst = auto_ptr<seed_histogram> (new seed_histogram (*query_seqs<_val>::data_, _val()));
		const pair<size_t,size_t> query_len_bounds = query_seqs<_val>::data_->len_bounds(shape_config::get().get_shape(0).length_);
		timer_mapping.stop();
		timer.finish();
		const bool long_addressing_query = query_seqs<_val>::data_->raw_len() > (size_t)std::numeric_limits<uint32_t>::max();

		if(query_len_bounds.second <= (size_t)std::numeric_limits<uint8_t>::max()) {
			if(long_addressing_query)
				run_query_chunk<_val,_locr,uint64_t,uint8_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds);
			else
				run_query_chunk<_val,_locr,uint32_t,uint8_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds);
		} else if(query_len_bounds.second <= (size_t)std::numeric_limits<uint16_t>::max()) {
			if(long_addressing_query)
				run_query_chunk<_val,_locr,uint64_t,uint16_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds);
			else
				run_query_chunk<_val,_locr,uint32_t,uint16_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds);
		} else {
			if(long_addressing_query)
				run_query_chunk<_val,_locr,uint64_t,uint32_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds);
			else
				run_query_chunk<_val,_locr,uint32_t,uint32_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds);
		}
	}

	timer.go("Closing the output files");
	timer_mapping.resume();
	Output_stack<_val>::get().close_all();
	timer_mapping.stop();

	timer.go("Closing the database file");
	db_file.close();

	timer.finish();
	verbose_stream << "Total time = " << boost::timer::format(total_timer.elapsed(), 1, "%ws\n");
	verbose_stream << "Mapping time = " << boost::timer::format(timer_mapping.elapsed(), 1, "%ws\n");
	statistics.print();
}

template<typename _val>
void master_thread()
{
	cpu_timer timer2, timer_mapping;
	timer_mapping.stop();

	task_timer timer ("Opening the database", 1);
	Database_file db_file;
	timer.finish();
	program_options::set_options<_val>(ref_header.block_size);
	verbose_stream << "Reference = " << program_options::database << endl;
	verbose_stream << "Sequences = " << ref_header.sequences << endl;
	verbose_stream << "Letters = " << ref_header.letters << endl;
	verbose_stream << "Block size = " << (size_t)(ref_header.block_size * 1e9) << endl;

	if(ref_header.long_addressing)
		master_thread<_val,uint64_t>(db_file, timer_mapping, timer2);
	else
		master_thread<_val,uint32_t>(db_file, timer_mapping, timer2);
}

#endif /* MASTER_THREAD_H_ */
