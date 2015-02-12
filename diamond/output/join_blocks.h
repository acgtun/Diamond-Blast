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

#ifndef JOIN_BLOCKS_H_
#define JOIN_BLOCKS_H_

#include <algorithm>
#include <vector>
#include "output_file.h"
#include "output_stack.h"

using std::endl;
using std::cout;
using std::vector;

void join_blocks(unsigned ref_blocks, Output_stream *master_out, unsigned filter)
{
	vector<Block_output*> files;
	vector<Block_output::Iterator> records;
	Block_output::Iterator r;
	for(unsigned i=0;i<ref_blocks;++i) {
		files.push_back(new Block_output (filter, i, program_options::tmpdir));
		if(files.back()->next(r))
			records.push_back(r);
	}
	std::make_heap(records.begin(), records.end());
	unsigned query, block, subject, n = 0;
	query = block = subject = std::numeric_limits<unsigned>::max();
	int top_score=0;
	Buffered_ostream master (*master_out);
	while(!records.empty()) {
		const Block_output::Iterator &next = records.front();
		const unsigned b = next.block_;

		if(next.info_.query_id != query) {
			query = next.info_.query_id;
			n = 0;
			top_score = next.info_.score;
			statistics.inc(Statistics::ALIGNED);
		}
		const bool same_subject = b == block && next.info_.subject_id == subject;
		if(program_options::output_range(n, next.info_.score, top_score) || same_subject) {
			files[b]->copy(master, next);
			statistics.inc(Statistics::MATCHES);
			if(!same_subject) {
				block = b;
				subject = next.info_.subject_id;
				++n;
			}
		} else
			files[b]->skip(next);

		std::pop_heap(records.begin(), records.end());
		records.pop_back();
		if(files[b]->next(r)) {
			records.push_back(r);
			std::push_heap(records.begin(), records.end());
		}
	}
	for(unsigned i=0;i<ref_blocks;++i) {
		files[i]->close();
		files[i]->remove();
		delete files[i];
	}
}

template<typename _val>
void join_blocks(unsigned ref_blocks, Output_stack<_val> &stack)
{
	for(unsigned i=0;i<stack.size();++i)
		join_blocks(ref_blocks, &stack[i].master_file, i);
}

#endif /* JOIN_BLOCKS_H_ */
