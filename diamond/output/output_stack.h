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

#ifndef OUTPUT_STACK_H_
#define OUTPUT_STACK_H_

#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>
#include "output_buffer.h"
#include "output_filter.h"

using std::vector;
using std::auto_ptr;
using std::string;
using boost::ptr_vector;

template<typename _val>
struct Output_stack : public ptr_vector<Output_filter<_val> >
{

	ptr_vector<Output_buffer<_val> > get_buffers() const
	{
		ptr_vector<Output_buffer<_val> > res;
		for(unsigned i=0;i<this->size();++i)
			res.push_back(ref_header.n_blocks > 1 ? new Temp_output_buffer<_val> () : new Output_buffer<_val> ());
		return res;
	}

	template<typename _fmt>
	void add(const string &file_name)
	{
		if(file_name.length() != 0)
			this->push_back(new Output_filter<_val> (_fmt::instance, file_name));
	}

	void close_all()
	{
		for(typename Output_stack::iterator i = this->begin(); i != this->end(); ++i)
			i->master_file.close();
	}

	void set_line_sizes(size_t max_query_id, size_t max_ref_id, size_t max_query_len)
	{
		for(typename Output_stack::iterator i = this->begin(); i != this->end(); ++i)
			i->max_line_size = i->format.max_line_length(max_query_id, max_ref_id, max_query_len);
	}

	static Output_stack& get()
	{ return instance; }

private:

	static Output_stack instance;

};

template<typename _val> Output_stack<_val> Output_stack<_val>::instance;

#endif /* OUTPUT_STACK_H_ */
