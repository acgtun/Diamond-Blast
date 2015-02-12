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

#ifndef OUTPUT_FILE_H_
#define OUTPUT_FILE_H_

#include <string>
#include "output_buffer.h"

using std::string;

struct Temp_output_file : public Output_stream
{

	Temp_output_file(unsigned filter, unsigned ref_block, const string &tmpdir):
		Output_stream (file_name(tmpdir, filter, ref_block), program_options::compress_temp == 1)
	{ }

	static string file_name(const string &tmpdir, unsigned filter, unsigned ref_block)
	{ return tmpdir + "/diamond_out_" + boost::to_string(program_options::magic_number) + "_" + boost::to_string(filter) + "_" + boost::to_string(ref_block) + ".tmp"; }

};

struct Block_output : public Buffered_file
{

	struct Iterator
	{
		unsigned block_;
		Segment_meta_info info_;
		bool operator<(const Iterator &rhs) const
		{ return info_.query_id > rhs.info_.query_id || (info_.query_id == rhs.info_.query_id && info_.score < rhs.info_.score); }
	};

	bool next(Iterator &it)
	{
		if(this->eof())
			return false;
		this->read(it.info_);
		it.block_ = block_;
		return true;
	}

	void skip(const Iterator &it)
	{
		char c;
		for(unsigned i=0;i<it.info_.len;++i)
			this->read(c);
	}

	void copy(Buffered_ostream &dest, const Iterator &it)
	{
		char c;
		for(unsigned i=0;i<it.info_.len;++i) {
			this->read(c);
			dest.write(c);
		}
	}

	Block_output(unsigned filter, unsigned ref_block, const string &tmpdir):
		Buffered_file (Temp_output_file::file_name(tmpdir, filter, ref_block), program_options::compress_temp == 1),
		block_ (ref_block)
	{ }

private:

	const unsigned block_;

};

#endif /* OUTPUT_FILE_H_ */
