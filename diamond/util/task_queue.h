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

template<typename _t, typename _callback>
struct Task_queue
{

	Task_queue(unsigned n, unsigned limit, _callback &callback):
		head_ (0),
		tail_ (0),
		limit_ (limit),
		end_ (n),
		callback_ (callback)
	{ }

	volatile bool waiting() const
	{ return tail_ - head_ >= limit_; }

	volatile bool get(size_t &n)
	{
		{
			boost::unique_lock<boost::mutex> lock (mtx_);
			if(tail_ >= end_)
				return false;
			while(waiting() && tail_ < end_)
				cond_.wait(lock);
			if(tail_ >= end_)
				return false;
			n = tail_++;
		}
		if(tail_ == end_)
			cond_.notify_all();
		return true;
	}

	volatile void push(unsigned n, const _t& v)
	{
		mtx_.lock();
		if(n == head_) {
			mtx_.unlock();
			flush(v);
		} else {
			queue_.push(Item (n, v));
			mtx_.unlock();
		}
	}

	volatile void flush(const _t &v)
	{
		callback_(v);
		bool notify, next=false;
		_t nextv;
		{
			boost::lock_guard<boost::mutex> lock(mtx_);
			notify = waiting();
			++head_;
			if(!queue_.empty() && queue_.top().n == head_) {
				nextv = queue_.top().value;
				queue_.pop();
				next = true;
			}
		}
		if(notify)
			cond_.notify_one();
		if(next)
			flush(nextv);
	}

private:

	struct Item
	{
		unsigned n;
		_t value;
		Item(unsigned n, const _t &value):
			n (n),
			value (value)
		{ }
		bool operator<(const Item& rhs) const
		{ return n > rhs.n; }
	};

	std::priority_queue<Item> queue_;
	boost::mutex mtx_;
	boost::condition_variable cond_;
	volatile unsigned head_, tail_, limit_, end_;
	_callback &callback_;

};
