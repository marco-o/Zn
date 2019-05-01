#ifndef znqueue_H
#define znqueue_H

#include <list>
#include <mutex>
#include <condition_variable>

namespace zn
{
	template <class T>
	class shared_list_t
	{
	public:
		bool empty(void) const
		{
			std::unique_lock<std::mutex> lock(mutex_);

			return list_.empty();
		}
		size_t size(void) const
		{
			std::unique_lock<std::mutex> lock(mutex_);

			return list_.size();
		}
		void clear(void)
		{
			std::unique_lock<std::mutex> lock(mutex_);

			list_.clear();
		}
		void push(const T &t)
		{
			std::unique_lock<std::mutex> lock(mutex_);

			list_.push_back(t);
			cond_.notify_one();
		}
		T pop(void)
		{
			std::unique_lock<std::mutex> lock(mutex_);
			cond_.wait(lock, [this](void) { return !list_.empty(); });
			T result = *list_.begin();
			list_.pop_front();
			return result;
        }
		bool pop(T &t)
		{
			std::unique_lock<std::mutex> lock(mutex_);
			if (list_.empty())
				return false;
			t = *list_.begin();
			list_.pop_front();
			return true;
		}
	private:
		std::list<T>			list_;
		mutable std::mutex		mutex_;
		std::condition_variable cond_;
	};
}
#endif
