#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <limits>
#include <vector>

template<class T, class U>
class Harmonic_Rec
	{
public:
	typedef T data_t;
	
	T operator()()
		{
		return T(0);
		}

	T operator()(T last, U curr)
		{
		return last + T(1)/T(curr);
		}
	};

template<class T, class U>
class Harmonic_2_Rec
	{
public:
	typedef T data_t;
	
	T operator()()
		{
		return T(0);
		}

	T operator()(T last, U curr)
		{
		return last + T(1)/T(curr*curr);
		}
	};

template<class OP>
class RecursiveSeries
	{
public:
	typedef typename OP::data_t data_t;
	
protected:
	std::vector<data_t> _data;
	OP _op;

public:
	RecursiveSeries()
		: _data(1)
		{
		_data[0] = _op();
		}
	
	template<class I>
	data_t operator()(I n)
		{
		for (size_t i=_data.size(); i<=n; i++)
			_data.push_back(_op(_data.back(), i));
		return _data[n];
		}
	};

typedef RecursiveSeries<Harmonic_Rec<double, int> > Harmonic;
typedef RecursiveSeries<Harmonic_2_Rec<double, int> > Harmonic_2;

inline double undefined()
	{
	return std::numeric_limits<double>::quiet_NaN();
	}

template<class I>
I num_pairs(I n) { return n*(n-1)/2; }

#endif	// MATHUTIL_H
