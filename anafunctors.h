#ifndef ANAFUNCTORS_H
#define ANAFUNCTORS_H

#include <cmath>

template<class SAMPLE>
struct AnalysisBase
	{
	typedef SAMPLE sample_t;

	virtual double analyse(const SAMPLE & s) const = 0;
	double operator()(const SAMPLE & s) const 
		{ return this->analyse(s); }
	};

template<class SAMPLE, class T>
struct AnalysisWrapper
	{
	typedef SAMPLE sample_t;

	typedef T (sample_t::* func_t)() const;

	template<func_t FUNC>
	struct Analysis : public AnalysisBase<SAMPLE>
		{
		double analyse(const SAMPLE & s) const
			{
			return (s.*FUNC)();
			}
		};

	template<func_t FUNC>
	static AnalysisBase<SAMPLE> * create()
		{
		return new Analysis<FUNC>();
		}
	};

struct Aggregate
	{
protected:
	size_t _n;
	double _mean, _sqr_sum, _std, _skew_sum, _kurt_sum, _skew, _kurt;

public:
	Aggregate()
		: _n(0), _mean(0), _sqr_sum(0), _std(0), _skew_sum(0), _kurt_sum(0)
		{}

	void operator()(double value)
		{
		if (value != value)
			return;

		_n++;
		_mean += value;
		_sqr_sum += value * value;
		//_skew_sum += pow(value,3);
		_kurt_sum += pow((value - _mean),4);


		_skew_sum += pow((value - _mean),3);

		}

	void analyse()
		{
		if (_n<2)
			return;
		
		const double ssq = _sqr_sum - _mean*_mean/_n;


		_std = sqrt(ssq / (_n-1));
		
		_mean /= _n;
		_skew = _skew_sum / _n;

		_kurt = (_kurt_sum / _n) - 3;

		}

	double mean() const
		{ return _mean; }
	double std() const
		{ return _std; }
	double skew() const
		{ return _skew; }
	double kurt() const
		{ return _kurt; }
	};


#endif	// ANAFUNCTORS_H
