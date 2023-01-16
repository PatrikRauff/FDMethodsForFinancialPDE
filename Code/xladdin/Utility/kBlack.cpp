#include "kBlack.h"
#include "kSolver.h"
#include "kFd1d.h"

class kBlackObj : public kSolverObjective
{
public:

	//	constructor
	kBlackObj(
		double	expiry,
		double	strike,
		double	price,
		double	forward)
		: kSolverObjective(),
		myExpiry(expiry),
		myStrike(strike),
		myPrice(price),
		myForward(forward)
	{}

	//	value
	virtual double	value(double x)
	{
		double res = kBlack::call(myExpiry, myStrike, myForward, x) - myPrice;

		//	done
		return res;
	}

	//	deriv
	virtual double	deriv(double x)
	{
		double res = kBlack::vega(myExpiry, myStrike, myForward, x);

		//	done
		return res;
	}

	//	private parts
private:

	//	expiry
	double	myExpiry;
	double	myStrike;
	double	myPrice;
	double	myForward;

};

// implied vol
double
kBlack::implied(
	double	expiry,
	double	strike,
	double	price,
	double	forward)
{
	//	calc intrinsic
	double intrinc = max(forward - strike,0.0);
	if (price <= intrinc) return 0.0;

	//	objective
	kBlackObj obj(expiry, strike, price, forward);

	//	start guess
	double volatility = 0.1;
	int    numIter = 10;
	double epsilon = (price - intrinc) * kConstants::epsilon();

	//	solve
	kSolver::newtonRapson(obj, volatility, numIter, epsilon, nullptr);

	//	bound
	volatility = max(0.0, volatility);

	//	done
	return volatility;
}

/// From Bachlier FD Runner ////
//	fd runner
bool kBlack::fdRunner (
	const double		s0,
	const double		r,
	const double		mu,
	const double		sigma,
	const double		expiry,
	const double		strike,
	const bool			dig,		// digital
	const int			pc,			//	put (-1) call (1)
	const int			ea,			//	european (0), american (1)
	const int			smooth,		//	smoothing
	const double		theta,
	const int			wind,
	const double		numStd,
	const int			numT,
	const int			numS,
	const bool			update,
	const int			numPr,
	double& res0,
	kVector<double>& s,
	//kVector<double>& time , 
	kVector<double>& res,
	kVector<double>& exBound,
	string& error)
{
	//	helps
	int h, i, p ;
	double x0;
	
	//Log S0 before grid
	x0 = log(s0);
	kVector<double> x;

	//	construct s axis
	/*
	Making the S grid initially for the logged spot (Log{S_0}).
	This seems to create more stable results, since it it
	possible to distance the grid exponentially equal from step to step
	*/
	double t = max(0.0, expiry);
	double std = sigma * sqrt(t);
	double xlow = x0 - numStd * std;
	double xhigh = x0 + numStd * std;
	int    nums = 2 * (max(0, numS) / 2 + 1) ;
	int div = nums - 2;
	if (numS <= 1 || xlow == xhigh)
	{
		nums = 1;
		div = 1;
	}
	double dx = (xhigh - xlow) / div ;
	double ds = exp(dx);
	s.resize(nums);
	x.resize(nums);
	x(0) = xlow;
	s(0) = exp(xlow);
	for (i = 1; i < nums; ++i)
	{
		x(i) = x(i - 1) + dx;
		s(i) = s(i - 1) * ds;
	}

	//	construct fd grid
	kFd1d<double> fd;
	fd.init(1, x, true);

	/*
	Help American Condition - maybe useful maybe not - Normal is call.
	Reverse is Put.
	Intuition: We want for Call/Put the Smallest/Largest value for direct exercise.
	Hence we loop from largest to smallest/Smallest to largest -  starting ITM
	till we reach the limit for beaing OTM.
	*/
	int start, stop, iter;
	start = nums - 1;
	stop = -1;
	iter = -1;
	if (pc < 0)
	{
		start = 0;
		stop = nums - 1;
		iter = 1;
	}

	//	set terminal result
	double sl, su;
	res.resize(nums);
	for (i = 0; i < nums; ++i)
	{
		if (smooth == 0 || i == 0 || i == nums - 1)
		{
			if (dig) res(i) = 0.5 * (kInlines::sign(s(i) - strike) + 1.0);
			else    res(i) = max(0.0, s(i) - strike);
		}
		else
		{
			sl = 0.5 * (s(i - 1) + s(i));
			su = 0.5 * (s(i) + s(i + 1));
			if (dig) res(i) = kFiniteDifference::smoothDigital(sl, su, strike);
			else	 res(i) = kFiniteDifference::smoothCall(sl, su, strike);
		}

		if (pc < 0)
		{
			if (dig) res(i) = 1.0 - res(i);
			else    res(i) -= (s(i) - strike);
		}
	}

	//	time steps
	int    numt = max(0, numT);
	double dt = t / max(1, numt);

	//Define ExBound and Time grid:
	exBound.resize(numt + 1 , 0.0);
	exBound(numt) = strike;

	//	repeat
	int nump = max(1, numPr);
	for (p = 0; p < nump; ++p)
	{
		//	set parameters
		for (i = 0; i < nums; ++i)
		{
			fd.r()(i) = r;
			fd.mu()(i) = mu;
			fd.var()(i) = sigma * sigma;
		}
		//	roll
		fd.res()(0) = res;
		for (h = numt - 1; h >= 0; --h) // Loops over number of time steps. 
		{
			//time(h) = t - (h) * dt;
			fd.rollBwd(dt, update || h == (numt - 1), theta, wind, fd.res()); //Solves Backward PDE System for
			if (ea > 0)	
			{	
				//Maybe useful - Not that many If statements outside the double loop - maybe more efficient?
				for (i = start; i != stop ; i += iter)
				{
					if (res(i) > abs( fd.res()(0)(i) ) )
					{
						fd.res()(0)(i) = res(i);
						exBound(h) = s(i); //s(i)5
					}
				}
			}
		}
	}

	//	set result
	res = fd.res()(0);
	res0 = fd.res()(0)(nums / 2 - 1);
	//	done
	return true;
}

//	fd runner
bool kBlack::fdRunnerBarrier(
	const double		s0,
	const double		r,
	const double		mu,
	const double		sigma,
	const double		expiry,
	const double		strike,
	const double		barrier,
	const bool			dig,		// digital
	const int			pc,			//	put (-1) call (1)
	const int			ea,			//	european (0), american (1)
	const int			smooth,		//	smoothing
	const double		theta,
	const int			wind,
	const double		numStd,
	const int			numT,
	const int			numS,
	const bool			update,
	const int			numPr,
	double& res0,
	kVector<double>& s,
	//kVector<double>& time , 
	kVector<double>& res,
	kVector<double>& exBound,
	string& error)
{
	//	helps
	int h, i, p;
	double x0;

	//Log S0 before grid
	x0 = log(s0);
	kVector<double> x;

	//	construct s axis
	/*
	Making the S grid initially for the logged spot (Log{S_0}).
	This seems to create more stable results, since it it
	possible to distance the grid exponentially equal from step to step
	*/
	double t = max(0.0, expiry);
	double std = sigma * sqrt(t);
	double xlow = x0 - numStd * std;
	double xhigh = x0 + numStd * std;
	int    nums = 2 * (max(0, numS) / 2 + 1);
	int div = nums - 2;
	if (numS <= 1 || xlow == xhigh)
	{
		nums = 1;
	}
	double dx = (xhigh - xlow) / div;
	double ds = exp(dx);
	s.resize(nums);
	x.resize(nums);
	x(0) = xlow;
	s(0) = exp(xlow);
	for (i = 1; i < nums; ++i)
	{
		x(i) = x(i - 1) + dx;
		s(i) = s(i - 1) * ds;
	}

	//	construct fd grid
	kFd1d<double> fd;
	fd.init(1, x, true);

	/*
	Help American Condition - maybe useful maybe not - Normal is call.
	Reverse is Put.
	Intuition: We want for Call/Put the Smallest/Largest value for direct exercise.
	Hence we loop from largest to smallest/Smallest to largest -  starting ITM
	till we reach the limit for beaing OTM.
	*/
	int start, stop, iter;
	start = nums - 1;
	stop = -1;
	iter = -1;
	if (pc < 0)
	{
		start = 0;
		stop = nums - 1;
		iter = 1;
	}

	//	set terminal result
	double sl, su;
	res.resize(nums);
	for (i = 0; i < nums; ++i)
	{
		if (smooth == 0 || i == 0 || i == nums - 1)
		{
			if (dig) res(i) = 0.5 * (kInlines::sign(s(i) - strike) + 1.0);
			else    res(i) = max(0.0, s(i) - strike);
		}
		else
		{
			sl = 0.5 * (s(i - 1) + s(i));
			su = 0.5 * (s(i) + s(i + 1));
			if (dig) res(i) = kFiniteDifference::smoothDigital(sl, su, strike);
			else	 res(i) = kFiniteDifference::smoothCall(sl, su, strike);
		}

		if (pc < 0)
		{
			if (dig) res(i) = 1.0 - res(i);
			else    res(i) -= (s(i) - strike);
		}
	}

	//	time steps
	int    numt = max(0, numT);
	double dt = t / max(1, numt);

	//Define ExBound and Time grid:
	exBound.resize(numt + 1, 0.0);
	exBound(numt) = strike;

	//	repeat
	int nump = max(1, numPr);
	for (p = 0; p < nump; ++p)
	{
		//	set parameters
		for (i = 0; i < nums; ++i)
		{
			fd.r()(i) = r;
			fd.mu()(i) = mu;
			fd.var()(i) = sigma * sigma;
		}
		//	roll
		fd.res()(0) = res;
		for (h = numt - 1; h >= 0; --h) // Loops over number of time steps. 
		{
			//time(h) = t - (h) * dt;
			fd.rollBwd(dt, update || h == (numt - 1), theta, wind, fd.res()); //Solves Backward PDE System for
			if (ea > 0)
			{
				//Maybe useful - Not that many If statements outside the double loop - maybe more efficient?
				for (i = start; i != stop; i += iter)
				{
					if (res(i) > abs(fd.res()(0)(i)))
					{
						fd.res()(0)(i) = res(i);
						exBound(h) = s(i); //s(i)5
					}
				}
			}
			for (i = 0; i < nums; ++i)
			{
				if (x(i) <= log(barrier)) fd.res()(0)(i) = 0;
			}
		}
	}

	//	set result
	res = fd.res()(0);
	res0 = fd.res()(0)(nums / 2 - 1);
	//	done
	return true;
}


////	fd runner
//bool kBlack::fdRunnerBarrier(
//	const double		s0,
//	const double		r,
//	const double		mu,
//	const double		sigma,
//	const double		expiry,
//	const double		strike,
//	const double		bar,		// barrier 
//	const bool			dig,		// digital
//	const int			pc,			//	put (-1) call (1)
//	const int			ea,			//	european (0), american (1)
//	const int			smooth,		//	smoothing
//	const double		theta,
//	const int			wind,
//	const double		numStd,
//	const int			numT,
//	const int			numS,
//	const bool			update,
//	const int			numPr,
//	double& res0,
//	kVector<double>& s,
//	kVector<double>& res,
//	kVector<double>& exBound,
//	string& error)
//{
//	//	helps
//	int h, i, p;
//	double x0;
//
//	//Log S0 before grid
//	x0 = log(s0);
//	kVector<double> x;
//
//	//	construct s axis
//	/*
//	Making the S grid initially for the logged spot (Log{S_0}).
//	This seems to create more stable results, since it it
//	possible to distance the grid exponentially equal from step to step
//	*/
//	double t = max(0.0, expiry);
//	double std = sigma * sqrt(t);
//	double xlow = x0 - numStd * std;
//	double xhigh = x0 + numStd * std;
//	int    nums = 2 * (max(0, numS) / 2 + 1);
//	//int nums = 2 * (numS / 2);
//	if (numS <= 1 || xlow == xhigh)
//	{
//		nums = 1;
//	}
//	double dx = (xhigh - xlow) / nums;
//	double ds = exp(dx);
//	s.resize(nums);
//	x.resize(nums);
//	x(0) = xlow;
//	s(0) = exp(xlow);
//	for (i = 1; i < nums; ++i)
//	{
//		x(i) = x(i - 1) + dx;
//		s(i) = s(i - 1) * ds;
//	}
//	//	construct fd grid
//	kFd1d<double> fd;
//	fd.init(1, x, true);
//
//	/*
//	Help American Condition - maybe useful maybe not - Normal is call.
//	Reverse is Put.
//	Intuition: We want for Call/Put the Smallest/Largest value for direct exercise.
//	Hence we loop from largest to smallest/Smallest to largest -  starting ITM
//	till we reach the limit for beaing OTM.
//	*/
//	int start, stop, iter;
//	start = nums - 1;
//	stop = -1;
//	iter = -1;
//	if (pc < 0)
//	{
//		start = 0;
//		stop = nums - 1;
//		iter = 1;
//	}
//
//	//	set terminal result
//	double sl, su;
//	res.resize(nums);
//	for (i = 0; i < nums; ++i)
//	{
//		if (smooth == 0 || i == 0 || i == nums - 1)
//		{
//			if (dig) res(i) = 0.5 * (kInlines::sign(s(i) - strike) + 1.0);
//			else    res(i) = max(0.0, s(i) - strike);
//		}
//		else
//		{
//			sl = 0.5 * (s(i - 1) + s(i));
//			su = 0.5 * (s(i) + s(i + 1));
//			if (dig) res(i) = kFiniteDifference::smoothDigital(sl, su, strike);
//			else	 res(i) = kFiniteDifference::smoothCall(sl, su, strike);
//		}
//
//		if (pc < 0)
//		{
//			if (dig) res(i) = 1.0 - res(i);
//			else    res(i) -= (s(i) - strike);
//		}
//	}
//
//	//	time steps
//	int    numt = max(0, numT);
//	double dt = t / max(1, numt);
//
//	//Define ExBound and Time grid:
//	exBound.resize(numt + 1, 0.0);
//	exBound(numt) = strike;
//	//time.resize(numt + 1, 0.0);
//
//	//	repeat
//	int nump = max(1, numPr);
//	for (p = 0; p < nump; ++p)
//	{
//		//	set parameters
//		for (i = 0; i < nums; ++i)
//		{
//			fd.r()(i) = r;
//			fd.mu()(i) = mu ;
//			fd.var()(i) = sigma * sigma ;
//		}
//
//		//	roll
//		fd.res()(0) = res;
//		for (h = numt - 1; h >= 0; --h) // Loops over number of time steps. 
//		{
//			fd.rollBwd(dt, update || h == (numt - 1), theta, wind, fd.res()); //Solves Backward PDE System for
//			if (ea > 0)
//			{
//				//Maybe useful - Not that many If statements outside the double loop - maybe more efficient?
//				for (i = start; i != stop; i += iter)
//				{
//					if (res(i) > abs(fd.res()(0)(i)))
//					{
//						fd.res()(0)(i) = res(i);
//						exBound(h) = s(i); //s(i)5
//					}
//				}
//			}
//			for (i = 0; i < nums; ++i)
//			{
//				if ( x(i) <= log(bar) ) fd.res()(0)(i) = 0;
//			}
//		}
//	}
//
//	//	set result
//	res = fd.res()(0);
//	res0 = fd.res()(0)(nums / 2);
//	//	done
//	return true;
//}


//	fd runner
bool kBlack::fdRunnerBarAbs(
	const double		s0,
	const double		r,
	const double		mu,
	const double		sigma,
	const double		expiry,
	const double		strike,
	const double		bar,		// Parrier 
	const bool			dig,		// digital
	const int			pc,			//	put (-1) call (1)
	const int			ea,			//	european (0), american (1)
	const int			smooth,		//	smoothing
	const double		theta,
	const int			wind,
	const double		numStd,
	const int			numT,
	const int			numS,
	const bool			update,
	const int			numPr,
	double& res0,
	kVector<double>& s,
	//kVector<double>& time,
	kVector<double>& res,
	kVector<double>& exBound,
	string& error)
{
	//	helps
	int h, i, p;
	double x0;
	bool barrier = (bar > 0.0001 && bar < s0);

	//Log S0 before grid
	x0 = log(s0);
	kVector<double> x;

	//	construct s axis
	/*
	Making the S grid initially for the logged spot (Log{S_0}).
	This seems to create more stable results, since it it
	possible to distance the grid exponentially equal from step to step
	*/
	double t = max(0.0, expiry);
	double std = sigma * sqrt(t);
	double xlow = x0 - numStd * std;
	double xhigh = x0 + numStd * std;
	int    nums = 2 * (max(0, numS) / 2 + 1);
	int div = nums - 2;
	if (numS <= 1 || xlow == xhigh)
	{
		nums = 1;
		div = 1;
	}
	double dx = (xhigh - xlow) / div;
	double ds = exp(dx);
	s.resize(nums);
	x.resize(nums);
	x(0) = xlow;
	s(0) = exp(xlow);
	for (i = 1; i < nums; ++i)
	{
		x(i) = x(i - 1) + dx;
		s(i) = s(i - 1) * ds;
	}

	// Barrier check and index the grid with the barrier
	int b = 0;
	if (barrier)
	{
		 b = (int)((log(bar) - xlow) / dx);
		s(b) = bar;
		x(b) = log(bar);
	}

	//	construct fd grid
	kFd1d<double> fd;
	fd.init(1, x, true);

	/*
	Help American Condition - maybe useful maybe not - Normal is call.
	Reverse is Put.
	Intuition: We want for Call/Put the Smallest/Largest value for direct exercise.
	Hence we loop from largest to smallest/Smallest to largest -  starting ITM
	till we reach the limit for beaing OTM.
	*/
	int start, stop, iter;
	start = nums - 1;
	stop = -1;
	iter = -1;
	if (pc < 0)
	{
		start = 0;
		stop = nums - 1;
		iter = 1;
	}

	//	set terminal result
	double sl, su;
	res.resize(nums);
	for (i = 0; i < nums; ++i)
	{
		if (smooth == 0 || i == 0 || i == nums - 1)
		{
			if (dig) res(i) = 0.5 * (kInlines::sign(s(i) - strike) + 1.0);
			else    res(i) = max(0.0, s(i) - strike);
		}
		else
		{
			sl = 0.5 * (s(i - 1) + s(i));
			su = 0.5 * (s(i) + s(i + 1));
			if (dig) res(i) = kFiniteDifference::smoothDigital(sl, su, strike);
			else	 res(i) = kFiniteDifference::smoothCall(sl, su, strike);
		}

		if (pc < 0)
		{
			if (dig) res(i) = 1.0 - res(i);
			else    res(i) -= (s(i) - strike);
		}
	}

	//	time steps
	int    numt = max(0, numT);
	double dt = t / max(1, numt);

	//Define ExBound and Time grid:
	exBound.resize(numt + 1, 0.0);
	exBound(numt) = strike;
	//time.resize(numt + 1, 0.0);

	//	repeat
	int nump = max(1, numPr);
	for (p = 0; p < nump; ++p)
	{
		//	set parameters
		for (i = 0; i < nums; ++i)
		{
			fd.r()(i) = r;
			fd.mu()(i) = mu;
			fd.var()(i) = sigma * sigma;
		}

		if (barrier)
		{
			fd.var()(b) = 0;
			fd.mu()(b) = 0;
		}

		//	roll
		fd.res()(0) = res;
		for (h = numt - 1; h >= 0; --h) // Loops over number of time steps. 
		{
			//time(h) = t - (h)*dt;
			fd.rollBwd(dt, update || h == (numt - 1), theta, wind, fd.res()); //Solves Backward PDE System for
			if (ea > 0)
			{
				//Maybe useful - Not that many If statements outside the double loop - maybe more efficient?
				for (i = start; i != stop; i += iter)
				{
					if (res(i) > abs(fd.res()(0)(i)))
					{
						fd.res()(0)(i) = res(i);
						exBound(h) = s(i); //s(i)5
					}
				}
			}
		}
	}

	//	set result
	res = fd.res()(0);
	res0 = fd.res()(0)(nums / 2 - 1);
	//	done
	return true;
}

/// From Bachlier FD Runner ////
//	fd runner
bool kBlack::fdFwdRunner(
	const double		s0,
	const double		r,
	const double		mu,
	const double		sigma,
	const double		expiry,
	const double		strike,
	const bool			dig,
	const int			pc,			//	put (-1) call (1)
	const int			ea,			//	european (0), american (1)
	const int			smooth,		//	smoothing
	const double		theta,
	const int			wind,
	const double		numStd,
	const int			numT,
	const int			numS,
	const bool			logFd,
	const bool			update,
	const int			numPr,
	kVector<double>& K,
	kVector<double>& s,
	//kVector<double>& time,
	kVector<kVector<double>>& res,
	kVector<kVector<double>>& prob,
	string& error)
{
	//	helps
	int h, i, p , j;
	double x0;

	//Log S0 before grid
	x0 = log(s0);
	kVector<double> x , payoff;

	//	construct s axis
	/*
	Making the S grid initially for the logged spot (Log{S_0}).
	This seems to create more stable results, since it it
	possible to distance the grid exponentially equal from step to step
	*/
	double t = max(0.0, expiry);
	double std = sigma * sqrt(t);
	double xlow = x0 - numStd * std;
	double xhigh = x0 + numStd * std;
	int    nums = 2 * (max(0, numS) / 2 + 1);
	int div = nums - 2;
	if (numS <= 1 || xlow == xhigh)
	{
		nums = 1;
		div = 1;
	}
	double dx = (xhigh - xlow) / div;
	double ds = exp(dx);
	s.resize(nums);
	x.resize(nums);
	x(0) = xlow;
	s(0) = exp(xlow);
	for (i = 1; i < nums; ++i)
	{
		x(i) = x(i - 1) + dx;
		s(i) = s(i - 1) * ds;
	}

	//	construct fd grid
	kFd1d<double> fd;
	if (logFd) {
		fd.init(1, x, true);

	}
	else
	{
		fd.init(1, s, false);
	}


	kVector<double> tempProb;
	//	set initial probabilities
	tempProb.resize(nums , 0.);
	tempProb(nums / 2) = 1.;

	//	time steps
	int    numt = max(0, numT);
	double dt = t / max(1, numt);

	//Define ExBound and Time grid:
	//time.resize(numt + 1, 0.0);
	res.resize(numt + 1); 
	prob.resize(numt + 1); //res containts numt + 1 of transition
	//	repeat
	int nump = max(1, numPr);
	for (p = 0; p < nump; ++p)
	{
		//	set parameters
		if (logFd) {
			for (i = 0; i < nums; ++i)
			{
				fd.r()(i) = r;
				fd.mu()(i) = mu ;
				fd.var()(i) = sigma * sigma ;
			}

		}
		else
		{
			for (i = 0; i < nums; ++i)
			{
				fd.r()(i) = r;
				fd.mu()(i) = mu * s(i);
				fd.var()(i) = sigma * sigma * s(i) * s(i);
			}
		}
		//	roll
		fd.res()(0) = tempProb;
		prob(0) = tempProb;
		for (h = 0; h  < numt ; ++h) // Loops over number of time steps. 
		{
			//time(h + 1) =  h*dt + dt;
			fd.rollFwd(dt, update || h == 0, theta, wind, fd.res()); //Solves Backward PDE System for
			prob(h+1) = fd.res()(0);
		}
	}

	//	set result from transition probabiliies
	kVector<double> tempRes;
	kMatrix<double> PO;
	PO.resize(nums, nums);
	tempRes.resize(nums, 0.);
	double sl, su;

	for (i = 0; i < nums; ++i)
	{
		for (j = i; j < nums; ++j)
		{
			if (smooth == 0 || i == 0 || i == nums - 1)
			{
				if (dig) PO(j, i) = 0.5 * (kInlines::sign(s(j) - s(i)) + 1.0);
				else    PO(j, i) = max(0.0, s(j) - s(i));
			}
			else
			{
				sl = 0.5 * (s(j - 1) + s(j));
				su = 0.5 * (s(j) + s(j + 1));
				if (dig) PO(j, i)= kFiniteDifference::smoothDigital(sl, su, s(i));
				else	 PO(j, i) = kFiniteDifference::smoothCall(sl, su, s(i));
			}

			if (pc < 0)
			{
				if (dig) PO(j, i) = 1.0 - PO(j, i);
				else    PO(j, i) -= (s(j) - s(i));
			}
		}
	}

	//Compute prices for (tau , K) Vector of Tau containing vector for K
	for (h = 0; h <= numt ; ++h) // time
	{
		res(h) = tempRes; // get size into vector(vector)
		for (j = 0; j < nums; ++j) // K grid
		{
			for (i = 0; i < nums; ++i) // underlying
			{
				res(h)(j) += prob(h)(i) * PO(i,j);
			}
		}
	}
	
	//prob = fd.res()(0);
	//	done
	return true;
}