/*
 * pnet.hh
 *
 *  Created on: Jan 3, 2013
 *      Author: tel
 */

#ifndef HAZARD_HH_
#define HAZARD_HH_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <functional>

using namespace boost::numeric::ublas;

class Hazard {

public:
	typedef matrix<double> state_type;

	/// constructors
	/// default
	Hazard():
		c(),
		H(),
		H0(0) {}

	/// empty, correctly dimensioned instance
	Hazard(int u, int v):
		c(v,1),
		H(v,1),
		H0(0) {}

	/// initialized from premade matrices
	Hazard(state_type Mi, state_type ci, state_type Si):
		c(ci),
		H(ci.size1(), 1),
		H0(0.0001),
		S(Si),
		Hfunc(ci.size1(), 1) {}

	const state_type c;
	const state_type S;
	state_type H;
	matrix<std::function<double(state_type)>> Hfunc;
	double H0;
	virtual void InitHfunc();
	void operator() ( const state_type &x , state_type &dxdt , const double /* t */ );
	void Update(state_type M);
};

#endif /* HAZARD_HH_ */
