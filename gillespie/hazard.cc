/*
 * hazard.cc
 *
 *  Created on: Jan 9, 2013
 *      Author: tel
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <cmath>
#include "hazard.hh"

using namespace boost::numeric::ublas;

void Hazard::Update(state_type M) {
	static const matrix<double> summer (scalar_matrix<double> (1, H.size1(), 1));
	for (int i = 0; i<Hfunc.size1(); ++i) {
		H(i,0) = Hfunc(i,0)(M);
	}
	H0 = std::abs(prod(summer, H)(0,0));
}
void Hazard::operator() ( const state_type &x , state_type &dxdt , const double /* t */ ) {
	double temp;
	for (int i = 0; i<x.size1(); ++i){
		temp = 0;
		for (int j = 0; j<S.size2(); ++j) {
			temp += S(i,j)*Hfunc(j,0)(x);
		}
		dxdt(i,0) = temp;
	}
}

void Hazard::InitHfunc() {
	Hfunc(0,0) = [this] (const state_type &M) {return (this->c)(0,0)*M(0,0);};
	Hfunc(1,0) = [this] (const state_type &M) {return (this->c)(1,0)*M(1,0);};
}
