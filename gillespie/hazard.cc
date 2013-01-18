/*
 * hazard.cc
 *
 *  Created on: Jan 9, 2013
 *      Author: tel
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/range/irange.hpp>
#include <iostream>
#include <cmath>
#include "hazard.hh"

using namespace boost::numeric::ublas;

typedef matrix<double> state_type;
typedef matrix<int>::iterator1 i1_t;
typedef matrix<int>::iterator2 i2_t;

void Hazard::Update(state_type M) {
	static const matrix<double> summer (scalar_matrix<double> (1, H.size1(), 1));
	for (int i = 0; i<Hfunc.size1(); ++i) {
		//std::cout << i << "\t" << H << "\t"  << Hfunc.size1() << "\t" << Hfunc.size2() << std::endl;
		std::cout << "pop H" << std::endl;
		H(i,0) = Hfunc(i,0)(M);
		std::cout << "push H" << std::endl;
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

matrix<std::function<double(matrix<double>)>> Hazard::InitHfunc() {
	matrix<std::function<double(state_type)>> hfunc(c.size1(), 1);
	double c_coeff;
	double rl = 1;
	for (matrix<int>::const_iterator1 it1 = Pre.begin1(); it1!=Pre.end1(); ++it1) {
		c_coeff = std::accumulate(it1.begin(), it1.end(), 1.0, [] (double c, int s) {
			if (s > 1) {
				for (int x: boost::irange(2, s+1))
					c *= x;
				return c;
			}
		});
		hfunc(it1.index1(), 0) = [=, this, &c_coeff, &rl] (state_type M) {
			static double c = (this->c)((int)it1.index1(), 0)*(1/c_coeff);
			rl = 1;
			for (matrix<int>::const_iterator2 it2 = it1.begin(); it2!=it1.end(); ++it2) {
				if (*it2 > 0) {
					boost::integer_range<int> ir = boost::irange(0, *it2);
					rl *= std::accumulate(ir.begin(), ir.end(), 1,[=, &M] (double rrl, int r) {
						std::cout << it1.index1() << "\t" << it2.index2() << "\t" << M << std::endl;
						std::cout << "butter" << std::endl;
						std::cout << M(it2.index2(), 0) << std::endl;
						return rrl*(M(it2.index2(), 0) - r);
					});
				}
			}
			return c*rl;
		};
	}
	return hfunc;
}
