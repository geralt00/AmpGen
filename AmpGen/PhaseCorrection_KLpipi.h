#ifndef AMPGEN_PHASECORRECTION_KLPIPI_H
#define AMPGEN_PHASECORRECTION_KLPIPI_H

#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Expression.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MyStructs.h"


#include <string>
#include <map>
#include <vector>

namespace AmpGen
{

	class PhaseCorrection_KLpipi
	{
		using CE = CompiledExpression<real_t(const real_t*, const real_t*)>;

		public:
			explicit PhaseCorrection_KLpipi(); 
			explicit PhaseCorrection_KLpipi(MinuitParameterSet &MPS);

			// create an expression for 1 term in polynomial phase correction, with i and j powers of each coordinate
			Expression polynomial(EventType& type, size_t i, size_t j);
			// return (s-mu)/sigma all squared for exponent of gaussian function
			real_t gaussianExponent(real_t& s, real_t& mu, real_t& sigma);
			// formula for a single bias
			real_t bias(dalitzPair<real_t>& coords, size_t& index);


			// evaluate the phase correction at a given event
			real_t eval(const Event& event);
			// as above for the gaussian bias rather than a polynomial correction
			real_t evalBias(Event& event);


			// two options for coordinate systems to use, calculating the normal s12/s13 + transforming them
			// std::pair<Expression, Expression> baseDalitzCoords(EventType& eventType);
			dalitzPair<Expression>  transformedDalitzCoords(EventType& eventType);
			dalitzPair<Expression>  squareDalitzCoords(EventType& eventType);

			// compile the CompiledExpressions for each of the P_i_j needed
			void compilePolynomialExpressions(EventType& eventType);

			// update the coefficients vector from the MPS
			void updateCoeffs(MinuitParameterSet& MPS);

			bool doPolynomial();
			bool doBias();

		private:
			real_t TEMPdeltaC_;
			size_t order_; // order of polynomial correction
			std::string correctionType_; // note can only cope with legendre at the moment
			Expression legendre(Expression& x, real_t n); // separated from polynomial so that more types of polys can be added if needed
			bool doPolynomial_{false}; // are we doing a polynomial phase correction? ( any order > 0)
			bool doBias_{false}; // are we doing a gaussian bias? ( signle or double )

			// Equivalent vectors to store values for CE's of the terms P_i_j, coeffs Ci_j, values of i and j			
			std::vector<CE> compiledExpressions_;
			std::vector<real_t> coefficients_;
			//Be careful with these, they are the indices of the coefficients in the MPS, not the actual values of i and j
			//Will they share the same index as the coefficients? I think so, but not sure
			std::vector<size_t> iIndices_;
			std::vector<size_t> jIndices_;
			size_t nTerms_{0};
			size_t nBias_{0};

			// coeffs to do with the gaussian bias - pair so first is the f1 and second f2 if using a double
			dalitzPair<real_t> emptyPair{};
			std::vector<dalitzPair<real_t>> mu_ = {emptyPair, emptyPair};
			std::vector<dalitzPair<real_t>> sigma_= {emptyPair, emptyPair};
			std::vector<real_t> epsilon_= {0, 0};
			std::vector<real_t> A_= {0, 0};

			// Jake manually found numbers to scale dalitz coords to be in range 0->1:
			real_t c1_ = -3.1171885586526695;
			real_t m1_ = 2.23407421671132946;
			real_t c2_ = -9.54231895051727e-05;
			real_t m2_ = 0.8051636393861085;
	};

}

#endif