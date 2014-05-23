#pragma warning ( disable : 4503 )
#include "dynamo_solver.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;
using namespace DRDSP;

DynamoSolver::DynamoSolver( const Dynamo& dynamo ) : Dynamo(dynamo) {
	double nTime = 6000;
	double hsn0 = ::sinh(eta0);
	dtMax = 1.0 / (6.0 * nTime * hsn0*hsn0);
	Init();
}

void DynamoSolver::Advance( double dt ) {
	double dta = dt;
	uint32_t steps = 1;
	if( dt > dtMax && dtMax > 0.0 ) {
		steps = uint32_t(dt / dtMax + 1.0);
		dta = dt / steps;
	}

	for(uint32_t j=0;j<nJ;++j)
		pZero.col(j) = state.segment(j*nI,nI);
	for(uint32_t j=0;j<nJ;++j)
		tZero.col(j) = state.segment(N+j*nI,nI);
	
	npCoff( dta );
	ntCoff( dta );

	for(uint32_t i=0;i<steps;++i) {
		Step( dta );
	}

	for(uint32_t j=0;j<nJ;++j)
		state.segment(j*nI,nI) = pZero.col(j);
	for(uint32_t j=0;j<nJ;++j)
		state.segment(N+j*nI,nI) = tZero.col(j);
}

void DynamoSolver::Step( double dt ) {
	npStep( dt );
	ntStep( dt );

	pMinus = pZero;
	pZero = pPlus;
	pPlus = pMinus;

	tMinus = tZero;
	tZero = tPlus;
	tPlus = tMinus;
}

void DynamoSolver::Init() {
	state.setZero(dimension);
	pMinus.setZero(nI,nJ);
	pPlus.setZero(nI,nJ);
	pZero.setZero(nI,nJ);
	pCoef1.setZero(nI,nJ);
	pCoef2.setZero(nI,nJ);
	pCoef3.setZero(nI,nJ);
	pCoef4.setZero(nI,nJ);
	pCoef5.setZero(nI,nJ);

	tMinus.setZero(nI,nJ);
	tPlus.setZero(nI,nJ);
	tZero.setZero(nI,nJ);
	tCoef1.setZero(nI,nJ);
	tCoef2.setZero(nI,nJ);
	tCoef3.setZero(nI,nJ);
	tCoef4.setZero(nI,nJ);
	tCoef5.setZero(nI,nJ);

	tpCof1.setZero(nI,nJ);
	tpCof2.setZero(nI,nJ);
	tpCof3.setZero(nI,nJ);
	tpCof4.setZero(nI,nJ);
	tpCof5.setZero(nI,nJ);
	ptCof5.setZero(nI,nJ);

	diffA.setZero(nI,nJ);
	under.setZero(nI,nJ);
	alpha.setZero(nI,nJ);
}

void DynamoSolver::npCoff( double dTime ) {
	double Q1 = dTime / (ds*ds);
	double Q2 = dTime / (dth*dth);
	double Q3 = dTime / (2.0*ds);
	double Q4 = dTime / (2.0*dth);

	double x, th, C, C2, ax, bx, ath, bth, cxth, denom;
	for(uint32_t j=0;j<nJ;++j) {
		th = theta(j);
		for(uint32_t i=1;i<nI-1;++i) {
			x = s(i);
			C = c(i,j);
			C2 = C * C;
			ax = C2 * x * x;
			bx = x * C * ( C*(1.0 - cotheta(i)) + sinheta(i) );
			ath = C2;
			bth = -C * sintheta(j);
			cxth = -C2 / (sinheta(i)*sinheta(i));

			denom = 1.0 / ( 1.0 + 2.0*Q1*ax + 2.0*Q2*ath - 0.5*cxth*dTime );
			pCoef1(i,j) = 2.0 * (Q1*ax - bx*Q3) * denom;
			pCoef2(i,j) = 2.0 * (Q2*ath - bth*Q4) * denom;
			pCoef3(i,j) = 2.0 * (Q2*ath + bth*Q4) * denom;
			pCoef4(i,j) = 2.0 * (Q1*ax + bx*Q3) * denom;
			pCoef5(i,j) = ( 1.0 - 2.0*Q1*ax - 2.0*Q2*ath + 0.5*cxth*dTime ) * denom;
			ptCof5(i,j) = (2.0*dTime) * denom;
		}
	}
}

void DynamoSolver::ntCoff( double dTime ) {
	double Q1 = dTime / (ds * ds);
	double Q2 = dTime / (dth * dth);
	double Q3 = dTime / (2.0 * ds);
	double Q4 = dTime / (2.0 * dth);
  
	double x, th, C, c2, mult, ax, bx, pbx, ath, bth, pbth, cxth, denom;
 
	for(uint32_t j=0;j<nJ;j++) {
		th = theta(j);
		for(uint32_t i=1;i<nI-1;i++) {
			x = s(i);
			C = c(i,j);
			c2 = C * C;
			mult = -1.5 * cOmega * pi32(i,j);
			ax = c2*x*x;
			bx = x*C*( C*(1.0-cotheta(i)) + sinheta(i) );
			pbx = -x * mult * sintheta(j)*sinheta(i);
			ath = c2;
			bth = -C * sintheta(j);
			pbth = mult * ( 1.0 - cosheta(i)*costheta(j) );
			cxth = -c2 / (sinheta(i)*sinheta(i));
    
			denom = 1.0 / (1.0 + 2.0*Q1*ax + 2.0*Q2*ath - 0.5*cxth*dTime);
			under(i,j) = denom;
			tCoef1(i,j) = 2.0*(Q1*ax - bx*Q3)*denom;
			tCoef2(i,j) = 2.0*(Q2*ath - bth*Q4)*denom;
			tCoef3(i,j) = 2.0*(Q2*ath + bth*Q4)*denom;
			tCoef4(i,j) = 2.0*(Q1*ax + bx*Q3)*denom;
			tCoef5(i,j) = (1.0 - 2.0*Q1*ax - 2.0*Q2*ath + 0.5*cxth*dTime) * denom;
			tpCof1(i,j) = -2.0*denom*Q3*pbx;
			tpCof2(i,j) = -2.0*denom*Q4*pbth;
			tpCof3(i,j) = 2.0*denom*Q4*pbth;
			tpCof4(i,j) = 2.0*denom*Q3*pbx;
			tpCof5(i,j) = 2.0*dTime*denom;
		}
	}
}

void DynamoSolver::npStep( double dTime ) {

	alpha.col(0).setZero();

	for(uint32_t j=1;j<nJ;j++)
		for(uint32_t i=1;i<nI-1;i++)
			alpha(i,j) = Alpha(pZero,tZero,i,j);

	alpha.row(0).fill(alpha.row(1).mean());

	// Calculate pPlus
	for(uint32_t j=0;j<nJ;j++) {
		for(uint32_t i=1;i<nI-1;i++) {
			pPlus(i,j) = pCoef1(i,j)*pZero(i-1,j)
						+ pCoef2(i,j)*pZero(i,jm1(j))
						+ pCoef3(i,j)*pZero(i,jp1(j))
						+ pCoef4(i,j)*pZero(i+1,j)
						+ pCoef5(i,j)*pMinus(i,j)
						+ ptCof5(i,j)*tZero(i,j)*alpha(i,j);
		}
	}

	// surface boundary condition 
	pPlus.row(nI-1) = pPlus.row(nI-2) * Bound1.transpose() + pPlus.row(nI-3) * Bound2.transpose();

	// boundary condition at centre
	pPlus.row(0).fill(pPlus.row(1).mean());

	// calculate the diffusion operator on a
	diffA = ( pPlus - pMinus - alpha.cwiseProduct(tZero) * (2.0*dTime) ).cwiseProduct(under);
}

void DynamoSolver::ntStep( double dTime ) {

	// Calculate tPlus
	for(uint32_t j=0;j<nJ;j++) {
		for(uint32_t i=1;i<nI-1;i++) {
			tPlus(i,j) = tCoef1(i,j)*tZero(i-1,j)
						+ tCoef2(i,j)*tZero(i,jm1(j))
						+ tCoef3(i,j)*tZero(i,jp1(j))
						+ tCoef4(i,j)*tZero(i+1,j)
						+ tCoef5(i,j)*tMinus(i,j)
						+ tpCof1(i,j)*pZero(i-1,j)
						+ tpCof2(i,j)*pZero(i,jm1(j))
						+ tpCof3(i,j)*pZero(i,jp1(j))
						+ tpCof4(i,j)*pZero(i+1,j)
						- alpha(i,j)*diffA(i,j);
			// non-linear term
			tPlus(i,j) -= c(i,j)*( Btheta(pZero,i,j)*s(i)*dsp(alpha,i,j) + Beta(pZero,i,j)*dtp(alpha,i,j) ) * under(i,j)*2.0*dTime;
		}
	}

	// surface boundary condition 
	tPlus.row(nI-1).setZero();

	// boundary condition at centre
	tPlus.row(0).fill(tPlus.row(1).mean());
}

double DynamoSolver::dsp( const MatrixXd& p, uint32_t i, uint32_t j ) const {
	return ( p(i+1,j) - p(i-1,j) ) / ( 2.0*ds );
}

double DynamoSolver::dtp( const MatrixXd& p, uint32_t i, uint32_t j ) const {
	return ( p(i,jp1(j)) - p(i,jm1(j)) ) / ( 2.0*dth );
}

double DynamoSolver::Alpha( const MatrixXd& a, const MatrixXd& b, uint32_t i, uint32_t j ) const {
	return ( cAlpha * sintheta(j) ) / ( 1.0 + alphaB * NormB2(a,b,i,j) );
}

double DynamoSolver::Beta( const MatrixXd& a, uint32_t i, uint32_t j ) const {
	return c(i,j)*dtp(a,i,j) - a(i,j)*sintheta(j);
}

double DynamoSolver::Btheta( const MatrixXd& a, uint32_t i, uint32_t j ) const {
	return a(i,j) * (cosheta(i)*costheta(j)-1.0)/sinheta(i) + c(i,j)*s(i) * dsp(a,i,j);
}

double DynamoSolver::NormB2( const MatrixXd& a, const MatrixXd& b, uint32_t i, uint32_t j ) const {
	double b1 = Beta(a,i,j);
	double b2 = Btheta(a,i,j);
	double b3 = b(i,j);
	return b1*b1 + b2*b2 + b3*b3;
}
