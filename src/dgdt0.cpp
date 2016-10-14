/*
 * evolution of clique size distribution p_k
 * truncated at max k
 * 
 * NOTE this version is for unstructured model only
 */


// [[Rcpp::depends(RcppArmadillo, BH)]]

#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <RcppArmadillo.h>

static const double MINY = 1e-12;

using namespace arma;
using namespace Rcpp; 
using namespace std; 

typedef std::vector<double> state_type; 

// globals
vec g_Fs, g_Ys; 
double g_Ne; 
double g_hres; 
double g_treeT; 
int g_maxk; 
double g_A; 


////////////////////////////////////////////////////////////////////////
// solve change in psi vec (transition probabilities), not counting accrual
// this is same as Q for the p_uk, but includes extra decay with Ne
class ODE_dgdt0{ 
public: 
	vec g0;
	double Lambda0; 
public:
	ODE_dgdt0(vec g0_, double Lambda0_) : g0(g0_),Lambda0(Lambda0_) {};
	void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
	{
		int i_t =  (int)std::max(0., std::min( g_hres * t / g_treeT , (double)(g_hres-1.))); 
		double f = g_Fs(i_t); 
		double y = g_Ys(i_t); 
				
		int k, i, j; 
		
		vec g2 = zeros<vec>(g_maxk); 
		vec ug2 = zeros<vec>(g_maxk); 
		//~ vec g(x); //  not quite 
		vec g = zeros<vec>(g_maxk); 
		vec ug = zeros<vec>(g_maxk); 
		for (k = 1; k <= g_maxk; k++){
			g(k-1) = max(0., x.at(k-1)); 
			ug(k-1) = max(0., x.at(k-1)); 
		}
		g = normalise( g, 1.);
		// B= A / g'(1)
		double m1 = 0.; 
		for (k = 1; k <= g_maxk; k++){
			m1 += g(k-1) * k; 
		}
		double B = g_A / m1;
		//~ double B = min( g_A / m1, y);
		
		// (dg / dt) = (g2 - g) * B * f / Y2
		for (i = 1; i <= g_maxk; i++){
			for ( j = 1; j <= g_maxk; j++){
				k = min( g_maxk, i + j); 
				g2(k-1) += g(i-1)*g(j-1);
				ug2(k-1) += ug(i-1)*ug(j-1);
			}
		}
		
		//~ double tr = B * f / (y*y);
		double tr = max(0.,(B-1.)) * f / (y*y);
		for (k = 1; k <= g_maxk; k++){
			if (x.at(k-1) < 0){
				dxdt.at(k - 1)= max(0., (g2(k-1) - g(k-1)) * tr); // 
			} else{
				dxdt.at(k - 1)=  (g2(k-1) - g(k-1)) * tr; // 
			}
			//~ dxdt.at(k - 1)= (g2(k-1) - g(k-1)) * tr - g(k-1) * (k*(k-1))/(2.*g_Ne); 
			//~ dxdt.at(k - 1) = (ug2(k-1) - ug(k-1)) * tr - ug(k-1) * (k*(k-1))/(2.*g_Ne); 
		}
		
		// d Lambda/dt = ...
		// lambda = (A/(2*Ne)) * g'' / g'
		double dLdt = 0.; 
		for (k = 2; k <= g_maxk ; k++){
			dLdt += B * g(k-1) * k * (k-1) / (2. * g_Ne);  // remember *B!
		}
		dxdt.at( i_Lambda() ) = dLdt; 
	}
	state_type get_statevec0(){
		state_type x(g_maxk+1, 0.); 
		for (int i = 0; i < g_maxk; i++){
			x.at(i) = g0.at(i); 
		}
		x.at(g_maxk) = Lambda0; 
		return x;
	}
	double get_Lambda(state_type x){
		return x.at(i_Lambda());
	}
	vec get_g(state_type x){
		vec g = zeros<vec>(g_maxk); 
		for (int i = 0; i < g_maxk; i++){
			g(i) = x.at(i); 
		}
		return g;
	}
private:
	int i_Lambda(){
		//~ return g_maxk-1; 
		return g_maxk; 
	}
};

//[[Rcpp::export()]]
List solve_dgdt0(vec times, vec Fs, vec Ys
 , vec g0 // initial state vec
 , double Lambda0 // initial cum hazard
 , double h0
 , double h1
 , double A
 , double Ne
 , double treeT // actually length of time axis //TODO deprecate
){
	::g_Fs = Fs; 
	::g_Ys = Ys; 
	::g_hres = times.size();
	//~ ::g_treeT = treeT; 
	::g_Ne = Ne; 
	::g_maxk = g0.size(); 
	::g_A = A; 
	::g_treeT = times(0)-times(times.size()-1);
	
	ODE_dgdt0 dgdt0(g0, Lambda0); 
	state_type x = dgdt0.get_statevec0();
	size_t steps = boost::numeric::odeint::integrate( dgdt0 ,  x, h0 , h1 , (h1-h0)/100. );  
	
	List o; 
	o["Lambda"] = dgdt0.get_Lambda(x);
	
	vec g = dgdt0.get_g(x);
	for (int k = 1; k <= g_maxk; k++){
		g(k-1) = max(0., x.at(k-1)); 
	}
	g = normalise( g, 1.);
	o["g"] = g;
	return o; 
}



////////////////////////////////////////////////////////////////////////
// CoM12: simple ODE solver to find cumulative coalescent hazard in each interval
class ODE_CoM12L{ 
public: 
	double Lambda0; 
public:
	ODE_CoM12L(double Lambda0_) : Lambda0(Lambda0_) {};
	void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
	{
		int i_t =  (int)std::max(0., std::min( g_hres * t / g_treeT , (double)(g_hres-1.))); 
		double f = g_Fs[i_t]; 
		double y = g_Ys[i_t]; 
		
		//~ dxdt.at(0) = (g_A * (g_A-1.)/2.) * 2. * f / ( y * y ); 
		dxdt.at(0) = (g_A * (g_A-1.)/2.) * 2. * f / std::max(1., ( y * y- y )); 
	}
	state_type get_statevec0(){
		state_type x(1, 0.); 
		x.at(0) = Lambda0; 
		return x;
	}
	double get_Lambda(state_type x){
		return x.at(0);
	}
};

//[[Rcpp::export()]]
double solve_CoM12L(vec times, vec Fs, vec Ys
 , double Lambda0 // initial cum hazard
 , double h0
 , double h1
 , double A
){
	::g_Fs = Fs; 
	::g_Ys = Ys; 
	::g_hres = times.size();
	::g_A = A; 
	::g_treeT = times(0)-times(times.size()-1); 
	
	ODE_CoM12L com12L( Lambda0); 
	state_type x = com12L.get_statevec0();
	size_t steps = boost::numeric::odeint::integrate( com12L ,  x, h0 , h1 , (h1-h0)/100. );  
	
	return com12L.get_Lambda(x);
}
