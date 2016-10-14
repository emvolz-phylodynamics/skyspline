// [[Rcpp::depends(RcppArmadillo, BH)]]

#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <RcppArmadillo.h>

static const double MINY = 1e-12;

using namespace arma;
using namespace Rcpp; 
using namespace std; 
using namespace boost::numeric::odeint;

typedef std::vector<double> state_type; 

// globals
vec g_betas; 
double g_gamma;

////////////////////////////////////////////////////////////////////////
// solve change in psi vec (transition probabilities), not counting accrual
// this is same as Q for the p_uk, but includes extra decay with Ne
class ODE_sp0{ 
private: 
	double t1;
	double t0;
	double dt; 
	int res;
public:
	ODE_sp0(double t0_, double t1_, int res_) : t0(t0_),t1(t1_),res(res_) {
		dt = (t1 - t0) / (double)(res-1); 
	};
	void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
	{
		int i_t =  (int)std::max(0., std::min( (double)res * (t-t0) / (t1-t0) , (double)(res-1))); 
		double beta = g_betas(i_t); 
		dxdt[0] = x[0] * ( beta - g_gamma);
	}
};


struct push_back_state_and_time
{
    std::vector< double >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< double > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x[0] );
        m_times.push_back( t );
    }
};

//~ solve_semiPar0(times,x0[1], betas, theta['gamma'] )
//[[Rcpp::export()]]
vec solve_semiPar0( double t0, double t1, int res, double y0, vec betas, double gamma, double eps_abs=1.0e-10, double eps_rel =1.0e-6)
{
	::g_betas = betas; 
	::g_gamma = gamma;
	
	double dt = (t1-t0)/(res-1); 
	
	ODE_sp0 sp0(t0, t1, res);
	std::vector<double> times; times.reserve(res);
	std::vector<double> Ys; Ys.reserve(res);
	push_back_state_and_time pbsat( Ys, times);
	
	//~ runge_kutta_fehlberg78
	state_type x(1, y0);
	size_t steps = integrate_const( make_controlled( eps_abs , eps_rel , runge_kutta_fehlberg78<state_type>() )
	   , sp0 , x , t0 , t1 +dt/2., dt , pbsat );
	
	//~ List o;
	//~ o["y"] = wrap( Ys ) ;
	//~ o["times"] = wrap( times );
	//~ return o;
	return vec(Ys);
}
