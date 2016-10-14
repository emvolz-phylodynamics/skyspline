/*
 * helpers 
 */


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>


using namespace arma;
using namespace Rcpp; 
using namespace std; 


//[[Rcpp::export()]]
List update_lambda_gx_node(double A, double Ne, vec g_x, double y)
{
	int mk = g_x.size();
	double gpp = 0.;
	double gp = 0.; 
	int i,j,k;
	for (k=1; k<=mk; k++){
		gpp += g_x(k-1) * k * (k-1); 
		gp += g_x(k-1) * k; 
	}
	double lambda = gpp * A / (2.*Ne*gp); 
	//~ 
	//~ double B = min( y, A / gp);
	//~ double lambda = gpp * B / (2.*Ne); 
	
	vec wk = zeros<vec>(mk); 
	double sumwk = 0.;
	for (k = 1; k <= mk; k++){
		wk(k-1) = max(0., g_x(k-1) * (gpp*A - k*(k-1)*gp) );
		sumwk += wk(k-1); 
	}
	if ( sumwk <= 0.){
		wk = zeros<vec>(mk); 
		wk(0) = 1.; 
	} else{
		wk = wk / sumwk; 
	}
	
	List o; 
	o["lambda"] = lambda;
	o["g_x"] = wk; 
	return o; 
}

//[[Rcpp::export()]]
List update_lambda_gx_node2(double A, double Ne, vec g_x, double y)
{ // see notes august 26
	int mk = g_x.size();
	double gpp = 0.;
	double gp = 0.; 
	int i,j,k;
	vec qk = zeros<vec>(mk);
	for (k=1; k<=mk; k++){
		gpp += g_x(k-1) * k * (k-1); 
		gp += g_x(k-1) * k; 
		qk(k-1) = k * (k-1) * g_x(k-1); 
	}
	qk /= gpp; 
	double lambda = gpp * A / (2.*Ne*gp); 
	//~ 
	//~ double B = min( y, A / gp);
	//~ double lambda = gpp * B / (2.*Ne); 
	
	vec wk = zeros<vec>(mk); 
	double sumwk = 0.;
	
	double B = A / gp; 
	
	vec nk = B * g_x; 
	vec nk2 = zeros<vec>(mk); // after change
	for (k = 1; k <= mk; k++){
		if (k< mk){
			nk2(k-1) = max(0., (nk(k-1)-1.)*qk(k-1) + (nk(k-1)+1.)*qk(k) + nk(k-1)*(1-qk(k-1)-qk(k))); 
		} else{
			nk2(k-1) = max(0., (nk(k-1)-1.)*qk(k-1) + nk(k-1)*(1-qk(k-1))); 
		}
	}
	wk = nk2 / sum(nk2); 
	
	List o; 
	o["lambda"] = lambda;
	o["g_x"] = wk; 
	return o; 
}


//[[Rcpp::export()]]
List eventTimes2extant( vec eventTimes, vec nodeheights, vec parentheights )
{
	List o(eventTimes.size());
	List nodesAtHeight( eventTimes.size());  
	double h, nodeheight, parentheight ; 
	for (int i =0  ; i < eventTimes.size(); i++){
		h = eventTimes(i); 
		std::vector<int> x; 
		x.reserve(nodeheights.size()); 
		std::vector<int> nah; 
		nah.reserve( nodeheights.size()); 
		for (int u = 0 ; u < nodeheights.size(); u++){
			nodeheight = nodeheights(u);
			parentheight = parentheights(u); 
			if (nodeheight==h){
				nah.push_back( u + 1 );
			}
			if (!NumericVector::is_na( parentheight)){
				if ( nodeheight<= h && parentheight > h ){
					x.push_back( u+1 ); 
				}	
			} else{
				if (nodeheight <= h ){
					x.push_back( u+1 ) ; 
				}
			}
		}
		o[i] = wrap(x); 
		nodesAtHeight[i] = wrap(nah); 
	} 
	List oo; 
	oo["extant"] = o;
	oo["nodesAtHeight"] = nodesAtHeight ; 
	return oo;
}
