/* 
msstats - read data from ms via stdin, calculate common summary statistics

Copyright (C) 2002-2007 Kevin Thornton

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include <iostream>
#include <vector>
#include <cstdio>
#include <numeric>
#include <utility>
#include <Sequence/SimParams.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <Sequence/SeqConstants.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/FST.hpp>
#include <otherstats.hpp>

using namespace std;
using namespace Sequence;

void calcstats(const SimData & d, const unsigned & mincount);

int main(int argc, char *argv[]) {
	std::vector<unsigned int> config;
	bool multipop = false;
	int mincount = 1;
	for(int arg = 1 ; arg < argc ; ++arg){
		if( string(argv[arg]) == "-I" ){
			multipop = true;
			int npop = atoi(argv[++arg]);
			for( int i=0;i<npop;++i ){
				config.push_back(atoi(argv[++arg]));
			}
		}
		else if (string(argv[arg]) == "-m"){
			mincount = atoi(argv[++arg]);
		}
	}
	int total = std::accumulate(config.begin(),config.end(),0,plus<int>());
	SimParams p;
	p.fromfile(stdin);
	SimData d;
	
	#if __GNUG__ && __GNUC__ >= 3
	std::ios_base::sync_with_stdio(true);
	#endif
	if(multipop){
		std::cout << "rep\tpop\t";
	}
	std::cout << "S\t"
	<< "n1\t"
	<< "next\t"
	<< "theta\t"
	<< "pi\t"
	<< "thetaH\t"
	<< "tajd\t"
	<< "fulif\t"
	<< "fulid\t"
	<< "fulifs\t"
	<< "fulids\t"
	<< "rm\t"
	<< "rmmg\t"
	<< "nhaps\t"
	<< "hdiv\t"
	<< "wallb\t"
	<< "wallq\t"
	<< "rosasrf\t"
	<< "rosasru\t"
	<< "zns\t"
	<< "FST"
	<< endl;
	int rv;
	int rep=0;
	
	while( (rv=d.fromfile(stdin)) != EOF ){

		if(!multipop){
			calcstats(d,mincount);
		}
		else{
			SimData d2;
			if(d.size() != total){
				std::cerr << "oh crap\n";
				exit(10);
			}
			int sum = 0;
			
/************START hack code for Weir & Cockerham Fst*/
			vector < vector < double > > ps (config.size());  
			vector < double > pbar (d.numsites(),0), 
				hbar (d.numsites(),0), 
				ssq (d.numsites(),0), 
				f (d.numsites(),0),
				a (d.numsites(),0), 
				b (d.numsites(),0), 
				c (d.numsites(),0);
			double sumnsq=0;
			double r = config.size(); 
		
			for(int i = 0 ; i < config.size() ; ++i){
				d2.assign(&*d.pbegin(),d.numsites(),
				&d[sum],config[i]);
				sum += config[i];						
				sumnsq+=config[i]*config[i];
				for (unsigned k = 0 ; k < d2.numsites() ; ++k){    // sites        
                 stateCounter Counts; 
                 for (unsigned j = 0 ; j < d2.size() ; ++j){ // inds
                 		Counts(d2[j][k]);              		
                 }
                  int crap = Counts.one; 
                	ps[i].push_back(double(Counts.one)/double(d2.size()));
           	}
			}
			
			double nbar=sum/config.size();
			for(int i = 0 ; i < config.size() ; ++i){
				for (unsigned k = 0 ; k < d.numsites() ; ++k){    // sites
					pbar[k]+=ps[i][k]*config[i]/(config.size()*nbar);
					double hi=(1-ps[i][k]*ps[i][k]-(1-ps[i][k])*(1-ps[i][k]));
					hbar[k]+=config[i]*hi/(r*nbar);
				}
			}
			for(int i = 0 ; i < config.size() ; ++i){
				for (unsigned k = 0 ; k < d2.numsites() ; ++k){    // sites
					ssq[k]+=config[i]*(ps[i][k]-pbar[k])*(ps[i][k]-pbar[k])/((r-1)*nbar);
				}
			}
			double nc=(r*nbar-sumnsq/(r*nbar))/(r-1);
			double fbar=0;
			for (unsigned k = 0 ; k < d.numsites() ; ++k){ 
				c[k] = 0.5*hbar[k];
				b[k]=(nbar/(nbar-1))*( pbar[k]*(1-pbar[k])-((r-1)/r)*ssq[k]-(2*nbar-1)*hbar[k]/(4*nbar) );
				a[k]=(nbar/nc)*ssq[k]-(1/(nbar-1))*( pbar[k]*(1-pbar[k])-((r-1)/r)*ssq[k]-hbar[k]/4 );
				f[k] = a[k]/(a[k]+b[k]+c[k]);
				fbar+=f[k]/d.numsites();
			}
/************END hack code for Weir & Cockerham Fst*/
			sum = 0;	// JRI mod
			for(int i = 0 ; i < config.size() ; ++i){
				d2.assign(&*d.pbegin(),d.numsites(),
				&d[sum],config[i]);
				sum += config[i];
				cout << rep << '\t' << i << '\t';
				RemoveInvariantColumns(&d2);
				calcstats(d2,mincount);
				if( i < config.size()-1 ){ cout << "nan\n"; }   // JRI mod
			}
			cout << fbar << endl; // JRI mod
		}
		++rep;
	} 
}

void calcstats(const SimData & d, const unsigned & mincount){
	PolySIM P(&d);

	cout << P.NumPoly()   << '\t' 
		<< P.NumSingletons() << '\t'
		<< P.NumExternalMutations() << '\t'
		<< P.ThetaW()    << '\t' 
		<< P.ThetaPi()   << '\t'
		<< P.ThetaH()    << '\t' 
		<< P.TajimasD()  << '\t'
		<< P.FuLiD()     << '\t'
		<< P.FuLiF()     << '\t'
		<< P.FuLiFStar() << '\t'
		<< P.FuLiDStar() << '\t';
	unsigned rm = P.Minrec(),nhaps=P.DandVK();
	cout << ((rm != SEQMAXUNSIGNED) ? rm : strtod("NAN",NULL)) << '\t'
		<< Rm_MG(d,P.NumPoly(),nhaps) << '\t'
		<< nhaps    << '\t'
		<< P.DandVH()    << '\t'
		<< P.WallsB()    << '\t'
		<< P.WallsQ()    << '\t';
	pair<double,double> r2 = RozasR(d,P.ThetaPi(),P.NumPoly());
	cout << r2.first << '\t' << r2.second << '\t';
	if(d.numsites() > 1){
		vector< vector<double> > ld = P.Disequilibrium(mincount);
		if(ld.empty()){  //no sites made the mincount cutoff
			cout << "nan\t";
		}
		else{
			double zns = 0.;
			for(unsigned i=0;i<ld.size();++i){
				zns += ld[i][2];
			}
			zns /= double(ld.size());
			cout << zns << '\t';
		}
	}
	else {
		cout << "nan\t";
	}
}
