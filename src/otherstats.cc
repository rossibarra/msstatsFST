#include <otherstats.hpp>
#include <cmath>

using namespace std;
using namespace Sequence;

std::pair<double,double> RozasR(const SimData & matrix, const double & thetapi, const unsigned & segsites ) 
/*
  Ramon-Osnin & Rozas' (2005, MBE 19: 2092) R2 and R2E statistics, returned in that order in a pair.
  Code by Murray Cox w/some modification by KRT
*/
{
  if(segsites == 0) return std::make_pair(strtod("NAN",NULL),strtod("NAN",NULL));

  unsigned int nsam = matrix.size();
  unsigned int sites = matrix.numsites();

  double sum_f=0.,sum_u=0.;
  unsigned folded_count,unfolded_count;

  for ( unsigned ind = 0; ind < nsam; ++ind ) 
    { 
      folded_count = unfolded_count = 0;
	    
      for ( unsigned site = 0; site < sites; ++site ) 
	{
	  unsigned zero_count = 0, one_count = 0;
	  switch ( matrix[ind][site] ) 
	    {
	    case '0':
	      for ( unsigned int loop0 = 0; loop0 < nsam; ++loop0 ) 
		{                      
		  if ( matrix[loop0][site] == '0' ) ++zero_count;
		}
	      if ( zero_count == 1 ) 
		{
		  ++folded_count;
		}
	      break;
	    case '1':
	      for ( unsigned int loop1 = 0; loop1 < nsam; ++loop1 ) 
		{            
		  if ( matrix[loop1][site] == '1' ) ++one_count;
		}
	      if ( one_count == 1 ) 
		{
		  ++folded_count;
		  ++unfolded_count;
		}
	      break;
	    }
	}
      sum_f += pow(( folded_count - ( thetapi / 2. )), 2); 
      sum_u += pow(( unfolded_count - ( thetapi / 2. )), 2); 
    }
  return make_pair( sqrt(sum_f/double(nsam))/double(segsites), sqrt(sum_u/double(nsam))/double(segsites) );
}

unsigned Rm_MG( const Sequence::SimData & matrix, unsigned segsites, unsigned nhaps )
/*
  Returns a simple bound on the minumum number of recombination events, calculated according
  to equation 4 of Myers & Griffiths paper
*/
{
  bool anc = ( std::find(matrix.begin(),matrix.end(),std::string(segsites,'0')) != matrix.end() );
  return (nhaps > segsites) ? ( (anc) ? nhaps-segsites-1 : nhaps-segsites ) : 0;
}
