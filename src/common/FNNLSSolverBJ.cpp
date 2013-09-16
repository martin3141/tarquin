#include "FNNLSSolverBJ.hpp"
#include <set>
#include <cmath>


void tarquin::FNNLSSolverBJ::solve(cvm::rmatrix& Z, cvm::rvector& x, cvm::rvector& d, bool use_starting_value)
{
    cvm::rmatrix ZTZ = ~Z * Z;
    cvm::rvector ZTx = ~Z * x;

    const integer L = Z.msize();
    const integer M = Z.nsize();

    std::set<integer> P;
    std::set<integer> R;

    treal tol = NUMERICAL_TOL_NNLS;

    bool outer_set_update = true;

    if( use_starting_value )
    {
        for(integer i = 1; i <= M; ++i) 
        {
            if( d[i] > tol )
                P.insert(i);
            else
                R.insert(i);
        }

        outer_set_update = false;
    }
    else
    {
        for(integer i = 1; i <= M; ++i) 
        {
            R.insert(i);
            d[i] = 0.0;
        }
    }

    cvm::rvector w = ZTx - ZTZ*d;

    cvm::rvector s(d.size());

    bool first_pass = true;

    while( R.size() > 0 || (use_starting_value && first_pass) )
    {
        if( outer_set_update )
        {
            // m = arg max_{n\in R} w
            treal wmax = -std::numeric_limits<treal>::max();
            integer m = 0;
            for( std::set<integer>::iterator itj = R.begin(); itj != R.end(); ++itj ) 
            {
                integer j = *itj;

                if( w(j) >= wmax ) 
                {
                    m = j;
                    wmax = w(j);
                }
            }

            if( wmax < tol )
                return;

            P.insert(m);
            R.erase(m);
        }

        for( int iter = 0; iter < 100*L; ++iter )
        {
            cvm::rmatrix ZTZP(P.size(), P.size());
            cvm::rvector ZTxP(P.size());
            {
               
               int i = 1;
               for( std::set<integer>::iterator itp = P.begin(); itp != P.end(); ++itp, ++i )	
                   ZTxP(i) = ZTx(*itp);

               int ii = 1;
               int jj = 1;
               for( std::set<integer>::iterator itpi = P.begin(); itpi != P.end(); ++itpi, ++ii )	
               {
                   jj = 1;
                   for( std::set<integer>::iterator itpj = P.begin(); itpj != P.end(); ++itpj, ++jj )	
                       ZTZP[ii][jj] = ZTZ[*itpi][*itpj];
               }
            }

            if ( ZTZP.msize() == 0 ) // added by MW 11 Sept 13 to prevent a crash
                break;
            
            cvm::rvector sp = ZTZP.pinv(tol) * ZTxP;

            int i = 1;
            for( std::set<integer>::iterator itp = P.begin(); itp != P.end(); ++itp, ++i )	
                s[*itp] = sp[i];

            for( std::set<integer>::iterator itr = R.begin(); itr != R.end(); ++itr )	
                s[*itr] = 0.0;

            treal min_sp = std::numeric_limits<treal>::max();
            for( int i = 1; i <= sp.size(); ++i )
            {
                if( sp[i] < min_sp )
                    min_sp = sp[i];
            }

            // no need to enter inner loop
            if( min_sp > 0.0 )
                break;

            treal alpha = std::numeric_limits<treal>::max();
            for( std::set<integer>::iterator itp = P.begin(); itp != P.end(); ++itp )	
            {
                if( s[*itp] <= tol )
                {
                    float z = d[*itp] / (d[*itp] - s[*itp]);
                    if( z < alpha )
                        alpha = z;
                }
            }

            d += alpha*(s - d);

            {
                std::set<integer> Pdel;
                for( std::set<integer>::iterator itp = P.begin(); itp != P.end(); ++itp )	
                {
                    if( std::abs(d(*itp)) < tol ) 
                    {
                        R.insert(*itp);
                        Pdel.insert(*itp);
                    }
                }

                for( std::set<integer>::iterator itp = Pdel.begin(); itp != Pdel.end(); ++itp )	
                    P.erase(*itp);
            }
        }

        d = s;
        w = ZTx - ZTZ*d;

        first_pass       = false;
        outer_set_update = true;
    }
}

