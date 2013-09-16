#include "NNLSSolverBJ.hpp"
#include <set>
#include <cmath>
#include <stdexcept>


void tarquin::NNLSSolverBJ::solve(cvm::rmatrix& Z, cvm::rvector& x, cvm::rvector& d, bool bUseStartingValue)
{
    const integer L = Z.msize();
    const integer M = Z.nsize();

    std::set<integer> P;
    std::set<integer> R;

    cvm::rmatrix Zp(Z.msize(), Z.nsize());

    treal tol = NUMERICAL_TOL_NNLS;

    for(integer i = 1; i <= M; i++) 
    {
        R.insert(i);
        d[i] = 0.0;
    }

    cvm::rvector w = ~Z * (x - Z*d);

    cvm::rvector s(d.size());

    while( R.size() > 0 )
    {
        // m = arg max_{n\in R} w
        treal wmax = -std::numeric_limits<treal>::max();
        integer m = 0;
        for( std::set<integer>::iterator itj = R.begin(); itj != R.end(); itj++ ) 
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

        for( ; ;  )
        {
            cvm::rmatrix ZP(L, P.size());
            {
                int i = 1;
                for( std::set<integer>::iterator itp = P.begin(); itp != P.end(); itp++, ++i )	
                    ZP(i) = Z(*itp);
            }

            cvm::rmatrix G  = (~ZP)*ZP;
            cvm::rvector sp = G.pinv(tol) * (~ZP *x);
            
            // update s
            int i = 1;
            for( std::set<integer>::iterator itp = P.begin(); itp != P.end(); itp++, ++i )	
                s[*itp] = sp[i];

            for( std::set<integer>::iterator itr = R.begin(); itr != R.end(); itr++ )	
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

            treal min_z = std::numeric_limits<treal>::max();
            for( std::set<integer>::iterator itp = P.begin(); itp != P.end(); itp++ )	
            {
                if( s[*itp] <= tol )
                {
                    float z = d[*itp] / (d[*itp] - s[*itp]);
                    if( z < min_z )
                        min_z = z;
                }
            }

            treal alpha = min_z;

            d = d + alpha*(s - d);


            {
                std::set<integer> Pdel;
                for( std::set<integer>::iterator itp = P.begin(); itp != P.end(); itp++ )	
                {
                    if( std::abs(d(*itp)) < tol ) 
                    {
                        R.insert(*itp);
                        Pdel.insert(*itp);
                    }
                }

                for( std::set<integer>::iterator itp = Pdel.begin(); itp != Pdel.end(); itp++ )	
                    P.erase(*itp);
            }
        }

        d = s;
        w = ~Z * (x - Z*d);
    }
}

