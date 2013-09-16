#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <stdio.h>
#include <fstream>
#include <boost/filesystem.hpp>

#include "CNNLSSolverLH.hpp"
#include "NNLSSolverBJ.hpp"
#include "FNNLSSolverBJ.hpp"

namespace fs = boost::filesystem;

BOOST_AUTO_TEST_SUITE( nnls )

BOOST_AUTO_TEST_CASE( nnls_lh )
{
    cvm::rmatrix Z(4, 3);
    cvm::rvector x(4);
    cvm::rvector d(3);

    Z[1][1] = 59.0f;
    Z[1][2] = 51.0f;
    Z[1][3] = 55.0f;

    Z[2][1] = 22.0f;
    Z[2][2] = 70.0f;
    Z[2][3] = 14.0f;

    Z[3][1] = 75.0f;
    Z[3][2] = 89.0f;
    Z[3][3] = 15.0f;

    Z[4][1] = 26.0f;
    Z[4][2] = 96.0f;
    Z[4][3] = 26.0f;

    x[1] = 84.0f;
    x[2] = 25.0f;
    x[3] = 81.0f;
    x[4] = 24.0f;

    d[1] = 1.0148;
    d[2] = 0.00;
    d[3] = 0.3337;


    cvm::rvector dhat(3);

    tarquin::CNNLSSolverLH nnls;
    nnls.solve(Z, x, dhat, false);

    for( int i = 0; i < dhat.size(); ++i )
    {
        tarquin::treal diff = std::fabs(dhat[1+i] - d[1+i]);
        BOOST_CHECK_LT(diff, 1e-3);
    }
}

BOOST_AUTO_TEST_CASE( nnls_bj )
{
    tarquin::NNLSSolverBJ nnls;

    cvm::rmatrix Z(4, 3);
    cvm::rvector x(4);
    cvm::rvector d(3);

    Z[1][1] = 59.0f;
    Z[1][2] = 51.0f;
    Z[1][3] = 55.0f;

    Z[2][1] = 22.0f;
    Z[2][2] = 70.0f;
    Z[2][3] = 14.0f;

    Z[3][1] = 75.0f;
    Z[3][2] = 89.0f;
    Z[3][3] = 15.0f;

    Z[4][1] = 26.0f;
    Z[4][2] = 96.0f;
    Z[4][3] = 26.0f;

    x[1] = 84.0f;
    x[2] = 25.0f;
    x[3] = 81.0f;
    x[4] = 24.0f;

    d[1] = 1.0148;
    d[2] = 0.00;
    d[3] = 0.3337;

    cvm::rvector dhat(3);

    nnls.solve(Z, x, dhat, false);

    for( int i = 0; i < dhat.size(); ++i )
    {
        tarquin::treal diff = std::fabs(dhat[1+i] - d[1+i]);
        BOOST_CHECK_LT(diff, 1e-3);
    }
}

BOOST_AUTO_TEST_CASE( fnnls_bj )
{
    tarquin::FNNLSSolverBJ nnls;

    cvm::rmatrix Z(4, 3);
    cvm::rvector x(4);
    cvm::rvector d(3);

    Z[1][1] = 73.0;
    Z[1][2] = 71.0;
    Z[1][3] = 52.0;

    Z[2][1] = 87.0;
    Z[2][2] = 74.0;
    Z[2][3] = 46.0;

    Z[3][1] = 72.0;
    Z[3][2] = 2.0;
    Z[3][3] = 7.0;

    Z[4][1] = 80.0;
    Z[4][2] = 89.0;
    Z[4][3] = 71.0;

    x[1] = 49.0;
    x[2] = 67.0;
    x[3] = 68.0;
    x[4] = 20.0;

    d[1] = 0.65;
    d[2] = 0.00;
    d[3] = 0.00;

    cvm::rvector dhat(3);

    nnls.solve(Z, x, dhat, false);

    for( int i = 0; i < dhat.size(); ++i )
    {
        tarquin::treal diff = std::fabs(dhat[1+i] - d[1+i]);
        BOOST_CHECK_LT(diff, 1e-3);
    }
}

BOOST_AUTO_TEST_SUITE_END()
