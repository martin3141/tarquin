#include <stdlib.h>
#include "../common/CBasis.hpp"
#include "../common/common.hpp"

#include <stdlib.h>
#include <algorithm>
#include <complex>
#include <sstream>
#include <iomanip>

using namespace tarquin;
using namespace std;

int main()
{
    /*
    for ( double echo = 30; echo < 400; echo = echo + 10 )
    {
        // create the dummy FID object
        CFID fid;

        fid.SetSamplingFrequency(2000.0);
        fid.SetTransmitterFrequency(2*6.3866e7);
        fid.SetEchoTime(echo*1e-3);
        fid.SetPPMRef(4.65);
        fid.SetNumberOfPoints(2048);
        //fid.SetNumberOfPoints(2*64);

        // create the basis object
        CBasis basis;

        // simulate a basis (internal) to match this FID
        tarquin::CBoswell log(tarquin::LOG_NOWHERE);
        basis.Simulate(fid, log);

        // define matrix of signals (no groups) S 
        const cvm::cmatrix& Q = basis.GetBasisMatrix();
        cvm::cmatrix S = Q;

        // a represents signal amplitudes, generate from gauss dist. rand numbers
        cvm::cvector a(S.nsize());
        cvm::cvector y(S.msize());

        
        // one to use
        //double ahat_array[] = {0.168526, 3.39376, 4.78156, 0, 2.57716, 0.251956, 1.54555, 4.59123, 0.855986, 0.923358, 3.56007, 0.0591376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.09216, 0.596475, 0.379632, 0.297313, 1};
        
        // no m-ins
        //double ahat_array[] = {0.168526, 3.39376, 4.78156, 0, 2.57716, 0.251956, 1.54555, 4.59123, 0.855986, 0.923358, 0, 0.0591376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.09216, 0.596475, 0.379632, 0.297313, 1};
        
        // no glu 
        //double ahat_array[] = {0.168526, 3.39376, 4.78156, 0, 2.57716, 0.251956, 1.54555, 0, 0.855986, 0.923358, 3.56007, 0.0591376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.09216, 0.596475, 0.379632, 0.297313, 1};
        
        // no gln
        //double ahat_array[] = {0.168526, 3.39376, 4.78156, 0, 2.57716, 0.251956, 0, 4.59123, 0.855986, 0.923358, 3.56007, 0.0591376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.09216, 0.596475, 0.379632, 0.297313, 1};
        
        // no NAA
        //double ahat_array[] = {0.168526, 3.39376, 4.78156, 0, 2.57716, 0.251956, 1.54555, 4.59123, 0.855986, 0.923358, 3.56007, 0.0591376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.596475, 0.379632, 0.297313, 1};
        
        // no GABA
        //double ahat_array[] = {0.168526, 3.39376, 4.78156, 0, 0, 0.251956, 1.54555, 4.59123, 0.855986, 0.923358, 3.56007, 0.0591376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.09216, 0.596475, 0.379632, 0.297313, 1};
        
        // no Asp
        double ahat_array[] = {0.168526, 0, 4.78156, 0, 2.57716, 0.251956, 1.54555, 4.59123, 0.855986, 0.923358, 3.56007, 0.0591376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.09216, 0.596475, 0.379632, 0.297313, 1};
        
        //double ahat_array[] = {0.168526, 3.39376, 0, 0, 2.57716, 0.251956, 1.54555, 4.59123, 0, 0.923358, 3.56007, 0.0591376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.297313, 1};

        //double ahat_array[] = {0.168526, 3.39376, 4.78156, 0, 2.57716, 0.251956, 1.54555, 4.59123, 0.855986, 0.923358, 3.56007, 0.0591376, 1.69248, 1.26724, 0, 0, 1.64341, 0.329271, 4.56396, 0.84077, 8.49015, 6.09216, 0.596475, 0.379632, 0.297313, 0};

        for ( int n = 0; n < S.nsize(); n++ )
        {
            a(n+1) = ahat_array[n];
            cout << basis.GetSignalName(n) << "\t" << a(n+1).real() << endl;
        }

        y = S*a;

        // now save y as .dpt format

        ostringstream strFilename;
        strFilename << "./sims/simulated_";
        strFilename << setw(3) << setfill('0');
        strFilename << echo << ".dpt";

        CFID dptout = fid;
        dptout.AppendFromVector(y);
        dptout.SaveToFile(strFilename.str());
    }
    
    */

    return 0;
}
