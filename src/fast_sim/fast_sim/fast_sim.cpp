#include <iostream>
#include "fast_sim.hpp"
#include <iostream>
#include <math.h>
#include <complex>
#include <fstream>

void get_spectrum(dcv& time_sig, drv& freq_sig, int type)
{
    int N = time_sig.size();
    Eigen::FFT<double> fft;
    dcv freq_sig_temp;
    fft.fwd(freq_sig_temp, time_sig);
    dcv freq_sig_shift(N);
    freq_sig_shift << freq_sig_temp.tail(N/2),freq_sig_temp.head(N/2);
    if ( type == 1 )
        freq_sig = freq_sig_shift.real();
    else if (type == 2 )
        freq_sig = freq_sig_shift.imag();
    else if (type == 3 )
        freq_sig = freq_sig_shift.array().abs();
    else
        freq_sig = freq_sig_shift.real();
}

void read_points(std::string fname, drv& x, drv& y)
{
    std::ifstream infile(fname.c_str());
    
    if (!infile.is_open())
        std::cout << "Error opening file." << std::endl;

    int N = 0;
    std::string line;
    while (std::getline(infile, line))
        ++N;

    x.resize(N);
    x = drv::Zero(N);
    y.resize(N);
    y = drv::Zero(N);

    infile.clear();
    infile.seekg(0,std::ios::beg);
   
    int n = 0;
    double a, b;
    while (infile >> a >> b)
    {
        x(n) = a;
        y(n) = b;
        n++;
    }
}

void write_array(std::string fname, drv& x)
{
    std::ofstream fout(fname.c_str(), std::ios::out);
    for ( int n = 0; n < x.size(); n++ )
        fout << x(n) << std::endl;
}

void write_points(std::string fname, drv& x, drv& y)
{
    if ( x.size() != y.size() )
    {
        std::cout << "Error : plot size mistmatch." << std::endl;
        return;
    }

    std::ofstream fout(fname.c_str(), std::ios::out);
    for ( int n = 0; n < x.size(); n++ )
        fout << x(n) << " " << y(n) << std::endl;
}

void plot(drv& x, drv& y, double xmin, double xmax, std::string title)
{
    if ( x.size() != y.size() )
    {
        std::cout << "Error : plot size mistmatch." << std::endl;
        return;
    }

    std::ofstream fout("file_out.csv", std::ios::out);
    for ( int n = 0; n < x.size(); n++ )
        fout << x(n) << ", " << y(n) << std::endl;

    std::ofstream gnuplot_outfile("spec.gnu");
	gnuplot_outfile << "set xrange [" << xmin << ":" << xmax << "]" << std::endl;
	gnuplot_outfile << "set nokey" << std::endl;
	gnuplot_outfile << "set title \"" << title << "\"" << std::endl;
    gnuplot_outfile	<< "plot 'file_out.csv' using 1:2 w l" << std::endl;
    gnuplot_outfile << std::endl;
    gnuplot_outfile << "pause -1";
    gnuplot_outfile << std::endl;
    std::string strRunMe;
    strRunMe = "gnuplot spec.gnu";
    int out = system(strRunMe.c_str());
}

void plot(drv& x, drv& y, std::string title)
{
    if ( x.size() != y.size() )
    {
        std::cout << "Error : plot size mistmatch." << std::endl;
        return;
    }

    std::ofstream fout("file_out.csv", std::ios::out);
    for ( int n = 0; n < x.size(); n++ )
        fout << x(n) << ", " << y(n) << std::endl;

    std::ofstream gnuplot_outfile("spec.gnu");
	gnuplot_outfile << "set nokey" << std::endl;
	gnuplot_outfile << "set title \"" << title << "\"" << std::endl;
    gnuplot_outfile	<< "plot 'file_out.csv' using 1:2 w l" << std::endl;
    gnuplot_outfile << std::endl;
    gnuplot_outfile << "pause -1";
    gnuplot_outfile << std::endl;
    std::string strRunMe;
    strRunMe = "gnuplot spec.gnu";
    int out = system(strRunMe.c_str());
}

void plot(drv& x, std::string title)
{
    std::ofstream fout("file_out.csv", std::ios::out);
    for ( int n = 0; n < x.size(); n++ )
        fout << x(n) << std::endl;

    std::ofstream gnuplot_outfile("spec.gnu");
	gnuplot_outfile << "set nokey" << std::endl;
	gnuplot_outfile << "set title \"" << title << "\"" << std::endl;
    gnuplot_outfile	<< "plot 'file_out.csv' w l" << std::endl;
    gnuplot_outfile << std::endl;
    gnuplot_outfile << "pause -1";
    gnuplot_outfile << std::endl;
    std::string strRunMe;
    strRunMe = "gnuplot spec.gnu";
    int out = system(strRunMe.c_str());
}


void plot_pdf(std::string fname, drv& x, drv& y, double xmin, double xmax, std::string title)
{
    if ( x.size() != y.size() )
    {
        std::cout << "Error : plot size mistmatch." << std::endl;
        return;
    }

    std::ofstream fout("file_out.csv", std::ios::out);
    for ( int n = 0; n < x.size(); n++ )
        fout << x(n) << ", " << y(n) << std::endl;

    std::ofstream gnuplot_outfile("spec.gnu");
	gnuplot_outfile << "set term pdf" << std::endl;
	gnuplot_outfile << "set output '" << fname << "'" << std::endl;
	gnuplot_outfile << "set xrange [" << xmin << ":" << xmax << "]" << std::endl;
	gnuplot_outfile << "set nokey" << std::endl;
	gnuplot_outfile << "set title \"" << title << "\"" << std::endl;
    gnuplot_outfile	<< "plot 'file_out.csv' w l" << std::endl;
    gnuplot_outfile << std::endl;
    std::string strRunMe;
    strRunMe = "gnuplot spec.gnu";
    int out = system(strRunMe.c_str());
}

void plot_pdf(std::string fname, drv& x, std::string title)
{
    std::ofstream fout("file_out.csv", std::ios::out);
    for ( int n = 0; n < x.size(); n++ )
        fout << x(n) << std::endl;

    std::ofstream gnuplot_outfile("spec.gnu");
	gnuplot_outfile << "set term pdf" << std::endl;
	gnuplot_outfile << "set output '" << fname << "'" << std::endl;
	gnuplot_outfile << "set nokey" << std::endl;
	gnuplot_outfile << "set title \"" << title << "\"" << std::endl;
    gnuplot_outfile	<< "plot 'file_out.csv' using 1 w l" << std::endl;
    gnuplot_outfile << std::endl;
    std::string strRunMe;
    strRunMe = "gnuplot spec.gnu";
    int out = system(strRunMe.c_str());
}


double spin_sys::get_tip_angle(drv& profile, drv& freq_scale, double angle, double frequency)
{
    // TODO linear interpolation

    drv ppm_scale = -freq_scale/(m_B0*1e-6);
    /*
    std::cout << "start:" << ppm_scale(0) << std::endl;
    std::cout << "end  :" << ppm_scale(ppm_scale.size()-1) << std::endl;
    std::cout << "freq :" << frequency << std::endl;
    */

    if ( frequency < ppm_scale(ppm_scale.size()-1) )
    {
        return 0;
    }
    else if ( frequency > ppm_scale(0) )
    {
        return 0;
    }
    else
    {
        //double max_profile = profile.maxCoeff();
        double max_profile = 1; // don't rescale
        int num;
        double profile_pt = (ppm_scale.array()-frequency).abs().minCoeff(&num);
        double ret_angle;
       
        if ( frequency == profile_pt )
        {
            ret_angle = profile(num)/max_profile*angle;
        }
        else
        {
            // linearly interpolate between points
            double x_left_pt, y_left_pt, x_right_pt, y_right_pt;
            if ( frequency > ppm_scale(num) )
            {
                x_left_pt = ppm_scale(num);
                x_right_pt = ppm_scale(num-1);
                y_left_pt = profile(num);
                y_right_pt = profile(num-1);
            }
            else
            {
                x_left_pt = ppm_scale(num-1);
                x_right_pt = ppm_scale(num);
                y_left_pt = profile(num-1);
                y_right_pt = profile(num);
            }
            
            //double y_interp = profile(num); // no intep
            
            double y_interp = y_left_pt + (frequency - x_left_pt) * (y_right_pt - y_left_pt) / (x_right_pt - x_left_pt);
                
            //std::cout << x_left_pt << " " << frequency << " " << x_right_pt << std::endl;
            //std::cout << y_left_pt << " " << y_interp << " " << y_right_pt << std::endl;

            ret_angle = y_interp/max_profile*angle;
        }

        return ret_angle;
    }
}

void gen_sinc(drv& waveform, drv& phase, int n, int pts, int zp_pts)
{
    if ( zp_pts < pts )
        zp_pts = pts;

    waveform.resize(zp_pts);
    waveform = drv::Zero(zp_pts);
    phase.resize(zp_pts);
    phase = drv::Zero(zp_pts);
    
    /*
    for ( int i = 0; i < pts; i++ )
    {
        double t = (i - (pts - 1) / 2.0) / (pts) * 2;
        if ( t == 0 )
            waveform(i) = 1;
        else
            waveform(i) = sin( n * 2.0 * M_PI * t) / ( n * 2.0 * M_PI * t);
    }*/

    // same as gamma
    /*
    for ( int i = 0; i < pts; i++ )
    {
        double t = (i - (pts - 1) / 2.0) / (pts + 1) * 2;
        if ( t == 0 )
            waveform(i) = 1;
        else
            waveform(i) = sin( n * M_PI * t) / ( n * M_PI * t);
    }
    */

    // less truncated version than gamma produces
    for ( int i = 0; i < pts; i++ )
    {
        double t = (i - (pts - 1) / 2.0) / (pts + 0) * 2;
        if ( t == 0 )
            waveform(i) = 1;
        else
            waveform(i) = sin( n * M_PI * t) / ( n * M_PI * t);
    }
}

void gen_gaus(drv& waveform, drv& phase, double trunc_fact, int pts, int zp_pts)
{
    if ( zp_pts < pts )
        zp_pts = pts;

    waveform.resize(zp_pts);
    phase.resize(zp_pts);
    phase = drv::Zero(zp_pts);

    double beta = log(trunc_fact);

    for ( int i = 0; i < pts; i++ )
    {
        //double t = (i - (pts - 1) / 2.0) / (pts + 1) * 2;

        double t = (i - (pts-1) / 2.0) / (pts-1) * 2;
        //std::cout << t << std::endl;
        waveform(i) = exp( beta * t * t );
    }
}

void gen_ham(drv& waveform, int pts, int zp_pts)
{
    if ( zp_pts < pts )
        zp_pts = pts;

    waveform.resize(zp_pts);
    for ( int i = 0; i < pts; i++ )
        waveform(i) = 0.54 - 0.46 * cos( 2.0 * M_PI * i / ( pts - 1 ) );
}

void gen_sinc_ham(drv& waveform, drv& phase, int n, int pts, int zp_pts)
{
    drv sinc_waveform;
    drv sinc_phase;
    gen_sinc(sinc_waveform, sinc_phase, n , pts, zp_pts);
    
    drv ham_waveform;
    drv ham_phase;
    gen_ham(ham_waveform, pts, zp_pts);

    waveform.resize(sinc_waveform.size());
    waveform = sinc_waveform.array() * ham_waveform.array();
    phase = sinc_phase;
}

void d2s(dcm& in, scm& out)
{
	double tol = 1e-19;
	out.resize(in.rows(), in.cols());
	for ( int n = 0; n < in.rows(); n++ )
		for ( int m = 0; m < in.cols(); m++ )
			if ( std::abs(in(m, n)) > tol )
				out.coeffRef(n, m) = in(n, m);
}

// generate the pauli spin matrices
void Iz_pauli(double I, scm& prod_op)
{
	size_t dim = size_t( 2 * I + 1 + 0.5 );
	prod_op.resize(dim, dim);
	for ( size_t n = 0; n < dim; n++ )
		prod_op.coeffRef(n,n) = I - n;
}

void Ip_pauli(double I, scm& prod_op)
{
	// generate Ip
	size_t dim = size_t(2 * I + 1 + 0.5);
	prod_op.resize(dim, dim);
	for ( size_t n = 0; n < dim - 1;  n++ )
		prod_op.coeffRef(n,n+1) = pow(I * ( I + 1 ) - ( I - n - 1) * ( I - n ), 0.5);
}

void Im_pauli(double I, scm& prod_op)
{
	// generate Im
	size_t dim = size_t(2 * I + 1 + 0.5);
	prod_op.resize(dim, dim);
	for ( size_t n = 0; n < dim - 1;  n++ )
		prod_op.coeffRef(n+1,n) = pow(I * ( I + 1 ) - ( I - n - 1) * ( I - n), 0.5);
}

void Ix_pauli(double I, scm& prod_op)
{
	// generate Ix
	size_t dim = size_t(2 * I + 1 + 0.5);
	prod_op.resize(dim, dim);
	scm Ip;
	Ip_pauli(I, Ip);
	scm Im;
	Im_pauli(I, Im);
	prod_op = ( Ip + Im ) * 0.5;
}

void Iy_pauli(double I, scm& prod_op)
{
	// generate Iy
	size_t dim = size_t(2 * I + 1 + 0.5);
	prod_op.resize(dim, dim);
	scm Ip;
	Ip_pauli(I, Ip);
	scm Im;
	Im_pauli(I, Im);
	std::complex<double> j (0,1);
	prod_op = - 0.5 * j * ( Ip - Im );
}

spin_sys::spin_sys( const drv& spin_vec, const drv& chem_shift_vec, 
	const drm& j_coupling_mat, double B0, double offset)
{
	// TODO - sanity checks on chem_shift_vec and j_coupling mat
	m_nspins = spin_vec.size();

	m_spin_vec = spin_vec;
    // convert ppm to Hz
	//m_chem_shift_vec = -chem_shift_vec.array() * B0 * 1e-6;
	m_chem_shift_vec = chem_shift_vec;
	m_j_coupling_mat = j_coupling_mat;
	m_B0 = B0;
    m_offset = offset;

	// calculate the sizes of the operators in the
    // Zeeman product basis
	m_fsize = 1;
	for ( int n = 0; n < m_spin_vec.size(); n++ )
		m_fsize = m_fsize * size_t(2 * m_spin_vec(n) + 1 + 0.5);

    // calculate the full Hamiltonian
    H(m_H); 
    
    // perform a symmetric eigenvalue decomposition
    Eigen::SelfAdjointEigenSolver<dcm> Hsolver(m_H);
    
    // assign the member variables 
    m_H_eig_vals = Hsolver.eigenvalues();
    m_H_eig_vals_mat.resize(m_fsize, m_fsize);
    m_H_eig_vals_mat.diagonal() = m_H_eig_vals;
    m_H_eig_vecs = Hsolver.eigenvectors();

	// for the weak coupling case
	m_H_eig_vals_diag = m_H.diagonal().real();
    m_H_eig_vals_diag_mat.resize(m_fsize, m_fsize);
	m_H_eig_vals_diag_mat.setZero(m_fsize, m_fsize);
    m_H_eig_vals_diag_mat.diagonal() = m_H_eig_vals_diag;
	m_H_eig_vecs_diag = dcv::Constant(m_fsize, 1).asDiagonal();

	// apply weak coupling approx
	/*m_H_eig_vals = m_H_eig_vals_diag;
	m_H_eig_vals_mat = m_H_eig_vals_diag_mat;
	m_H_eig_vecs = m_H_eig_vecs_diag;
	*/
}

void spin_sys::gen_fft_profile(drv& waveform, int pul_N, double pul_t, drv& pul_freq, drv& pul_profile)
{
    double pulse_fs = pul_t / pul_N;
    double total_bw = 1/pulse_fs;
    int N = waveform.size();
    pul_freq.resize(N);
    
    for ( int f = 0; f < waveform.size(); f++ )
        pul_freq(f) = (-total_bw/2 + total_bw*f/waveform.size());

    //for ( int f = 0; f < waveform.size(); f++ )
    //    pul_freq(f) = -(-total_bw/2 + total_bw*f/waveform.size())/(m_B0*1e-6);

    pul_profile.resize(N);
    
    dcv freq_profile(N);
    Eigen::FFT<double> fft_prof;
    fft_prof.fwd(freq_profile, waveform);
    dcv freq_profile_shift(waveform.size());
    freq_profile_shift << freq_profile.tail(N/2),freq_profile.head(N/2);

    pul_profile = freq_profile_shift.array().abs();
}

void spin_sys::gen_freq_scale(drv& freq_scale, int N, double fs)
{
    double bw = fs;
    freq_scale.resize(N);
    
    for ( int n = 0; n < N; n++ )
        freq_scale(n) = -bw/2 + bw*n/N;
}

void spin_sys::gen_ppm_scale(drv& ppm_scale, int N, double fs, double offset)
{
    gen_freq_scale(ppm_scale, N, fs);
    ppm_scale /= -m_B0*1e-6;
    ppm_scale = ppm_scale.array() + offset; 
}

void spin_sys::crushgrad_shaped_ypulse(drv& waveform, drv& phase, double angle, double t, double offset)
{
    dcm rho_init = m_rho;

	dcm rho_combined;
    for ( int n = 0; n < 4; n++ )
    {
        if ( n > 0 )
            m_rho = rho_init;

        zpulse(n*90);

        shaped_ypulse(waveform, phase, angle, t, offset);

        zpulse(n*90);

		// add the total
		if ( n == 0 )
			rho_combined = m_rho;
		else
			rho_combined += m_rho;
    }

    rho_combined /= 4;
	m_rho = rho_combined;
	zero_multqcs();
}

void spin_sys::crushgrad_xpulse(double angle)
{
    dcm rho_init = m_rho;

	dcm rho_combined;
    for ( int n = 0; n < 4; n++ )
    {
        if ( n > 0 )
            m_rho = rho_init;

        zpulse(n*90);

        xpulse(angle);

        zpulse(n*90);

		// add the total
		if ( n == 0 )
			rho_combined = m_rho;
		else
			rho_combined += m_rho;
    }

    rho_combined /= 4;
	m_rho = rho_combined;
	zero_multqcs();
}

void spin_sys::crushgrad_ypulse(double angle)
{
    dcm rho_init = m_rho;

	dcm rho_combined;
    for ( int n = 0; n < 4; n++ )
    {
        if ( n > 0 )
            m_rho = rho_init;

        zpulse(n*90);

        ypulse(angle);

        zpulse(n*90);

		// add the total
		if ( n == 0 )
			rho_combined = m_rho;
		else
			rho_combined += m_rho;
    }

    rho_combined /= 4;
	m_rho = rho_combined;
	zero_multqcs();
}

void spin_sys::crushgrad_shaped_xpulse(drv& waveform, drv& phase, double angle, double t, double offset)
{
    dcm rho_init = m_rho;

	dcm rho_combined;
    for ( int n = 0; n < 4; n++ )
    {
        if ( n > 0 )
            m_rho = rho_init;

        zpulse(n*90);

        shaped_xpulse(waveform, phase, angle, t, offset);

        zpulse(n*90);

		// add the total
		if ( n == 0 )
			rho_combined = m_rho;
		else
			rho_combined += m_rho;
    }

    rho_combined /= 4;
	m_rho = rho_combined;
	zero_multqcs();
}

// offset is the center of the pulse in ppm
// TODO test angle threshold values
void spin_sys::profile_xpulse(drv& profile, drv& frequency, double angle, double offset)
{
    for ( int n = 0; n < m_chem_shift_vec.size(); n++ )
    {
        double tip_angle = get_tip_angle( profile, frequency, angle, m_chem_shift_vec(n) - offset );
        //std::cout << tip_angle << std::endl;
        //std::cout <<  m_chem_shift_vec(n) - offset << std::endl;
        if ( tip_angle > 1.0 )
        {
	        pulse(n, tip_angle, "x");
        }
    }
}
void spin_sys::profile_ypulse(drv& profile, drv& frequency, double angle, double offset)
{
    for ( int n = 0; n < m_chem_shift_vec.size(); n++ )
    {
        double tip_angle = get_tip_angle( profile, frequency, angle, m_chem_shift_vec(n) - offset );
        if ( tip_angle > 1.0 )
        {
	        pulse(n, tip_angle, "y");
        }
    }
}

void spin_sys::crushgrad_profile_ypulse(drv& profile, drv& frequency, double angle, double offset)
{
    dcm rho_init = m_rho;

	dcm rho_combined;
    for ( int n = 0; n < 4; n++ )
    {
        if ( n > 0 )
            m_rho = rho_init;

        zpulse(n*90);

        profile_ypulse(profile, frequency, angle, offset);

        zpulse(n*90);

		// add the total
		if ( n == 0 )
			rho_combined = m_rho;
		else
			rho_combined += m_rho;
    }

    rho_combined /= 4.0;
	m_rho = rho_combined;
	zero_multqcs();
}

void spin_sys::crushgrad_profile_xpulse(drv& profile, drv& frequency, double angle, double offset)
{
    dcm rho_init = m_rho;

	dcm rho_combined;
    for ( int n = 0; n < 4; n++ )
    {
        if ( n > 0 )
            m_rho = rho_init;

        zpulse(n*90);

        profile_xpulse(profile, frequency, angle, offset);

        zpulse(n*90);

		// add the total
		if ( n == 0 )
			rho_combined = m_rho;
		else
			rho_combined += m_rho;
    }

    rho_combined /= 4.0;
	m_rho = rho_combined;
	zero_multqcs();
}


void spin_sys::shaped_xpulse(drv& waveform, drv& phase, double angle, double t, double offset)
{
    int pts = waveform.size();

    double offset_freq = (offset - m_offset) * m_B0 / 1.0e6;
    if ( offset_freq > pts/t/2.0 )
        std::cout << "Warning, pulse offset freq exceeds pulse BW" << std::endl;

    double phase_factor = 180.0 / M_PI * (2.0 * M_PI * offset_freq );
    double phase_offset = t * offset_freq * 180.0; // so that phase is zero in middle of pulse

    double sum = 0;
    for ( int n = 0; n < pts; n++ )
        sum += waveform(n) * cos(phase(n) * M_PI / 180.0);

    dcm delay_a;
	dcm delay_b;
	get_delay_ops(t/pts, delay_a, delay_b);

    dcm half_delay_a;
	dcm half_delay_b;
	get_delay_ops(t/pts/2, half_delay_a, half_delay_b);

	apply_op(half_delay_a, half_delay_b);
    for ( int n = 0; n < pts ; n++ )
    {
	    xypulse(waveform(n)/sum*angle, phase(n) + phase_factor * t * n / ( pts - 1.0 ) - phase_offset );
	    apply_op(delay_a, delay_b);
    }
	apply_op(half_delay_a, half_delay_b);
}

void spin_sys::shaped_ypulse(drv& waveform, drv& phase, double angle, double t, double offset)
{
    int pts = waveform.size();
    
    double offset_freq = (offset - m_offset) * m_B0 / 1e6;
    if ( offset_freq > pts/t/2.0 )
        std::cout << "Warning, pulse offset freq exceeds pulse BW" << std::endl;

    double phase_factor = 180.0 / M_PI * (2.0 * M_PI * offset_freq );
    double phase_offset = t * offset_freq * 180.0; // so that phase is zero in middle of pulse

    double sum = 0;
    for ( int n = 0; n < pts; n++ )
        sum += waveform(n) * cos(phase(n) * M_PI / 180.0);
    
    dcm delay_a;
	dcm delay_b;
	get_delay_ops(t/pts, delay_a, delay_b);

    dcm half_delay_a;
	dcm half_delay_b;
	get_delay_ops(t/pts/2, half_delay_a, half_delay_b);

	apply_op(half_delay_a, half_delay_b);
    for ( int n = 0; n < pts ; n++ )
    {
	    xypulse(waveform(n)/sum*angle, 90 + phase(n) + phase_factor * t * n / ( pts - 1 ) - phase_offset);
	    apply_op(delay_a, delay_b);
    }
	apply_op(half_delay_a, half_delay_b);
}

void spin_sys::xpulse(double angle)
{
	pulse(angle, "x");
}

void spin_sys::ypulse(double angle)
{
	pulse(angle, "y");
}

void spin_sys::zpulse(double angle)
{
	pulse(angle, "z");
}

void spin_sys::xpulse(size_t n, double angle)
{
	pulse(n, angle, "x");
}

void spin_sys::ypulse(size_t n, double angle)
{
	pulse(n, angle, "y");
}

void spin_sys::zpulse(size_t n, double angle)
{
	pulse(n, angle, "z");
}

void spin_sys::xypulse(double angle, double phase)
{
	// get the operators
	dcm a;
	dcm b;
	get_xypulse_ops(angle, phase, a, b);
	
	// apply pulse to the spin system
	apply_op(a, b);
}

void spin_sys::pulse(double angle, std::string dirn)
{
	// get the operators
	dcm a;
	dcm b;
	get_pulse_ops(angle, dirn, a, b);
	
	// apply pulse to the spin system
	apply_op(a, b);
}

void spin_sys::pulse(size_t n, double angle, std::string dirn)
{
	// get the operators
	dcm a;
	dcm b;
	get_pulse_ops(n, angle, dirn, a, b);
	
	// apply pulse to the spin system
	apply_op(a, b);
}

void spin_sys::apply_op(dcm& a, dcm& b)
{
	m_rho = a * m_rho * b;
}

void spin_sys::get_pulse_ops(double angle, std::string dirn, dcm& a, dcm& b)
{
	// generate the F operator
	scm F_op;
	if ( dirn == "x" )
		Fx(F_op);
	if ( dirn == "y" )
		Fy(F_op);
	if ( dirn == "z" )
		Fz(F_op);

	// get a dense version
	dcm F_op_dense(F_op);

	std::complex<double> j(0,1); 

	a = (-F_op_dense * j * angle *M_PI/180.0).exp();
	b = (F_op_dense * j * angle *M_PI/180.0).exp();
}

void spin_sys::get_pulse_ops(size_t n, double angle, std::string dirn, dcm& a, dcm& b)
{
	// generate the F operator
	scm I_op;
	if ( dirn == "x" )
	    Ix(I_op, n);
	if ( dirn == "y" )
	    Iy(I_op, n);
	if ( dirn == "z" )
	    Iz(I_op, n);

	// get a dense version
	dcm I_op_dense(I_op);

	std::complex<double> j(0,1); 

	a = (-I_op_dense * j * angle *M_PI/180.0).exp();
	b = (I_op_dense * j * angle *M_PI/180.0).exp();
}

// convention is - phase of 0 is an x-pulse
void spin_sys::get_xypulse_ops(double angle, double phase, dcm& a, dcm& b)
{
	// generate the Fx and Fy operators
	scm F_op_xy;
    if ( cos( phase * M_PI/180.0 ) == 1 )       // pure x pulse
    {
        Fx(F_op_xy);
    }
    else if ( sin( phase * M_PI/180.0 ) == 1 )  // pure y pulse
    {
        Fy(F_op_xy);
    }
    else                                        // mixed xy pulse
    {
        scm F_op_x; scm F_op_y;
        Fx(F_op_x);
        Fy(F_op_y);

        // combine the operators depending on the phase of the pulse
        F_op_xy = F_op_x * cos(phase * M_PI/180.0) + 
                    F_op_y * sin(phase * M_PI/180.0);
    }

	// get a dense version
	dcm F_op_xy_dense(F_op_xy);
	std::complex<double> j(0,1); 

	a = (-F_op_xy_dense * j * angle *M_PI/180.0).exp();
	b = (F_op_xy_dense * j * angle *M_PI/180.0).exp();
}


void spin_sys::pulse_eig(double angle, std::string dirn)
{
	// generate the F operator
	scm F_op;
	if ( dirn == "x" )
		Fx(F_op);
	if ( dirn == "y" )
		Fy(F_op);

	// get a dense version
	dcm F_op_dense(F_op);

	// find the eign watsits
	Eigen::SelfAdjointEigenSolver<dcm> FSolver(F_op_dense);

	// find the inverse of the eigenvector matrix
	dcm eig_vec_inv = FSolver.eigenvectors().inverse();

	// create an eigenvalue matrix
	// left part of multiplication
	scm eig_val_mat_left(m_fsize, m_fsize);
	// construct the matrix exponential
	std::complex<double> j(0,1); 
	for ( size_t n = 0; n < m_fsize; n++ )
		eig_val_mat_left.coeffRef(n,n) = exp(-FSolver.eigenvalues()(n) * j * angle *M_PI/180.0);

	dcm mat_exp_left = FSolver.eigenvectors() * eig_val_mat_left * eig_vec_inv;

	// right part of multiplication
	scm eig_val_mat_right(m_fsize, m_fsize);

	// construct the matrix exponential
	for ( size_t n = 0; n < m_fsize; n++ )
		eig_val_mat_right.coeffRef(n,n) = exp(FSolver.eigenvalues()(n) * j * angle *M_PI/180.0);
	dcm mat_exp_right = FSolver.eigenvectors() * eig_val_mat_right * eig_vec_inv;

	// apply rotation to the system
	m_rho = mat_exp_left * m_rho * mat_exp_right;
}

void spin_sys::get_delay_ops(double time, dcm& a, dcm& b)
{
    if ( time < 0 )
    {
        std::cout << "Warning, negative delay requsted!!" << std::endl;
    }

	// find the inverse of the eigenvector matrix
	dcm eig_vec_inv = m_H_eig_vecs.inverse();

    // create an eigenvalue matrix
	// left part of multiplication
	scm eig_val_mat_left(m_fsize, m_fsize);

	// construct the matrix exponential
	std::complex<double> j(0,1); 
	for ( size_t n = 0; n < m_fsize; n++ )
		eig_val_mat_left.coeffRef(n,n) = exp(m_H_eig_vals(n) * j * time * 2.0 * M_PI );
	a = m_H_eig_vecs * eig_val_mat_left * eig_vec_inv;

	// right part of multiplication
	scm eig_val_mat_right(m_fsize, m_fsize);

	// construct the matrix exponential
	for ( size_t n = 0; n < m_fsize; n++ )
		eig_val_mat_right.coeffRef(n,n) = exp(-m_H_eig_vals(n) * j * time * 2.0 * M_PI );
	b = m_H_eig_vecs * eig_val_mat_right * eig_vec_inv;
}

void spin_sys::delay(double time)
{
	// get the operators
	dcm a;
	dcm b;
	get_delay_ops(time, a, b);
	
	// apply evolution to the spin system
	apply_op(a, b);
}

void spin_sys::set_state(const scm& prod_op)
{
    dcm prod_op_dense(prod_op);
    m_rho = prod_op_dense;
}

void spin_sys::set_state(const dcm& prod_op)
{
    m_rho = prod_op;
}

void spin_sys::get_state(dcm& prod_op)
{
	prod_op = m_rho;
}

void spin_sys::get_eig_vecs(dcm& eig_vecs)
{
    eig_vecs = m_H_eig_vecs;
}

void spin_sys::get_H(dcm& H)
{
    H = m_H;
}

void spin_sys::qn_states(drv& qn_states_vec)
{
    scm Fz_mat;
    Fz(Fz_mat);
	dcm Fz_mat_dense(Fz_mat);
	qn_states_vec = Fz_mat_dense.real().diagonal();
}

void spin_sys::qn_states(drm& qn_states_mat)
{
	drv qn_states_vec;
	qn_states(qn_states_vec);
	qn_states_mat.resize(m_fsize, m_fsize);
	qn_states_mat.setZero(m_fsize, m_fsize);
	for ( size_t n = 0; n < m_fsize; n++ )
		for ( size_t m = 0; m < m_fsize; m++ )
			qn_states_mat(n, m) = qn_states_vec(n) - qn_states_vec(m);

}

void spin_sys::zero_pqcs()
{
	drv qn_states_vec;
	qn_states(qn_states_vec);
	for ( size_t n = 0; n < m_fsize; n++ )
		for ( size_t m = 0; m < m_fsize; m++ )
			if ( qn_states_vec(n) - qn_states_vec(m) >= 1 )
				m_rho(n, m) = 0;
}

void spin_sys::zero_multqcs()
{
	drv qn_states_vec;
	qn_states(qn_states_vec);
	for ( size_t n = 0; n < m_fsize; n++ )
		for ( size_t m = 0; m < m_fsize; m++ )
			if ( qn_states_vec(n) - qn_states_vec(m) >= 2 )
				m_rho(n, m) = 0;
}


void spin_sys::get_peak_groups(std::vector<double>& freqs, std::vector<int>& freq_id, drv& group_vec)
{
	// first make sure the vectors are empty
	freqs.resize(0);
	freq_id.resize(0);

	dcm coherence = m_H_eig_vecs_diag.adjoint() * m_rho * m_H_eig_vecs_diag;
	
	for ( int n = 0; n < group_vec.maxCoeff(); n++ )
	{
		scm Ip_op(m_fsize, m_fsize);
		scm temp;
		for ( int m = 0; m < group_vec.size(); m++ )
        {
			if ( group_vec(m) - 1 == n )
			{
				Ip(temp, m);
				Ip_op += temp;
			}
        }
		
		//std::cout << Ip_op;
		dcm coupled_coherence = m_H_eig_vecs_diag.adjoint() * Ip_op * m_H_eig_vecs_diag;
		//std::cout << m_rho;

		std::complex<double> j(0,1); 
		std::complex<double> amp;
		double freq;

		// following can be used when all ppms are +ve
		// for (int p = 0+m+1; n < coherence.rows(); p++)

		for (int m = 0; m < coherence.rows(); m++)
        {
			for (int p = 0; p < coherence.rows(); p++)
			{
				// calculate resonance amplitudes
				amp = 2.0 * j * coherence(m,p) * coupled_coherence(m,p);
				amp /= coherence.rows();
				//std::cout << amp << std::endl;
				// ignore very small resonances
				if ( abs(amp) > 0.0001 )
				{
					// calculate resonance frequencies
					freq = -m_H_eig_vals_diag(p) + m_H_eig_vals_diag(m);
                    // following line was wrong in versions < 4.2.4
					//freq = -m_H_eig_vals(p) + m_H_eig_vals(m);
					freqs.push_back(freq);
					freq_id.push_back(n);
				}
			}
        }
	}
}

void spin_sys::acquire(dcm& time_sig, double fs, size_t N, double ref, double lambda, drv group_vec, double rec_phase, double delay)
{
	std::vector<double> freqs;
	std::vector<int> freq_id;	
	get_peak_groups(freqs, freq_id, group_vec);
    
    /*std::cout << freqs.size() << std::endl;
    for ( int n = 0; n < freqs.size(); n++ )
    {
        std::cout << freqs[n] << std::endl;
        std::cout << freq_id[n] << std::endl;
    }*/

    //std::cout << group_vec << std::endl;
    
    // resize time_sig according to the N parameter
    time_sig.resize(N,  size_t(group_vec.maxCoeff()));
	time_sig.setZero(N, size_t(group_vec.maxCoeff()));
    
    // check for no peaks
    if ( freq_id.size() == 0 )
        return;
	
	Eigen::ArrayXd freqs_eig = Eigen::ArrayXd::Map(&freqs[0], freqs.size());

	// find the diagonalised frequences for each group
    dcm coherence = m_H_eig_vecs.adjoint() * m_rho * m_H_eig_vecs;
    
    // diagonalised
	//dcm coherence = m_H_eig_vecs_diag.adjoint() * m_rho * m_H_eig_vecs_diag;

    scm Fp_op;
    Fp(Fp_op);
    dcm coupled_coherence = m_H_eig_vecs.adjoint() * Fp_op * m_H_eig_vecs;
    
    // diagonalised
    //dcm coupled_coherence = m_H_eig_vecs_diag.adjoint() * Fp_op * m_H_eig_vecs_diag;

    drv t(N);
    for (int n = 0; n < t.size(); n++)
        t(n) = n/fs + delay;
    
	// used for filtering out coherence orders
	//drv qn_states_vec;
	//qn_states(qn_states_vec);

    std::complex<double> j(0,1); 
    std::complex<double> amp;
    double freq;
	
	// following can be used when all ppms are +ve
    // for (int n = 0+m+1; n < coherence.rows(); n++)

    for (int m = 0; m < coherence.rows(); m++)
        for (int n = 0; n < coherence.rows(); n++)
        {
            // calculate resonance amplitudes
            amp = 2.0 * j * coherence(m,n) * coupled_coherence(m,n);
            amp /= coherence.rows();
            // ignore very small resonances
            // if ( ( abs(amp) > 0.0001 ) && ( qn_states_vec(m) - qn_states_vec(n) == 2.0 ) )
            if ( abs(amp) > 0.0001 )
            {
				//std::cout << qn_states_vec(m) - qn_states_vec(n) << "," << abs(amp) << std::endl;

                // calculate resonance frequencies
                freq = -m_H_eig_vals(n) + m_H_eig_vals(m);

                //freq = -m_H_eig_vals_diag(n) + m_H_eig_vals_diag(m);

				// find out which group it belongs to
				drv res = ( freqs_eig - freq ).abs();
				drv::Index min_ele;
				//double min = 
                res.minCoeff(&min_ele);
				int ind = freq_id[min_ele];

                // update the time_sig vector
                for (int p = 0; p < t.size(); p++)
                    time_sig(p,ind) += amp*exp(j*rec_phase*M_PI/180.0) * exp(t(p) * (j*2.0*M_PI*(freq + ref*m_B0/1e6)-lambda));
            }
        }

    // first data point correction
    time_sig(0) = time_sig(0) * 0.5;
}

void spin_sys::acquire(dcv& time_sig, double fs, size_t N, double ref, double lambda, double rec_phase, double delay)
{
    // resize time_sig according to the N parameter
    time_sig.resize(N);
	time_sig.setZero(N);
	
	// TODO the following two matrix multiplications are slooow
	// should we be using sparse matrices at this point?
	//scm m_H_eig_vecs_sp;
	//d2s(m_H_eig_vecs, m_H_eig_vecs_sp);

    dcm coherence = m_H_eig_vecs.adjoint() * m_rho * m_H_eig_vecs;

    scm Fp_op;
    Fp(Fp_op);

    dcm coupled_coherence = m_H_eig_vecs.adjoint() * Fp_op * m_H_eig_vecs;
	//Eigen::IOFormat CleanFmt(2, 0, ", ", "\n", "[", "]");
	//std::cout << m_H_eig_vecs.format(CleanFmt) << std::cout;

    drv t(N);
    for (int n = 0; n < t.size(); n++)
        t(n) = n/fs + delay;
    
    std::complex<double> j(0,1); 
    std::complex<double> amp;
    double freq;

	// following can be used when all ppms are +ve
    // for (int n = 0+m+1; n < coherence.rows(); n++)

    for (int m = 0; m < coherence.rows(); m++)
        for (int n = 0; n < coherence.rows(); n++)
        {
            // calculate resonance amplitudes
            amp = 2.0 * j * coherence(m,n) * coupled_coherence(m,n);
            amp /= coherence.rows();
            // ignore very small resonances
            if ( abs(amp) > 0.0001 )
            {
                // calculate resonance frequencies
                freq = -m_H_eig_vals(n) + m_H_eig_vals(m);
                // update the time_sig vector
                for (int n = 0; n < t.size(); n++)
                    time_sig(n) += amp*exp(j*rec_phase*M_PI/180.0) * exp(t(n) * (j*2.0*M_PI*(freq + ref*m_B0/1e6)-lambda));
            }
        }
    
    // first data point correction
    time_sig(0) = time_sig(0) * 0.5;
}

void spin_sys::Iz_spin(size_t n, scm& prod_op)
{
	Iz_pauli(m_spin_vec(n), prod_op);
}

void spin_sys::Ip_spin(size_t n, scm& prod_op)
{
	Ip_pauli(m_spin_vec(n), prod_op);
}

void spin_sys::Im_spin(size_t n, scm& prod_op)
{
	Im_pauli(m_spin_vec(n), prod_op);
}

void spin_sys::Ix_spin(size_t n, scm& prod_op)
{
	Ix_pauli(m_spin_vec(n), prod_op);
}

void spin_sys::Iy_spin(size_t n, scm& prod_op)
{
	Iy_pauli(m_spin_vec(n), prod_op);
}

void spin_sys::F(scm& prod_op, std::string op)
{
	I(prod_op, 0, op);
	scm temp;
	for ( int n = 1; n < m_spin_vec.size(); n++ )
	{
		I(temp, n, op);
        prod_op += temp;
	}
}

void spin_sys::F(scm& prod_op, std::string op, dcv& spin_num_vec)
{
    // TODO have an if statement here to create an empty matrix if spin_num_vec == 0

	I(prod_op, 0, op);
    prod_op *= spin_num_vec(0);

	for ( int n = 1; n < m_spin_vec.size(); n++ )
	{
        if ( std::abs(spin_num_vec(n)) != 0 )
        {
	        scm temp;
            I(temp, n, op);
            prod_op += temp * spin_num_vec(n);
        }
	}
}

void spin_sys::Fz(scm& prod_op)
{
	F(prod_op, "z");
}

void spin_sys::Fx(scm& prod_op)
{
	F(prod_op, "x");
}

void spin_sys::Fy(scm& prod_op)
{
	F(prod_op, "y");
}

void spin_sys::Fp(scm& prod_op)
{
	F(prod_op, "p");
}

void spin_sys::Fm(scm& prod_op)
{
	F(prod_op, "m");
}

void spin_sys::Fz(scm& prod_op, dcv& spin_num_vec)
{
	F(prod_op, "z", spin_num_vec);
}

void spin_sys::Fx(scm& prod_op, dcv& spin_num_vec)
{
	F(prod_op, "x", spin_num_vec);
}

void spin_sys::Fy(scm& prod_op, dcv& spin_num_vec)
{
	F(prod_op, "y", spin_num_vec);
}

void spin_sys::Fp(scm& prod_op, dcv& spin_num_vec)
{
	F(prod_op, "p", spin_num_vec);
}

void spin_sys::Fm(scm& prod_op, dcv& spin_num_vec)
{
	F(prod_op, "m", spin_num_vec);
}

void spin_sys::Iz(scm& prod_op, size_t n)
{
	I(prod_op, n, "z");
}

void spin_sys::Ix(scm& prod_op, size_t n)
{
	I(prod_op, n, "x");
}

void spin_sys::Iy(scm& prod_op, size_t n)
{
	I(prod_op, n, "y");
}

void spin_sys::Ip(scm& prod_op, size_t n)
{
	I(prod_op, n, "p");
}

void spin_sys::Im(scm& prod_op, size_t n)
{
	I(prod_op, n, "m");
}

void spin_sys::I(scm& prod_op, size_t n, std::string op)
{
	scm In;
	if ( op == "z" )
		Iz_spin(n, In);
	else if ( op == "x" )
		Ix_spin(n, In);
	else if ( op == "y" )
		Iy_spin(n, In);
	else if ( op == "p" )
		Ip_spin(n, In);
	else if ( op == "m" )
		Im_spin(n, In);
	
	// if this is the left most spin only one kron is required
	if (n == 0)
	{
		size_t dim = size_t( m_fsize / (floor(2 * m_spin_vec(n) + 1 + 0.5)) );
		scm I(dim, dim);
		for ( int m = 0; m < I.rows(); m++ )
			I.coeffRef(m,m) = 1;
		kron(In, I, prod_op);
	}
	// if this is the right most spin only one kron is required
	else if ( n == size_t( m_spin_vec.size() - 1 ) )
	{
		size_t dim = size_t( m_fsize / (floor(2 * m_spin_vec(n) + 1 + 0.5)) );
		scm I(dim, dim);
		for ( int m = 0; m < I.rows(); m++ )
			I.coeffRef(m,m) = 1;
		kron(I, In, prod_op);
	}
	// otherwise two krons are required
	else
	{
		size_t lsize = 1;
		for ( size_t m = 0; m < n; m++ )
			lsize = lsize * size_t(2 * m_spin_vec(m) + 1 + 0.5);

		size_t rsize = 1;
		for ( int m = n + 1; m < m_spin_vec.size(); m++ )
			rsize = rsize * size_t(2 * m_spin_vec(m) + 1 + 0.5);

		scm l_I(lsize, lsize);
		for ( int m = 0; m < l_I.rows(); m++ )
			l_I.coeffRef(m,m) = 1;

		scm r_I(rsize, rsize);
		for ( int m = 0; m < r_I.rows(); m++ )
			r_I.coeffRef(m,m) = 1;

		// do the left part of the kron
		scm temp; 
		kron(l_I, In, temp);

		// do the right part of the kron
		kron(temp, r_I, prod_op); 
	}
}

void spin_sys::H(dcm& prod_op)
{
    scm H_op;	
	H(H_op);
    dcm H_dense(H_op);
    //std::ofstream mat_file("H_mat.txt");
    //mat_file << "Sparse matrix" << std::endl;
    //mat_file << H_op << std::endl;
    //mat_file << "Dense matrix" << std::endl;
    //mat_file << H_dense << std::endl;
    prod_op = H_dense;
}

void spin_sys::change_offset(double offset)
{
    scm Fz_mat;
    Fz(Fz_mat);
    m_H += Fz_mat * ( -m_offset + offset) * m_B0 * 1e-6;
    m_offset = offset;
    
    // perform a symmetric eigenvalue decomposition
    Eigen::SelfAdjointEigenSolver<dcm> Hsolver(m_H);

    // assign the member variables 
    m_H_eig_vals = Hsolver.eigenvalues();
    m_H_eig_vals_mat.resize(m_fsize, m_fsize);
	m_H_eig_vals_mat.setZero(m_fsize, m_fsize);
    m_H_eig_vals_mat.diagonal() = m_H_eig_vals;
    m_H_eig_vecs = Hsolver.eigenvectors();

	// for the weak coupling case
	m_H_eig_vals_diag = m_H.diagonal().real();
    m_H_eig_vals_diag_mat.resize(m_fsize, m_fsize);
	m_H_eig_vals_diag_mat.setZero(m_fsize, m_fsize);
    m_H_eig_vals_diag_mat.diagonal() = m_H_eig_vals_diag;
	m_H_eig_vecs_diag = dcv::Constant(m_fsize, 1).asDiagonal();
}

void spin_sys::H(scm& prod_op)
{
	// set prod_op to be zero and the correct size
	prod_op = scm(m_fsize, m_fsize);

	// chemical shift part
	for ( int n = 0; n < m_spin_vec.size(); n++ )
	{
		scm Izn;
		Iz(Izn, n);
		prod_op += Izn * ( (-m_chem_shift_vec(n) + m_offset ) * m_B0 * 1e-6);

        //m_chem_shift_vec = -chem_shift_vec.array() * B0 * 1e-6;
	}

	// J-coupling part (note this involved matrix matrix multiply
	// which can be costly if not done with sparse matrices
	scm tempn;
	scm tempm;

	for ( int n = 0; n < m_spin_vec.size(); n++ )
	{
		for ( int m = n; m < m_spin_vec.size(); m++ )
		{
			if ( m_j_coupling_mat(n,m) != 0.0 )
			{
				std::complex<double> j_coup (m_j_coupling_mat(n,m), 0);
				// x part
				I(tempn, n, "x");
				I(tempm, m, "x");
				prod_op += j_coup * tempn * tempm;
				// y part
				I(tempn, n, "y");
				I(tempm, m, "y");
				prod_op += j_coup * tempn * tempm;
				// z part
				I(tempn, n, "z");
				I(tempm, m, "z");
				prod_op += j_coup * tempn * tempm;
			}
		}
	}
}

// Sparse Kronecker product function
// A kron B = C is calculated using the following maps, for 1 based indices:
// A(r1, c1) * B(r2, c2) = C( (r1-1)*B.rows() + r2, (c1-1)*B.cols() + c2 )
// for 0 based indices:
// A(r1, c1) * B(r2, c2) = C( r1*B.rows() + r2, c1*B.cols() + c2 )
void kron(const scm& A, const scm& B, scm& C)
{
	// resize C based on A and B1
	size_t rows = A.rows() * B.rows();
	size_t cols = A.cols() * B.cols();
	
	// clear C and set it to the correct size
	C = scm(rows, cols);

	// cycle though non-zero elements of A and B
	for (int k = 0; k < A.outerSize(); k++)
		for (scm::InnerIterator itA(A,k); itA; ++itA)
			for (int l = 0; l < B.outerSize(); l++)
				for (scm::InnerIterator itB(B,l); itB; ++itB)
				{
					// and calculate corresponding element of C
					C.coeffRef( itA.row()*B.rows() + itB.row(), itA.col()*B.cols() + itB.col() )
						= itA.value() * itB.value();
				}
}
