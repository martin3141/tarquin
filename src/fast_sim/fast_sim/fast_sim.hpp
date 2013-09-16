#ifndef __FASTSIM__
#define __FASTSIM__

//#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
//#include <cvm.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions> // needed for matrix exponential
#include <unsupported/Eigen/FFT>
#include <complex>

typedef Eigen::MatrixXcd dcm;
typedef Eigen::MatrixXd drm;
typedef Eigen::VectorXcd dcv;
typedef Eigen::VectorXd drv;
typedef Eigen::SparseMatrix<std::complex<double> > scm;
//typedef Eigen::DynamicSparseMatrix<std::complex<double> > scm;

void d2s(dcm& in, scm& out);

// generate the pauli spin matrices
void Iz_pauli(double I, scm& prod_op);
void Ip_pauli(double I, scm& prod_op);
void Im_pauli(double I, scm& prod_op);
void Ix_pauli(double I, scm& prod_op);
void Iy_pauli(double I, scm& prod_op);
void kron(const scm& A, const scm& B, scm& C);

// pulses
void gen_sinc(drv& waveform, drv& phase, int n, int pts, int zp_pts = 0);
void gen_gaus(drv& waveform, drv& phase, double trunc_fact, int pts, int zp_pts = 0);
void gen_ham(drv& waveform, int pts, int zp_pts);
void gen_sinc_ham(drv& waveform, drv& phase, int n, int pts, int zp_pts = 0);

void read_points(std::string fname, drv& x, drv& y);
void write_points(std::string fname, drv& x, drv& y);
void write_array(std::string fname, drv& x);

void plot(drv& x, drv& y, double xmin, double xmax, std::string title = "");
void plot(drv& x, drv& y, std::string title = "");
void plot(drv& x, std::string title = "");
void plot_pdf(std::string fname, drv& x, drv& y, double xmin, double xmax, std::string title = "");
void plot_pdf(std::string fname, drv& x, std::string title = "");

void get_spectrum(dcv& time_sig, drv& freq_sig, int type = 1);

class spin_sys {
	public:

	// constructor
	spin_sys( const drv& spin_vec, const drv& chem_shift_vec, 
		const drm& j_coupling_mat, double B0, double offset = 0);
	
    void Iz_spin(size_t n, scm& prod_op);
	void Ip_spin(size_t n, scm& prod_op);
	void Im_spin(size_t n, scm& prod_op);
	void Ix_spin(size_t n, scm& prod_op);
	void Iy_spin(size_t n, scm& prod_op);
    
	void I(scm& prod_op, size_t n, std::string op);
	void Iz(scm& prod_op, size_t n);
	void Ip(scm& prod_op, size_t n);
	void Im(scm& prod_op, size_t n);
	void Ix(scm& prod_op, size_t n);
	void Iy(scm& prod_op, size_t n);

	void F(scm& prod_op, std::string op);
    void Fz(scm& prod_op);
    void Fp(scm& prod_op);
    void Fm(scm& prod_op);
    void Fx(scm& prod_op);
    void Fy(scm& prod_op);

	void F(scm& prod_op, std::string op, drv& spin_num_vec);
	void Fz(scm& prod_op, drv& spin_num_vec);
    void Fp(scm& prod_op, drv& spin_num_vec);
    void Fm(scm& prod_op, drv& spin_num_vec);
    void Fx(scm& prod_op, drv& spin_num_vec);
    void Fy(scm& prod_op, drv& spin_num_vec);

	void H(scm& prod_op);
    void H(dcm& prod_op);

    void set_state(const scm& prod_op);
    void set_state(const dcm& prod_op);

    void get_state(dcm& prod_op);
    void get_eig_vecs(dcm& eig_vecs);
    void get_H(dcm& H);

	void qn_states(drv& qn_states_vec);
	void qn_states(drm& qn_states_mat);
	void zero_pqcs();
	void zero_multqcs();
	void apply_op(dcm& a, dcm& b);

    void acquire(dcv& time_sig, double fs, size_t N, double ref, double lambda, double rec_phase = 180);
	void get_peak_groups(std::vector<double>& freqs, std::vector<int>& freq_id, drv& group_vec);

    void acquire(dcm& time_sig, double fs, size_t N, double ref, double lambda, drv group_vec, double rec_phase = 180);

    void shaped_xpulse(drv& waveform, drv& phase, double angle, double t, double offset = 0);
    void shaped_ypulse(drv& waveform, drv& phase, double angle, double t, double offset = 0);

    void profile_xpulse(drv& profile, drv& frequency, double angle, double offset);
    void profile_ypulse(drv& profile, drv& frequency, double angle, double offset);

    void crushgrad_shaped_xpulse(drv& waveform, drv& phase, double angle, double t, double offset = 0);
    void crushgrad_shaped_ypulse(drv& waveform, drv& phase, double angle, double t, double offset = 0);

    void crushgrad_xpulse(double angle);
    void crushgrad_ypulse(double angle);

    void crushgrad_profile_xpulse(drv& profile, drv& frequency, double angle, double offset);
    void crushgrad_profile_ypulse(drv& profile, drv& frequency, double angle, double offset);

    void gen_fft_profile(drv& waveform, int pul_N, double pul_t, drv& pul_freq, drv& pul_profile);

    void gen_freq_scale(drv& freq_scale, int N, double fs);

    void gen_ppm_scale(drv& ppm_scale, int N, double fs, double offset);

    double get_tip_angle(drv& profile, drv& freq_scale, double angle, double frequency);

	void xpulse(double angle);
    void xpulse(size_t n, double angle);
	void ypulse(double angle);
	void ypulse(size_t n, double angle);
	void zpulse(double angle);
	void zpulse(size_t n, double angle);
	void pulse_eig(double angle, std::string dirn);
	void pulse(double angle, std::string dirn);
    void pulse(size_t n, double angle, std::string dirn);
	void xypulse(double angle, double phase);
	void get_pulse_ops(double angle, std::string dirn, dcm& a, dcm& b);
	void get_xypulse_ops(double angle, double phase, dcm& a, dcm& b);
    void get_pulse_ops(size_t n, double angle, std::string dirn, dcm& a, dcm& b);

    void delay(double time);
	void get_delay_ops(double time, dcm& a, dcm& b);
    void change_offset(double offset);

	private:
	// Number of spins in the system
	size_t m_nspins;
    // Vector of nuclear spin numbers
	drv m_spin_vec;
    // Vector of chemical shift values
	drv m_chem_shift_vec;
    // Matrix of j-coupling network
	drm m_j_coupling_mat;
    // Magnetic field strength
	double m_B0;
    // Offset in PPM	
    double m_offset;
    // Matrix size of the spin ops
    size_t m_fsize;
    // Dense representation of the Hamiltonian
    dcm m_H;
    // H eigen values
    drv m_H_eig_vals;
    drm m_H_eig_vals_mat;
    // H eigen vectors
    dcm m_H_eig_vecs;
    // H eigen values weak coupling
    drv m_H_eig_vals_diag;
    drm m_H_eig_vals_diag_mat;
	// H eigen vectors weak coupling
    dcm m_H_eig_vecs_diag;
	// the current state of the system
    dcm m_rho;
};

#endif
