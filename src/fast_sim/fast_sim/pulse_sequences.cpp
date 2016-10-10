#include <iostream>
#include <fstream>
#include "pulse_sequences.hpp"
#include "fast_sim.hpp"

void pulse_acquire(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double delay)
{
    // initialise the spin system
	spin_sys sys(spin_vec, chem_shift_vec, j_coupling_mat, B0);

    // set the state of the system to be -Fy
    scm rho;
    sys.Fy(rho, spin_num_vec);
    rho *= -1;
	sys.set_state(rho);
	
	// acquire the result	
	sys.acquire(time_sig_mat, fs, N, ref, lambda, group_vec, 180, delay);
}

void spin_echo(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double tau)
{
    // initialise the spin system
	spin_sys sys(spin_vec, chem_shift_vec, j_coupling_mat, B0);

    // set the state of the system to be -Fy
    scm rho;
    sys.Fy(rho, spin_num_vec);
    rho *= -1;
	sys.set_state(rho);
	
	dcm a;
	dcm b;
	sys.get_delay_ops(tau/2.0, a, b);
	sys.apply_op(a, b);
    sys.ypulse(180);
	sys.apply_op(a, b);
	
	// acquire the result	
	sys.acquire(time_sig_mat, fs, N, ref, lambda, group_vec);
}

void press(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double TE1, double TE2)
{
    // initialise the spin system
	spin_sys sys(spin_vec, chem_shift_vec, j_coupling_mat, B0);

    // set the state of the system to be -Fy
    scm rho;
    sys.Fy(rho, spin_num_vec);
    rho *= -1;
	sys.set_state(rho);

	sys.delay(TE1/2.0);
	dcm a;
	dcm b;
    sys.get_pulse_ops(180, "y", a, b);
	sys.apply_op(a, b);
	sys.delay((TE1+TE2)/2.0);
	sys.apply_op(a, b);
	sys.delay(TE2/2.0);

	// acquire the result	
	sys.acquire(time_sig_mat, fs, N, ref, lambda, group_vec);
}

void semi_laser(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double t1, double t2, double t3, double t4)
{
    // initialise the spin system
	spin_sys sys(spin_vec, chem_shift_vec, j_coupling_mat, B0);

    // set the state of the system to be -Fy
    scm rho;
    sys.Fy(rho, spin_num_vec);
    rho *= -1;
	sys.set_state(rho);

	sys.delay(t1);
	dcm a;
	dcm b;
    sys.get_pulse_ops(180, "y", a, b);
	sys.apply_op(a, b);
	sys.delay(t2);
	sys.apply_op(a, b);
	sys.delay(t3);
	sys.apply_op(a, b);
	sys.delay(t4);
	sys.apply_op(a, b);
	sys.delay(t4-t3+t2-t1);

	// acquire the result	
	sys.acquire(time_sig_mat, fs, N, ref, lambda, group_vec);
}


void shaped_press(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double TE1, double TE2)
{
    // initialise the spin system
	spin_sys sys(spin_vec, chem_shift_vec, j_coupling_mat, B0);

    // set the state of the system to be -Fy
    scm rho;
    sys.Fy(rho, spin_num_vec);
    rho *= -1;
	sys.set_state(rho);

    drv waveform;
    drv phase;
    read_points("./Philips_style_180/pulse_shape.txt", waveform, phase);
    double pul_t = 6.912e-3;
    double bw_ppm = 10.0;
    double pl180 = pul_t;
    double pl90 = 0;

    dcm rho_init;
    sys.get_state(rho_init);
    dcm rho_combined;
    std::vector<dcm> rho_vec;
    
    double max_ppm = 4.1;
    double min_ppm = 1.3;
    double carrier_freq = ( max_ppm + min_ppm ) / 2.0;
    double full_bw = bw_ppm * 1.25;

    size_t x_pts = 32;
    double x_factor = (max_ppm + full_bw/2.0)*2.0/x_pts;
    size_t y_pts = 32;
    double y_factor = (max_ppm + full_bw/2.0)*2.0/y_pts;

    bool first_pass = true;

    double x_offset;
    double y_offset;

    for ( size_t y_pos = 0; y_pos < y_pts ; y_pos++ )
    {
        for ( size_t x_pos = 0; x_pos < x_pts ; x_pos++ )
        {
            if (!first_pass)
                sys.set_state(rho_init);

            x_offset = carrier_freq+(x_pos-(x_pts-1.0)/2.0)*x_factor;
            y_offset = carrier_freq+(y_pos-(y_pts-1.0)/2.0)*y_factor;

            //sys.xpulse(90);
            sys.delay(TE1/2.0-pl90/2.0-pl180/2.0);
            sys.crushgrad_shaped_xpulse(waveform, phase, 180, pl180, x_offset);
            sys.delay(TE1/2.0+TE2/2.0-pl180);
            sys.crushgrad_shaped_xpulse(waveform, phase, 180, pl180, y_offset);
            sys.delay(TE2/2.0-pl180/2.0);

            dcm rho_temp;
            sys.get_state(rho_temp);

            rho_vec.push_back(rho_temp);

            if ( first_pass )
                rho_combined = rho_temp;
            else
                rho_combined += rho_temp;
            
            first_pass = false;

            if ( x_pos == 0 ) 
                std::cout << "Location " << x_pos+y_pos*(y_pts) + 1 << " of " << x_pts*y_pts << std::endl;
        }
    }

    sys.set_state(rho_combined);

	// acquire the result	
	sys.acquire(time_sig_mat, fs, N, ref, lambda, group_vec);

    time_sig_mat *= 0.25/57.6686; // rescale

    /*
    // and plot
    dcv time_sig = time_sig_mat.col(0);
    // get the spectrum
    drv freq_sig;
    get_spectrum(time_sig, freq_sig);
    // get the frequency scale
    drv freq_ppm;
    sys.gen_ppm_scale(freq_ppm,N,fs,ref);
    plot(freq_ppm, freq_sig, 4.5, 0.5);
    */
}


void profile_press(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double TE1, double TE2)
{
    // initialise the spin system
	spin_sys sys(spin_vec, chem_shift_vec, j_coupling_mat, B0);

    // set the state of the system to be -Fy
    scm rho;
    sys.Fy(rho, spin_num_vec);
    rho *= -1;
	sys.set_state(rho);

    drv pul_freq;
    drv pul_profile;
    read_points("./Philips_style_180/pulse_profile.txt", pul_freq, pul_profile);
    double bw_ppm = 10.0;

    dcm rho_init;
    sys.get_state(rho_init);
    dcm rho_combined;
    std::vector<dcm> rho_vec;
    
    double max_ppm = 4.1;
    double min_ppm = 1.3;
    double carrier_freq = ( max_ppm + min_ppm ) / 2.0;
    double full_bw = bw_ppm * 1.25;

    size_t x_pts = 32;
    double x_factor = (max_ppm + full_bw/2.0)*2.0/x_pts;
    size_t y_pts = 32;
    double y_factor = (max_ppm + full_bw/2.0)*2.0/y_pts;

    drv projection(x_pts*y_pts);
    drm projection_mat(x_pts,y_pts);

    bool first_pass = true;

    double x_offset;
    double y_offset;

    for ( size_t y_pos = 0; y_pos < y_pts ; y_pos++ )
    {
        for ( size_t x_pos = 0; x_pos < x_pts ; x_pos++ )
        {
            if (!first_pass)
                sys.set_state(rho_init);

            x_offset = carrier_freq+(x_pos-(x_pts-1.0)/2.0)*x_factor;
            y_offset = carrier_freq+(y_pos-(y_pts-1.0)/2.0)*y_factor;

            //sys.xpulse(90);
            sys.delay(TE1/2.0);
            sys.crushgrad_profile_xpulse(pul_profile,pul_freq,180,x_offset);
            sys.delay((TE1+TE2)/2.0);
            sys.crushgrad_profile_xpulse(pul_profile,pul_freq,180,y_offset);
            sys.delay(TE2/2.0);

            dcm rho_temp;
            sys.get_state(rho_temp);

            rho_vec.push_back(rho_temp);

            if ( first_pass )
                rho_combined = rho_temp;
            else
                rho_combined += rho_temp;
            
            first_pass = false;
        }
    }

    sys.set_state(rho_combined);

	// acquire the result	
	sys.acquire(time_sig_mat, fs, N, ref, lambda, group_vec);
    
    time_sig_mat *= 0.25/59.3828; // rescale

    /*
    // and plot
    dcv time_sig = time_sig_mat.col(0);
    // get the spectrum
    drv freq_sig;
    get_spectrum(time_sig, freq_sig);
    // get the frequency scale
    drv freq_ppm;
    sys.gen_ppm_scale(freq_ppm,N,fs,ref);
    plot(freq_ppm, freq_sig, 4.5, 0.5);
    */
}

void mega_press(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double TE1, double TE2)
{
    // initialise the spin system
	spin_sys sys_ed_off(spin_vec, chem_shift_vec, j_coupling_mat, B0);

    // set the state of the system to be -Fy
    scm rho;
    sys_ed_off.Fy(rho, spin_num_vec);
    rho *= -1;
	sys_ed_off.set_state(rho);

	sys_ed_off.delay(TE1/2.0);
	dcm a;
	dcm b;
    sys_ed_off.get_pulse_ops(180, "y", a, b);
    // apply 180
	sys_ed_off.apply_op(a, b);
	sys_ed_off.delay((TE1+TE2)/2.0);
    // apply 180
	sys_ed_off.apply_op(a, b);
	sys_ed_off.delay(TE2/2.0);

	spin_sys sys_ed_on(spin_vec, chem_shift_vec, j_coupling_mat, B0);

    // set the state of the system to be -Fy
	sys_ed_on.set_state(rho);
	sys_ed_on.delay(TE1/2.0);
    // apply 180
	sys_ed_on.apply_op(a, b);
	sys_ed_on.delay((TE1+TE2)/4.0);
    for ( int n = 0; n < spin_vec.size(); n++ )
        if ( ( chem_shift_vec(n) > 1.6 ) && ( chem_shift_vec(n) < 2.2 ) )
            sys_ed_on.ypulse(n, 180);

	sys_ed_on.delay((TE1+TE2)/4.0);
    // apply 180
	sys_ed_on.apply_op(a, b);
	sys_ed_on.delay(TE2/4.0);
    // edit pulse
	for ( int n = 0; n < spin_vec.size(); n++ )
        if ( ( chem_shift_vec(n) > 1.6 ) && ( chem_shift_vec(n) < 2.2 ) )
            sys_ed_on.ypulse(n, 180);

    sys_ed_on.delay(TE2/4.0);
	// acquire the result	
    dcm edit_off;
	sys_ed_off.acquire(edit_off, fs, N, ref, lambda, group_vec);
    dcm edit_on;
	sys_ed_on.acquire(edit_on, fs, N, ref, lambda, group_vec);
    time_sig_mat = edit_on-edit_off;
}


void steam(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double TE, double TM)
{
    // initialise the spin system
	spin_sys sys(spin_vec, chem_shift_vec, j_coupling_mat, B0);
	
    // set the state of the system to be Fz
    scm rho;
    sys.Fz(rho, spin_num_vec);
	sys.set_state(rho);

	dcm rho_init;
	sys.get_state(rho_init);
	// set the number of dephasers
	int dephasers = 4;
	dcm rho_combined;

	dcm delay_half_TE_a;
	dcm delay_half_TE_b;
	sys.get_delay_ops(TE/2.0, delay_half_TE_a, delay_half_TE_b);
	
	dcm delay_TM_a;
	dcm delay_TM_b;
	sys.get_delay_ops(TM, delay_TM_a, delay_TM_b);

	dcm pulse_90x_a;
	dcm pulse_90x_b;
	sys.get_pulse_ops(90, "x", pulse_90x_a, pulse_90x_b);

	for ( int n = 0; n < dephasers; n++ )
	{
		sys.set_state(rho_init);
		double phase = n * 360.0 / dephasers;

		// generate the first and third 90 pulse ops
		dcm pulse_phase_a;
		dcm pulse_phase_b;
		sys.get_xypulse_ops(90, phase, pulse_phase_a, pulse_phase_b);

		// first 90
		sys.apply_op(pulse_phase_a, pulse_phase_b);
		// evolve TM/2
		sys.apply_op(delay_half_TE_a, delay_half_TE_b);
		// second 90
		sys.apply_op(pulse_90x_a, pulse_90x_b);
		// evolve TM/2
		sys.apply_op(delay_TM_a, delay_TM_b);
		// zero coherences >= 1
		sys.zero_pqcs();
		// third 90
		sys.apply_op(pulse_phase_a, pulse_phase_b);
		// evolve TE/2
		sys.apply_op(delay_half_TE_a, delay_half_TE_b);

		dcm rho_temp;
		sys.get_state(rho_temp);
		// add the the total
		if ( n == 0 )
			rho_combined = rho_temp;
		else
			rho_combined += rho_temp;
	}
	// scale depending on the number of dephasers
	rho_combined /= dephasers;
	sys.set_state(rho_combined);

	// acquire the result	
	sys.acquire(time_sig_mat, fs, N, ref, lambda, group_vec);

	// correct the phase
	time_sig_mat = -time_sig_mat;
}

void laser(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double tau)
{
    // initialise the spin system
	spin_sys sys(spin_vec, chem_shift_vec, j_coupling_mat, B0);

    // set the state of the system to be -Fy
    scm rho;
    sys.Fy(rho, spin_num_vec);
    rho *= -1;
    
	sys.set_state(rho);

	dcm delay12_a;
	dcm delay12_b;
	sys.get_delay_ops(tau/12.0, delay12_a, delay12_b);
	dcm delay6_a;
	dcm delay6_b;
	sys.get_delay_ops(tau/6.0, delay6_a, delay6_b);
	dcm pulse180_a;
	dcm pulse180_b;
    sys.get_pulse_ops(180, "y", pulse180_a, pulse180_b);
	
	// evolve tau/12
	sys.apply_op(delay12_a, delay12_b);
	// 180
	sys.apply_op(pulse180_a, pulse180_b);
	for ( int n = 0; n < 5; n++ )
	{
		// evolve tau/6
		sys.apply_op(delay6_a, delay6_b);
		// 180
		sys.apply_op(pulse180_a, pulse180_b);
	}
	// evolve tau/12
	sys.apply_op(delay12_a, delay12_b);

	// acquire the result	
	sys.acquire(time_sig_mat, fs, N, ref, lambda, group_vec);
}

void cpmg(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double tau, int cpmg_pulses)
{
    // initialise the spin system
	spin_sys sys(spin_vec, chem_shift_vec, j_coupling_mat, B0);

    // set the state of the system to be -Fy
    scm rho;
    sys.Fy(rho, spin_num_vec);
    rho *= -1;
	sys.set_state(rho);

	dcm half_delay_a;
	dcm half_delay_b;
	sys.get_delay_ops(tau/(2*cpmg_pulses), half_delay_a, half_delay_b);
	dcm delay_a;
	dcm delay_b;
	sys.get_delay_ops(tau/cpmg_pulses, delay_a, delay_b);
	dcm pulse180_a;
	dcm pulse180_b;
    sys.get_pulse_ops(180, "y", pulse180_a, pulse180_b);
	
	// evolve tau/(2*cpmg_pulses)
	sys.apply_op(half_delay_a, half_delay_b);
	// 180 pulse
	sys.apply_op(pulse180_a, pulse180_b);
	for ( int n = 0; n < cpmg_pulses - 1; n++ )
	{
		// evolve tau/(cpmg_pulses)
		sys.apply_op(delay_a, delay_b);
		// 180 pulse
		sys.apply_op(pulse180_a, pulse180_b);
	}
	// evolve tau/(2*cpmg_pulses)
	sys.apply_op(half_delay_a, half_delay_b);
	
	// acquire the result	
	sys.acquire(time_sig_mat, fs, N, ref, lambda, group_vec);
}

