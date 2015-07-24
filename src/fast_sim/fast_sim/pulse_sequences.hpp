#ifndef __PULSESEQ__
#define __PULSESEQ__

#include "fast_sim.hpp"

void pulse_acquire(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat);

void spin_echo(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double tau);

void press(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double TE1, double TE2);

void semi_laser(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double t1, double t2, double t3, double t4);

void profile_press(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double TE1, double TE2);

void shaped_press(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double TE1, double TE2);

void mega_press(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double TE1, double TE2);

void steam(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double TE, double TM);

void laser(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double tau);

void cpmg(drv& spin_vec, drv& chem_shift_vec, drm& j_coupling_mat, drv& group_vec, dcv& spin_num_vec, double B0, double fs, size_t N, double ref, double lambda, dcm& time_sig_mat, double tau, int cpmg_pulses);

#endif
