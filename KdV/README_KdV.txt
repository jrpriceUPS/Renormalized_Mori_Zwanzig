This folder contains Matlab files to derive and simulate reduced order models for the Kortweg-de Vries equation with small dispersion:

u_t + epsilon u_{xxx} + alpha * uu_x = 0.

%%%%%%%%%%%%%%%%%%%%%%
%simulation_functions% 
%%%%%%%%%%%%%%%%%%%%%%
Files necessary for running a KdV ROM simulation

full_init_KdV: takes simulation_params structure and updates it to include necessary parameters for running a full model of KdV

complete_init_KdV: takes simulation_params structure and updates it to include necessary parameters for running a complete memory approximation ROM of KdV

fft_norm: computes the FFT of a vector using the standard normalization

ifft_norm: computes the IFFT of a vector using the standard normalization

convolution_sum_KdV: computes the nonlinear part of KdV in real space and transforms back

PDE_solve: solves KdV with given parameters

RK4_stiff_nonstiff_step: computes a single step in an integration scheme known to function well for KdV with small dispersion

BKdV_files: files used for simulating Burgers-KdV (a new problem I have not yet completed)

old_files: out-of-date files used in past iterations of the repository



%%%%%%%%%%%
%nonlinear% 
%%%%%%%%%%%
Files associated with the RHS of KdV full and ROM simulations

nonlinear_full_KdV: computes the RHS of KdV for a full model

renormalized_complete_2nd_order_KdV: computes the RHS for an ROM with only first and second order terms of the complete memory approximation

renormalized_complete_4th_order_KdV: computes the RHS for an ROM with the first four terms of the complete memory approximation

markov_term_KdV: computes the Markov term of KdV

tmodel_term_KdV: computes the first order complete memory approximation term for KdV

t2model_term_complete_KdV: computes the second order complete memory approximation term for KdV

t3model_term_complete_KdV: computes the third order complete memory approximation term for KdV

t4model_term_complete_KdV: computes the fourth order complete memory approximation term for KdV

BKdV_files: files used for simulating Burgers-KdV (a new problem I have not yet completed)

old_files: out-of-date files used in past iterations of the repository



%%%%%%%%%%
%analysis% 
%%%%%%%%%%
Files associated with analysis of KdV results

get_energy: computes the energy in a subset of modes at all times from a solution array

make_real_space: converts a Fourier solution into real space for plotting and error calculations

make_plots: creates an animation comparing several solutions in real space

single_plot: creates an animation of a single solution in real space

generate_deriv_data_4func: computes data about a potential fourth order complete memory approximation from an exact simulation (for fitting optimal coefficients)

generate_deriv_data_4func2: computes data about a potential second order complete memory approximation from an exact simulation (for fitting optimal coefficients)

no_k_dependence_coeffs: computes optimal renormalization coefficients for a complete memory approximation ROM through fourth order given data from an “exact” simulation

no_k_dependence_coeffs2: computes optimal renormalization coefficients for a complete memory approximation ROM through second order given data from an “exact” simulation

find_scaling_law: takes optimal coefficients for a variety of resolutions N and degrees of dispersion epsilon and finds scaling laws to describe them

old_files: out-of-date files used in past iterations of the repository



%%%%%%%%%
%Figures% 
%%%%%%%%%
Files associated with figures generated for paper submissions

energy_flow: compares the net energy derivative in the first N = 20 modes according to a full model and the net energy derivatives associated with the second and fourth order terms of the complete memory approximation
produces: compare_energy.png

energy_presentation: creates plots for a presentation. Notably, presents the in and out flow of energy in the exact solution (compared to the Markov model). Compares Markov model, exact solution, unrenormalized 4th order ROM, and renormalized 4th order ROMs for N = 20
produces: blowup.png, kdv_energy.png, markov_energy.png

energy20: compares energy inflow and outflow for N = 20 modes according to exact, Markov, 2nd order renormalized complete memory approximation, 4th order renormalized complete memory approximation, and 4th order unrenormalized complete memory approximation
produces: energy_evo20.png

extrapolation_case: extrapolates our results to epsilon = 0.01, N = 256 4th order ROM (error in real space compared to Markov model)
produces: extrap_energy.png, extrap_error.png

log_err_presentation: produces plots indicating the error at t = 100 for Markov and 4th order models at a variety of resolutions (similar to long_time_error_PNAS)
produces: logerr.png

long_time_error_PNAS_SI: produces plots indicating the error at t = 100 for Markov, 2nd order, and 4th order models at a variety of resolutions
produces: logerr_SI.png

long_time_error_PNAS: produces plots indicating the error at t = 100 for Markov and 4th order models at a variety of resolutions
produces: logerr.png (in color)

long_time_error20: produces plots of the error for N = 20 Markov models, and second and fourth order ROMS for a long simulation
produces: real_err20.png

non_renormalized_blow_up: demonstrates that the unrenormalized models are unstable
produces: blowup.png

random_initial_conditions_2: computes optimal renormalization coefficients with random initial conditions of the form a*sin(x) + (1-a)*sin(2x) and plots the results
produces: prefactor2.png, N_exp2.png, eps_exp2.png

random_initial_conditions_3: computes optimal renormalization coefficients with random initial conditions of the form a*sin(x) + b*sin(2x) + (1-a-b)*sin(3x)
produces: prefactor3.png, N_exp3.png, eps_exp3.png

real_space_plots: displays examples of real space solutions according to the exact model, Markov, and renormalized fourth order complete memory approximation
produces: real_space20.png

ROM_sim_plots: simulates N = 20 model exactly, with Markov, and with renormalized 4th order model, then plots snapshots of the real space solution in order to produce an animation
produces: error.png, ROM_anim#.png

some_scaling_laws: computes optimal renormalization coefficients and the scaling laws describing them, then produces appealing log-log plots demonstrating the fit
produces: t2_eps.png, t4_eps.png, t2_N.png, t4_N.png

scaling_laws_presentation_plots: computes optimal renormalization coefficients and the scaling laws describing them, then produces appealing log-log plots demonstrating the fit for a presentation
produces: t_eps_pres.png, t_N_pres.png

some_scaling_laws_paper: computes optimal renormalization coefficients for fourth order complete memory approximations of KdV in using non dimensional quantities and plots them in an appealing way
produces: t_eps_N.png

some_scaling_laws_paper2: computes optimal renormalization coefficients for second order complete memory approximations of KdV in using non dimensional quantities and plots them in an appealing way
produces: t_eps_N2.png

non_ROM_anims: simulates N = 12 model exactly and with Markov and plots snapshots of the solutions to demonstrate the breakdown of accuracy for the Markov model.
produces: markov_anim#.png, exact_anim#.png

energy_deriv_pres: computes the exact energy derivative for N = 20 and compares it to the t^2 and t^4 terms for the presentation.
produces: energy_deriv.png

locked_in: computes optimal renormalization coefficients fit over different time ranges to show the “locking in” of the t^2 and t^4 model coefficients
produces: coeffs_convergence.png



%%%%%%%%%%%%%
%Simulations% 
%%%%%%%%%%%%%
Subdirectories associated with numerical experiments (all out of date now)