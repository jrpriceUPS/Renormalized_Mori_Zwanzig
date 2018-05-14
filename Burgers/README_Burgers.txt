This folder contains Matlab files to derive and simulate reduced order models for the Euler’s equations:

u_t + uu_x = 0.

%%%%%%%%%%%%%%%%%%%%%%
%simulation_functions% 
%%%%%%%%%%%%%%%%%%%%%%
Files necessary for running a Burgers ROM simulation

upwind_burgers: solves Burgers’ equation exactly using the upwind method

full_init_Burgers: takes simulation_params structure and updates it to include necessary parameters for running a full model of Burgers

ROM_init_Burgers: takes simulation_params structure and updates it to include necessary parameters for running a complete memory approximation ROM of Burgers

convolution_sum_Burgers: computes the nonlinear part of Burgers in real space and transforms back

fft_norm: computes the FFT of a vector using the standard normalization

ifft_norm: computes the IFFT of a vector using the standard normalization

PDE_solve: solves KdV with given parameters

RK4_stiff_nonstiff_step: computes a single step in an integration scheme known to function well for KdV with small dispersion



%%%%%%%%%%%
%nonlinear% 
%%%%%%%%%%%
Files associated with the RHS of Burgers full and ROM simulations

scaling_law: contains optimal renormalization coefficient laws (both constant and algebraically decaying) as a function of number of included degrees and resolution

renormalized_complete_Burgers: RHS of Burgers (computes only up through the degree included in the model)

markov_term_Burgers: computes the Markov term of Burgers

tmodel_term_Burgers: computes the first order complete memory approximation term for Burgers

t2model_term_Burgers: computes the second order complete memory approximation term for Burgers

t3model_term_Burgers: computes the third order complete memory approximation term for Burgers

t4model_term_Burgers: computes the fourth order complete memory approximation term for Burgers



%%%%%%%%%%
%analysis% 
%%%%%%%%%%
Files associated with analysis of Burgers results

linspecer: a function to produce visibly differentiated colors for plotting

get_energy: computes the energy in a subset of modes at all times from a solution array

make_real_space: converts a Fourier solution into real space for plotting and error calculations

make_plots: creates an animation comparing several solutions in real space

single_plot: creates an animation of a single solution in real space

energy_derivative: computes the exact energy derivative of every mode at every time from an exact trajectory

renormalize: computes optimal renormalization coefficients of both types for a variety of resolutions given an exact solution (uses all times)

renormalize_mult: computes optimal renormalization coefficients of both types for a variety of resolutions given an exact solution (uses a growing time window of [0,t] for all t to demonstrate “locking in” of coefficients as the window grows)

create_scaling_laws: computes the scaling laws and goodness of fit given optimal coefficient data



%%%%%%%%%
%Figures% 
%%%%%%%%%
Files associated with figures generated for the dissertation

create_data: creates exact solution data for fitting renormalization coefficients

scaling_law_plots: compute optimal constant renormalization coefficients for Burgers’ equation and plot the results on a log-log plot
produces: burgers_coeffs_constant.png

scaling_law_plots_KdV: compute optimal algebraically decaying renormalization coefficients for Burgers’ equation and plot the results on a log-log plot
produces: burgers_coeffs_algebraic.png

error_test: for a given resolution N, computes ROMs of all varieties and compares their errors and energy evolutions

generate_comparisons: executes error_test for a number of different resolutions N and plots and saves the result

dissertation_images: creates images using generate_comparisons for presentation in my dissertation
produces: energy_burgers.png, error_burgers.png

burgers_energy_plot:  produces exact plot of Burgers energy evolution for presentation
produces: burger_energy.png

coeff_test: compares different types of ROM solutions to Burgers’ equation (from testing phase, not final draft figures)


%%%%%%%%%%%%%
%Simulations% 
%%%%%%%%%%%%%
Subdirectories associated with numerical experiments (largely out of date)