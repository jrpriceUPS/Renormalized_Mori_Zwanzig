This folder contains Matlab files to derive and simulate reduced order models for the Euler’s equations:

*FILL IN HERE*.

%%%%%%%%%%%%%%%%%%%%%%
%simulation_functions% 
%%%%%%%%%%%%%%%%%%%%%%
Files necessary for running an Euler ROM simulation

u_squishify: compresses a full-sized state array into a compressed form that does not include duplicate data from conjugates

u_fullify: reconstructs a full-sized state array from a compressed form

taylor_green: creates an array corresponding to the Taylor-Green vortex initial condition

fftn_norm: computes 3-dimensional FFT of array with 3x1 vector at each spatial point using the standard normalization

ifftn_norm: computes 3-dimensional IFFT of array with 3x1 vector at each spatial point using the standard normalization

fftn_norm_conv: computes 3-dimensional FFT of array with 3x3 matrices at each spatial point



%%%%%%%%%%%
%nonlinear% 
%%%%%%%%%%%
Files associated with the RHS of Euler full and ROM simulations

Ck: computes the convolution of two arrays by transforming to real space and back

convolve: computes the outer product in real space of two arrays (used in Ck convolution)

Ck_fill: auxiliary function used in filling in the correct entries of Ck

mode_clearer: auxiliary function that zeros out resolved or unresolved modes (for making resolved and unresolved versions of convolutions)

Dk: computes the (u,v) convolution and the (v,u) convolution and sums them - a common occurrence in Euler’s equation complete memory approximation terms

scaling_law: contains optimal algebraically-decaying renormalization coefficient laws as a function of number of included degrees and resolution

scaling_law_time: contains optimal constant renormalization coefficient laws as a function of number of included degrees and resolution

markov_term: computes the Markov term for Euler’s equations

tmodel_term: computes the first order complete memory approximation for Euler’s equations

t2model_term: computes the second order complete memory approximation for Euler’s equations

t3model_term: computes the third order complete memory approximation for Euler’s equations

t4model_term: computes the fourth order complete memory approximation for Euler’s equations

RHS: a function to do pre-processing and then call the requested RHS function, used as input to ode45

full_RHS: RHS of full simulation of Euler’s equations

markov_RHS: RHS of Markov approximation of Euler’s equations

tmodel_RHS: RHS of first order complete memory approximation of Euler’s equations

t2model_RHS: RHS of second order complete memory approximation of Euler’s equations

t3model_RHS: RHS of third order complete memory approximation of Euler’s equations

t4model_RHS: RHS of fourth order complete memory approximation of Euler’s equations




%%%%%%%%%%
%analysis% 
%%%%%%%%%%
Files associated with analysis of Euler results

linspecer: a function to produce visibly differentiated colors for plotting

create_data: simulates the full Euler’s to create reference data for renormalization

create_data48: specialized script for simulating and analyzing full Euler’s equations of resolution N = 48 (the default version produces file sizes that are too large to save)

resolve_array: identifies timesteps from a full simulation that we trust for fitting purposes

renormalize: computes optimal renormalization coefficients from a reference data set for a number of different resolutions

get_3D_energy: computes the energy in a subset of modes for a simulation at every timestep

helicity: computes the helicity of a simulation at every timestep

vorticity: computes the maximum of the vorticity of a simulation at every timestep

enstrophy: computes the enstrophy (2-norm of vorticity) of a simulation at every timestep



%%%%%%%%%
%Figures% 
%%%%%%%%%
Files associated with figures generated for the dissertation

burgers_style_instability: demonstrates the instability of constant renormalization coefficient complete memory approximation ROMs
produces: unstable_burgers_12.png

perturbative_figure: produces plots demonstrating the convergence of the complete memory approximation as the number of terms included increases
produces: perturbative_euler.png

coeff_plot: load and plot the scaling law behavior of optimal renormalization coefficients
produces: coeff_plot48_ROM.png, coeff_plot48_ROM_t.png

renormalized_multiple_res: produces energy, enstrophy, and maximum vorticity plots for algebraically decaying renormalization coefficient complete memory approximation ROMs
produces: energy_mult_24.png, energy_mult_slopes_24.png, enstrophy_mult_24.png, enstropy_mult_trim24.png, vorticity_mult_24.png, vorticity_mult_trim24.png, vorticity_int_24.png, vorticity_int_trim24.png

singularity_plots: produces the plots from renormalized_multiple_res from saved data (in case of crashed code)
produces: turn_times.png, slopes.png, slopes2.png, enstrophy.png, vorticity.png

vort_integral: computes the integral of the vorticity at all times in a simulation

turn_time_plot: produces a plot to help find what the time of energy drain converges to
produces: turn_time_plot.png

vort_int_short: produces a shortened plot of the vorticity integral
produces: vorticity_int_24.png

turn_times: time at which ROMs started draining energy

slopes: initial slope (log-log plot) at which ROMs started draining energy

slopes2: secondary slope (log-log plot) at which ROMs eventually drain energy

ens_max: maximum enstrophy for different ROMs

ens_max_time: time at which those maxima were achieved

vort_max: maximum vorticity for different ROMs

vort_max_time: time at which this maxima were achieved



%%%%%%%%%%%%%
%Simulations% 
%%%%%%%%%%%%%
Subdirectories associated with numerical experiments (all out of date now)