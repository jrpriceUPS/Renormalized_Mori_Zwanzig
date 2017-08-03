function single_plot(t_list,u_list,N)
%
%Generates animation of solution for a single run
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  t_list  =  vector of time
%
%  u_list  =  array of Fourier modes at all time steps
%
%  N       =  number of modes to use

[x,u_real] = make_real_space(u_list,N);

make_plots(t_list,{x},{u_real});