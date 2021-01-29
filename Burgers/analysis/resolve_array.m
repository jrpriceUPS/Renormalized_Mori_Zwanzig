function [tsnaps,tmodel_size_list,tmodel_size_list_full] = resolve_array(u_list,t,alpha)
% 
% A function to take an output array of a full simulation and discard the
% point beyond which it can be taken to be unresolved.
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%    u  =  N length(t) solution array
%
%    t  =  list of times associated with the solution
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  tsnaps                 =  listing of index values for which the
%                            tolerance condition in satisfied
%
%  tmodel_size_list       =  a list of the magnitude of the t-model energy
%                            derivatives at each time (to be used for slicing
%                            arrays to different tolerances)
%
%  tmodel_size_list_full  =  a list of the magnitude of the t-model energy
%                            derivatives at each time of the full model

% initialize loop variables
tmodel_size_list = zeros(length(t),1);
tmodel_size_list_full = zeros(length(t),1);

% compute the size of the array and extract the number of resolved modes as
% well as the magnitude of M needed for dealiasing
s = size(u_list);
N = s(1)/2;
M = 3*N;
N_f = s(1);
M_f = 3*N_f;

% construct index lists
% Wave numbers updated
F_modes = [1:N,2*N+1:4*N+1,5*N+1:6*N];
G_modes = [N+1:5*N+1];
F_modes_f = [1:N_f,2*N_f+1:4*N_f+1,5*N_f+1:6*N_f];
G_modes_f = [N_f+1:5*N_f+1];

for i = 1:length(t)
   
    % display and save the current time
    time = t(i)
    u = u_list(1:N,i);
    u_f = u_list(1:N_f,i);
    
    %compute the model terms
    [~,~,t0tilde,u_full] = markov_term_Burgers(u,M,N,alpha,F_modes,G_modes);
    [t1,~,~] = tmodel_term_Burgers(u_full,t0tilde,alpha,F_modes,G_modes);
    t_model_size = time*sum(t1(1:N).*conj(u(:)) + conj(t1(1:N)).*u(:));
    tmodel_size_list(i) = abs(t_model_size);
    
    %compute the t-model term for the full model
    [~,~,t0tilde_f,u_full_f] = markov_term_Burgers(u_f,M_f,N_f,alpha,F_modes_f,G_modes_f);
    [t1_f,~,~] = tmodel_term_Burgers(u_full_f,t0tilde_f,alpha,F_modes_f,G_modes_f);
    t_model_size_full = time*sum(t1_f(1:N_f).*conj(u_f(:)) + conj(t1_f(1:N_f)).*u_f(:));
    tmodel_size_list_full(i) = abs(t_model_size_full);
    
end

max_tol = 1e-10;

tsnaps = find(tmodel_size_list < max_tol);
