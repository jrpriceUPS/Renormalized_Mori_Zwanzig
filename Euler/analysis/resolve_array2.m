function [memory_size_list,tmodel_size_list_full] = resolve_array2(u,t)
%
% A function to take an output array of a full simulation and discard the
% point beyond which it can be taken to be unresolved. Uses full memory and
% not just t-model to make decision.
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%    u  =  NxNxNx3x4xlength(t) solution array
%
%    t  =  list of times associated with the solution
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  memory_size_list       =  a list of the magnitude of the memory energy
%                            derivatives at each time (to be used for slicing
%                            arrays to different tolerances)
%
%  tmodel_size_list_full  =  a list of the magnitude of the t-model energy
%                            derivatives at each time of the full model

% compute the size of the array and extract the number of resolved modes as
% well as the magnitude of M needed for dealiasing
s = size(u);
N = s(1)/2;
M = 3*N;
M_full = 6*N;

% construct index lists
a = 2:M;
b = 2*M:-1:M+2;
a_tilde = N+1:M;
a_tilde2 = 2*N+1:M;

% construct index lists
a_f = 2:M_full;
b_f = 2*M_full:-1:M_full+2;
a_tilde_f = 2*N+1:M_full;
a_tilde2_f = 4*N+1:M_full;

% make k array
k_vec = [0:M-1,-M:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k = zeros(2*M,2*M,2*M,3);
k(:,:,:,1) = kx;
k(:,:,:,2) = ky;
k(:,:,:,3) = kz;

% make k array
k_vec = [0:M_full-1,-M_full:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k_f = zeros(2*M_full,2*M_full,2*M_full,3);
k_f(:,:,:,1) = kx;
k_f(:,:,:,2) = ky;
k_f(:,:,:,3) = kz;

% initialize loop variables
memory_size_list = zeros(length(t),1);
tmodel_size_list_full = zeros(length(t),1);


% as long as the t-model energy derivative is less than tol, compute the
% t-model energy derivative of the next timestep
for i = 1:length(t)
    
    % extract the current u array
    temp_u = squeeze(u(:,:,:,:,:,i));
    temp_u_full = u_fullify(temp_u,M_full);
    u_current = u_squishify(temp_u_full,N);
    u_full = u_fullify(u_current,M);
    
    % display and save the current time
    time = t(i)
    
    % compute the t-model term and reshape into the same shape as u
    [~,t0hat,~] = markov_term(u_full,a,b,k,a_tilde,a_tilde2);
    t0 = u_squishify(t0hat,N);
    
    
    
    
    
    
    % compute the t-model term and reshape into the same shape as u
    [~,t0hat_f,t0tilde_f] = markov_term(temp_u_full,a_f,b_f,k_f,a_tilde_f,a_tilde2_f);
    [~,t1hat_f,~] = tmodel_term(temp_u_full,t0tilde_f,a_f,b_f,k_f,a_tilde_f,a_tilde2_f);
    t1_f = u_squishify(t1hat_f,2*N);
    memory = u_squishify(t0hat_f,N);
    
    difference = t0 - memory;
    
    % compute the energy derivative due to the t-model and record it
    memory_size = sum(difference(:).*conj(u_current(:))+conj(difference(:)).*u_current(:))
    memory_size_list(i) = abs(memory_size);
    t_model_size_full = time*sum(t1_f(:).*conj(temp_u(:))+conj(t1_f(:)).*temp_u(:))
    tmodel_size_list_full(i) = abs(t_model_size_full);
end