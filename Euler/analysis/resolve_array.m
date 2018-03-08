function tmodel_size_list = resolve_array(u,t)
%
% A function to take an output array of a full simulation and discard the
% point beyond which it can be taken to be unresolved.
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
%  tmodel_size_list  =  a list of the magnitude of the t-model energy
%                       derivatives at each time (to be used for slicing
%                       arrays to different tolerances)

% compute the size of the array and extract the number of resolved modes as
% well as the magnitude of M needed for dealiasing
s = size(u);
N = s(1);
M = 2*N;

% construct index lists
a = 2:M;
b = 2*M:-1:M+2;
a_tilde = N+1:M;

% make k array
k_vec = [0:M-1,-M:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k = zeros(2*M,2*M,2*M,3);
k(:,:,:,1) = kx;
k(:,:,:,2) = ky;
k(:,:,:,3) = kz;

% initialize loop variables
tmodel_size_list = zeros(length(t),1);


% as long as the t-model energy derivative is less than tol, compute the
% t-model energy derivative of the next timestep
for i = 1:length(t)
    
    % extract the current u array
    current_u_temp = squeeze(u(:,:,:,:,:,i));
    
    % create the full version of the current u
    u_full = u_fullify(current_u_temp,M);
    
    % display and save the current time
    time = t(i)
    
    % compute the t-model term and reshape into the same shape as u
    [~,~,t0tilde] = markov_term(u_full,a,b,k,a_tilde);
    [~,t1hat,~] = tmodel_term(u_full,t0tilde,a,b,k,a_tilde);
    t1 = u_squishify(t1hat,N);
    
    % compute the energy derivative due to the t-model and record it
    t_model_size = time*sum(t1(:).*conj(current_u_temp(:))+conj(t1(:)).*current_u_temp(:))
    tmodel_size_list(i) = abs(t_model_size);
end