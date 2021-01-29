function du_dt = t2model_RHS_tau(params)
%
% Computes the t^2-model for every mode in 3D Euler
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%    u_full  =  full array of current Fourier state (2Mx2Mx2Mx3)
%
%         a  =  indices of positive resolved modes 1:N
%
%         b  =  indices of negative resolved modes -N:-1
%
%         k  =  array of wavenumbers (2Mx2Mx2Mx3)
%
%   a_tilde  =  indices of positive unresolved modes
%
%  a_tilde2  =  indices corresponding to modes included only for
%               dealiasing
%
%         N  =  maximal mode of reduced model
%
%      time  =  current time in simulation
%
%     coeff  =  array of renormalization coefficients
%
%       tau  =  degree of non-linearity of time dependence of
%              renormalization coefficients
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  du_dt  =  derivative of each mode

u_full = params.u_full;
a = params.a;
b = params.b;
k = params.k;
a_tilde = params.a_tilde;
a_tilde2 = params.a_tilde2;
N = params.N;
time = params.time;
coeff = params.coeff;
tau = params.tau;

% compute the full model term
[t0,t0hat,t0tilde] = markov_term(u_full,a,b,k,a_tilde,a_tilde2);

% compute the t-model term
[t1,~,t1tilde] = tmodel_term(u_full,t0tilde,a,b,k,a_tilde,a_tilde2);

% compute the t^2-model term
t2 = t2model_term(u_full,t0hat,t0tilde,t1tilde,a,b,k,a_tilde,a_tilde2);

t0 = u_squishify(t0,N);
t1 = u_squishify(t1,N);
t2 = u_squishify(t2,N);

% compute the derivative
du_dt = t0 + t1 * time^(1-1*tau) * coeff(1) + t2 * time^(2-2*tau) * coeff(2);