clear all;close all;

N = 8;
epsilon = 0.1;

x=linspace(0,2*pi*(2*N-1)/(2*N),2*N);
u_real = sin(x).';

U = sqrt(1/(2*pi)*integral(@(x) sin(x).^2,0,2*pi));
%U = 5;

x_scaling = epsilon/sqrt(U);
t_scaling = epsilon/(U^(3/2));

u = fft_norm(u_real);

x_prime = x / x_scaling;

u_real_prime = u_real/U;

u_prime = fft_norm(u_real_prime);