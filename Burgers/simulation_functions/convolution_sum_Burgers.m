function convo=convolution_sum_Burgers(u,v,alpha)
%
%convo = convolution_sum_Burgers(u,v,alpha)
%
%computes the convolution of two vectors in frequency space C(u,v)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u, v  = frequency space vectors to be convolved
%
%  alpha  =  degree of nonlinearity
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  convo  = the convolution sum of u and v C(u,v)


%compute length of input vectors for use in generating wavenumber vector
L = length(u);

%compute wavenumber vector
k = [0:L/2-1,-L/2:-1].';

%compute the convolution sum by multiplying in real space, and returning to
%Fourier space
convo=fft_norm(ifft_norm(u).*ifft_norm(v));

%implement the derivative and the constant
convo = -alpha/2*1j*k.*convo;