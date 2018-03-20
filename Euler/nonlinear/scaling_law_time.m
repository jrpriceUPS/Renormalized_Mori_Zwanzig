function coefficients = scaling_law_time(N,degree)
%
% Produces renormalization coefficients for a renormalized model of
% resolution N and degree "degree" that does not cancel out the algebraic
% time dependence of the Taylor series
%
% Based upon simulations done on March 19 from N = 32 data with no time
% dependence and a non-variable windwo (N = 16 t-model between 1e-16 and
% 1e-10)
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%       N  =  resolution of ROM
%
%  degree  =  maximal degree of ROM term to include
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  coefficients  =  a degree x 1 array of the renormalization coefficients
%  to use

coefficients = zeros(degree,1);

if degree == 1
    coefficients = exp(0.491389062034627) * N^-1.156688850923880;
end

if degree == 2
    coefficients(1) = exp(0.420497563039746) * N^-0.905315378418248;
    coefficients(2) = exp(-0.036830897117974) * N^-1.704712914110429;
end

if degree == 3
    coefficients(1) = exp(0.645964779956798) * N^-0.909872993797799;
    coefficients(2) = exp(1.133169252422726) * N^-1.921277398375949;
    coefficients(3) = exp(1.508851623575983) * N^-3.215107947530718;
end

if degree == 4
    coefficients(1) = exp(0.861780641199717) * N^-0.976499344747932;
    coefficients(2) = exp(1.734115592499790) * N^-2.107802909318264;
    coefficients(3) = exp(2.474312865174682) * N^-3.480207961963250;
    coefficients(4) = exp(3.748006978495073) * N^-5.709839893517779;
end
