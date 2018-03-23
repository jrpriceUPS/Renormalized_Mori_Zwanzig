function coefficients = scaling_law(N,degree);
%
% Produces renormalization coefficients for a renormalized model of
% resolution N and degree "degree"
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

if degree == 1 % updated 3/22/2018
    coefficients = exp(0.216042496107435) * N^-1.093111558750439;
end

if degree == 2 % updated 3/22/2018
    coefficients(1) = exp(0.807455983097315) * N^-1.095574434238564;
    coefficients(2) = exp(1.497163921272081) * N^-2.388816163431132;
end

if degree == 3 % updated 3/22/2018
    coefficients(1) = exp(0.770209345246030) * N^-0.996189505137105;
    coefficients(2) = exp(1.484542760836709) * N^-2.132733905239233;
    coefficients(3) = exp(1.438554055652650) * N^-3.280679641515915;
end

if degree == 4 % updated 3/22/2018
    coefficients(1) = exp(0.881697469455647) * N^-1.016137379094870;
    coefficients(2) = exp(1.810959536989056) * N^-2.199277103212345;
    coefficients(3) = exp(2.387434303594024) * N^-3.530254865250746;
    coefficients(4) = exp(4.114595061560149) * N^-5.992497764994379;
end