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

if degree == 1
    coefficients = exp(0.216042538948145) * N^-1.093111582148421;
end

if degree == 2
    coefficients(1) = exp(0.807464414266468) * N^-1.095579025575996;
    coefficients(2) = exp(1.497182389514764) * N^-2.388826221971421;
end

if degree == 3
    coefficients(1) = exp(0.770162592959785) * N^-0.996165830603617;
    coefficients(2) = exp(1.484422307729853) * N^-2.132672239666827;
    coefficients(3) = exp(1.438271947042133) * N^-3.280533408078200;
end

if degree == 4
    coefficients(1) = exp(0.881648270582790) * N^-1.016111918252141;
    coefficients(2) = exp(1.810852500536821) * N^-2.199218806267398;
    coefficients(3) = exp(2.387290431279766) * N^-3.530161767618130;
    coefficients(4) = exp(4.114478164325384) * N^-5.992383509914074;
end
