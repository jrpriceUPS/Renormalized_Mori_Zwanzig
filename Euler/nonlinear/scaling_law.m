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

if degree == 1 % updated 3/23/2018 N=48
    coefficients = exp(0.346140204693076) * N^-1.100982456346863;
end

if degree == 2 % updated 3/23/2018 N=48
    coefficients(1) = exp(0.688731579483296) * N^-0.997313455364361;
    coefficients(2) = exp(1.115048630933493) * N^-2.126965489206877;
    end115048630933493

if degree == 3 % updated 3/23/2018 N=48
    coefficients(1) = exp(0.870581358602445) * N^-0.997205857074403;
    coefficients(2) = exp(1.638663408564890) * N^-2.114753283686792;
    coefficients(3) = exp(1.722885571602311) * N^-3.289326488342230;
end

if degree == 4 % updated 3/23/2018 N=48
    coefficients(1) = exp(1.052501710233077) * N^-1.043950293373851;
    coefficients(2) = exp(2.094977342459796) * N^-2.226203385845042;
    coefficients(3) = exp(2.734550827044432) * N^-3.538251864483979;
    coefficients(4) = exp(3.058650506815392) * N^-5.187724919959430;
end


