function coefficients = scaling_law(N,degree)
%
% Produces renormalization coefficients for a renormalized model of
% resolution N and degree "degree"
%
% Based upon simulations done on April 27 from N = 48 data with 
% algebraically decaying renormalization coefficients and a non-variable 
% window (N = 24 t-model between 1e-16 and 1e-10)
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

if degree == 1 % updated 4/27/2018 N=48
    coefficients = exp(0.464178051931995) * N^-1.076854649494771;
end

if degree == 2 % updated 4/27/2018 N=48
    coefficients(1) = exp(0.895358015405417) * N^-0.998464070418845;
    coefficients(2) = -exp(0.850572528271061) * N^-2.135828666672579;
end

if degree == 3 % updated 4/27/2018 N=48
    coefficients(1) = exp(0.974697906494128) * N^-0.962070195642592;
    coefficients(2) = -exp(1.129553143083177) * N^-2.041713874431020;
    coefficients(3) = exp(0.065345288512282) * N^-3.148049560058902;
end

if degree == 4 % updated 4/127/2018 N=48
    coefficients(1) = exp(1.134515716447412) * N^-0.923605415992834;
    coefficients(2) = -exp(1.610715982059490) * N^-1.958746049349437;
    coefficients(3) = exp(1.594061554873856) * N^-3.114959739896249;
    coefficients(4) = -exp(0.603079335580313) * N^-4.311867871313719;
end

% if degree == 1 % updated 4/13/2018 N=48
%     coefficients = exp(0.464178051931995) * N^-1.076854649494771;
% end
% 
% if degree == 2 % updated 4/13/2018 N=48
%     coefficients(1) = exp(0.895358015405417) * N^-0.998464070418845;
%     coefficients(2) = exp(1.543719708831006) * N^-2.135828666672578;
% end
% 
% if degree == 3 % updated 4/13/2018 N=48
%     coefficients(1) = exp(0.974697906498827) * N^-0.962070195645239;
%     coefficients(2) = exp(1.822700323654759) * N^-2.041713874437572;
%     coefficients(3) = exp(1.857104757767303) * N^-3.148049560074087;
% end
% 
% if degree == 4 % updated 4/13/2018 N=48
%     coefficients(1) = exp(1.135943524780088) * N^-0.985718175635504;
%     coefficients(2) = exp(2.222251934246004) * N^-2.089820199929696;
%     coefficients(3) = exp(2.826994626808396) * N^-3.289350267350072;
%     coefficients(4) = exp(2.868938248275644) * N^-4.623321597576260;
% end

% if degree == 1 % updated 3/23/2018 N=48
%     coefficients = exp(0.346140204693076) * N^-1.100982456346863;
% end
% 
% if degree == 2 % updated 3/23/2018 N=48
%     coefficients(1) = exp(0.688731579483296) * N^-0.997313455364361;
%     coefficients(2) = exp(1.115048630933493) * N^-2.126965489206877;
% end
% 
% if degree == 3 % updated 3/23/2018 N=48
%     coefficients(1) = exp(0.870581358602445) * N^-0.997205857074403;
%     coefficients(2) = exp(1.638663408564890) * N^-2.114753283686792;
%     coefficients(3) = exp(1.722885571602311) * N^-3.289326488342230;
% end
% 
% if degree == 4 % updated 3/23/2018 N=48
%     coefficients(1) = exp(1.052501710233077) * N^-1.043950293373851;
%     coefficients(2) = exp(2.094977342459796) * N^-2.226203385845042;
%     coefficients(3) = exp(2.734550827044432) * N^-3.538251864483979;
%     coefficients(4) = exp(3.058650506815392) * N^-5.187724919959430;
% end


