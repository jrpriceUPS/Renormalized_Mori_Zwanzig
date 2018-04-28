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

if degree == 1 % updated 4/27/2018 N=48
    coefficients = exp(0.453286097514251) * N^-1.132132124518987;
end

if degree == 2 % updated 4/27/2018 N=48
    coefficients(1) = exp(0.517051066875360) * N^-0.925593175897768;
    coefficients(2) = -exp(-0.217225712380335) * N^-1.878787319537266;
end

if degree == 3 % updated 4/27/2018 N=48
    coefficients(1) = exp(0.670401224026936) * N^-0.915282161281184;
    coefficients(2) = -exp(0.451207931290370) * N^-1.924454108237361;
    coefficients(3) = exp(-0.556384942764219) * N^-3.121926846333570;
end

if degree == 4 % updated 4/27/2018 N=48
    coefficients(1) = exp(0.897702864150992) * N^-0.905755331418320;
    coefficients(2) = -exp(1.139284643357042) * N^-1.924193989726993;
    coefficients(3) = exp(0.966195600572123) * N^-3.091445690727842;
    coefficients(4) = -exp(-0.150995180299051) * N^-4.306096893189370;
end

% if degree == 1 % updated 4/14/2018 N=48
%     coefficients = exp(0.453286097514251) * N^-1.132132124518987;
% end
% 
% if degree == 2 % updated 4/14/2018 N=48
%     coefficients(1) = exp(0.517051066875360) * N^-0.925593175897768;
%     coefficients(2) = exp(0.475921468179610) * N^-1.878787319537266;
% end
% 
% if degree == 3 % updated 4/14/2018 N=48
%     coefficients(1) = exp(0.670401224027101) * N^-0.915282161281344;
%     coefficients(2) = exp(1.144355111850448) * N^-1.924454108237625;
%     coefficients(3) = exp(1.235374526465223) * N^-3.121926846334537;
% end
% 
% if degree == 4 % updated 4/14/2018 N=48
%     coefficients(1) = exp(0.845453094356771) * N^-0.946192807992691;
%     coefficients(2) = exp(1.614834267785675) * N^-2.001809193435187;
%     coefficients(3) = exp(2.122659874790261) * N^-3.233176429294022;
%     coefficients(4) = exp(1.762476498429390) * N^-4.484502875116744;
% end


% if degree == 1 % updated 3/24/2018 N=48
%     coefficients = exp(0.346140204693076) * N^-1.100982456346863;
% end
% 
% if degree == 2 % updated 3/24/2018 N=48
%     coefficients(1) = exp(0.450485290588850) * N^-0.910866820408359;
%     coefficients(2) = exp(0.241869704035654) * N^-1.812532078606838;
% end
% 
% if degree == 3 % updated 3/24/2018 N=48
%     coefficients(1) = exp(0.700418137617023) * N^-0.935138605878111;
%     coefficients(2) = exp(1.217147112806629) * N^-1.962577643260228;
%     coefficients(3) = exp(1.431162474047888) * N^-3.188011710973634;
% end
% 
% if degree == 4 % updated 3/24/2018 N=48
%     coefficients(1) = exp(0.898823727373543) * N^-0.988857926218814;
%     coefficients(2) = exp(1.753386619806672) * N^-2.103964030340424;
%     coefficients(3) = exp(2.388790403108812) * N^-3.415262749214611;
%     coefficients(4) = exp(2.412309528581516) * N^-4.945061502758477;
% end