function [times,energies,errors] = error_test(N,alpha,endtime)

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

if ~(exist(sprintf('u_list%i.mat',endtime),'file') == 2)
    
    [t_list,u_list] = upwind_burgers(1,10000,endtime,1e-4,1000);
    save(sprintf('u_list%i',endtime),'u_list');
    save(sprintf('t_list%i',endtime),'t_list');
    
else
    load(sprintf('u_list%i',endtime));
    load(sprintf('t_list%i',endtime));
end

leg = {'Exact','n = 1, constant coefficients','n = 2, constant coefficients','n = 3, constant coefficients','n = 4, constant coefficients','n = 1, decaying coefficients','n = 2, decaying coefficients','n = 3, decaying coefficients','n = 4, decaying coefficients'};




simulation_params.N = N;
simulation_params.alpha = alpha;
simulation_params.endtime = endtime;
simulation_params.print_time = 1;
simulation_params.time_range = t_list;

% Burgers-style 1 ROM
type = 'c1B';
simulation_params.time_dependence = 1;
simulation_params.degree = 4;
simulation_params.coeffs = scaling_law(N,type);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
simulation_params.time_dependence = 1;
[tc1B,uc1B] = PDE_solve(simulation_params);



% Burgers-style 1+2 ROM
type = 'c2B';
simulation_params.time_dependence = 1;
simulation_params.degree = 4;
simulation_params.coeffs = scaling_law(N,type);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
simulation_params.time_dependence = 1;
[tc2B,uc2B] = PDE_solve(simulation_params);


% Burgers-style 1+2+3 ROM
type = 'c3B';
simulation_params.time_dependence = 1;
simulation_params.degree = 4;
simulation_params.coeffs = scaling_law(N,type);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
simulation_params.time_dependence = 1;
[tc3B,uc3B] = PDE_solve(simulation_params);


% Burgers-style 1+2+3+4 ROM
type = 'c4B';
simulation_params.time_dependence = 1;
simulation_params.degree = 4;
simulation_params.coeffs = scaling_law(N,type);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
simulation_params.time_dependence = 1;
[tc4B,uc4B] = PDE_solve(simulation_params);




% KdV-style 1 ROM
type = 'c1KdV';
simulation_params.time_dependence = 1;
simulation_params.degree = 4;
simulation_params.coeffs = scaling_law(N,type);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
simulation_params.time_dependence = 0;
[tc1KdV,uc1KdV] = PDE_solve(simulation_params);


% KdV-style 1+2 ROM
type = 'c2KdV';
simulation_params.time_dependence = 1;
simulation_params.degree = 4;
simulation_params.coeffs = scaling_law(N,type);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
simulation_params.time_dependence = 0;
[tc2KdV,uc2KdV] = PDE_solve(simulation_params);


% KdV-style 1+2+3 ROM
type = 'c3KdV';
simulation_params.time_dependence = 1;
simulation_params.degree = 4;
simulation_params.coeffs = scaling_law(N,type);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
simulation_params.time_dependence = 0;
[tc3KdV,uc3KdV] = PDE_solve(simulation_params);


% KdV-style 1+2+3+4 ROM
type = 'c4KdV';
simulation_params.time_dependence = 1;
simulation_params.degree = 4;
simulation_params.coeffs = scaling_law(N,type);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
simulation_params.time_dependence = 0;
[tc4KdV,uc4KdV] = PDE_solve(simulation_params);







energy_exact = get_energy(u_list,N);
energyc1B = get_energy(uc1B,N);
energyc2B = get_energy(uc2B,N);
errc3B = get_energy(uc3B,N);
energyc4B = get_energy(uc4B,N);
errc1KdV = get_energy(uc1KdV,N);
errc2KdV = get_energy(uc2KdV,N);
errc3KdV = get_energy(uc3KdV,N);
erc4KdV = get_energy(uc4KdV,N);

figure(1)
hold off
plot(log(t_list),log(energy_exact),'linewidth',2)
hold on
plot(log(tc1B),log(energyc1B),'r')
plot(log(tc2B),log(energyc2B),'k')
plot(log(tc3B),log(errc3B),'c')
plot(log(tc4B),log(energyc4B),'m')

plot(log(tc1KdV),log(errc1KdV),'r--','linewidth',1.2)
plot(log(tc2KdV),log(errc2KdV),'k--','linewidth',1.2)
plot(log(tc3KdV),log(errc3KdV),'c--','linewidth',1.2)
plot(log(tc4KdV),log(erc4KdV),'m--','linewidth',1.2)
legend(leg{:},'location','southwest')

title(sprintf('N = %i',N),'fontsize',16)
xlabel('log(t)')
ylabel('log(energy)')
saveas(gcf,sprintf('Burgers_energy%i_%i',N,endtime),'png')
close

u_exact = u_list(1:N,:);
errc1B = sum((uc1B(:,1:length(tc1B)) - u_exact(:,1:length(tc1B))).*conj(uc1B(:,1:length(tc1B)) - u_exact(:,1:length(tc1B))),1)./sum(u_exact(:,1:length(tc1B)).*conj(u_exact(:,1:length(tc1B))),1);
errc2B = sum((uc2B(:,1:length(tc2B)) - u_exact(:,1:length(tc2B))).*conj(uc2B(:,1:length(tc2B)) - u_exact(:,1:length(tc2B))),1)./sum(u_exact(:,1:length(tc2B)).*conj(u_exact(:,1:length(tc2B))),1);
errc3B = sum((uc3B(:,1:length(tc3B)) - u_exact(:,1:length(tc3B))).*conj(uc3B(:,1:length(tc3B)) - u_exact(:,1:length(tc3B))),1)./sum(u_exact(:,1:length(tc3B)).*conj(u_exact(:,1:length(tc3B))),1);
errc4B = sum((uc4B(:,1:length(tc4B)) - u_exact(:,1:length(tc4B))).*conj(uc4B(:,1:length(tc4B)) - u_exact(:,1:length(tc4B))),1)./sum(u_exact(:,1:length(tc4B)).*conj(u_exact(:,1:length(tc4B))),1);

figure(1)
hold off
plot(tc1B,errc1B,'r','linewidth',1.5)
hold on
plot(tc2B,errc2B,'k','linewidth',1.5)
plot(tc3B,errc3B,'c','linewidth',1.5)
plot(tc4B,errc4B,'m','linewidth',1.5)
axis([0,endtime,0,2])
legend(leg{2:5},'location','northeast')

title(sprintf('N = %i',N),'fontsize',16)
xlabel('t')
ylabel('error')
saveas(gcf,sprintf('Burgers_err%i_%i',N,endtime),'png')
close


energies.c1B = energyc1B;
energies.c2B = energyc2B;
energies.c3B = energyc3B;
energies.c4B = energyc4B;
energies.c1KdV = energyc1KdV;
energies.c2KdV = energyc2KdV;
energies.c3KdV = energyc3KdV;
energies.c4KdV = energyc4KdV;

times.c1B = tc1B;
times.c2B = tc2B;
times.c3B = tc3B;
times.c4B = tc4B;
times.c1KdV = tc1KdV;
times.c2KdV = tc2KdV;
times.c3KdV = tc3KdV;
times.c4KdV = tc4KdV;

errors.c1B = errc1B;
errors.c2B = errc2B;
errors.c3B = errc3B;
errors.c4B = errc4B;

