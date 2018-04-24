function time_array = error_test(N,alpha,endtime)

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

leg = {'n = 1, constant coefficients','n = 2, constant coefficients','n = 3, constant coefficients','n = 4, constant coefficients','n = 1, decaying coefficients','n = 2, decaying coefficients','n = 3, decaying coefficients','n = 4, decaying coefficients'};


time_array = zeros(4,2);

simulation_params.N = N;
simulation_params.alpha = alpha;
simulation_params.endtime = endtime;
simulation_params.print_time = 1;
simulation_params.time_range = t_list;

% Burgers-style 1 ROM
if ~(exist(sprintf('uc1B_err%i.mat',endtime),'file') == 2)
    type = 'c1B';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 1;
    [tc1B,uc1B] = PDE_solve(simulation_params);
    save(sprintf('uc1B_err%i.mat',endtime),'uc1B')
    save(sprintf('tc1B_err%i.mat',endtime),'tc1B')
else
    load(sprintf('uc1B_err%i.mat',endtime))
    load(sprintf('tc1B_err%i.mat',endtime))
end
time_array(1,1) = min(diff(tc1B));


% Burgers-style 1+2 ROM
if ~(exist(sprintf('uc2B_err%i.mat',endtime),'file') == 2)
    type = 'c2B';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 1;
    [tc2B,uc2B] = PDE_solve(simulation_params);
    save(sprintf('uc2B_err%i.mat',endtime),'uc2B')
    save(sprintf('tc2B_err%i.mat',endtime),'tc2B')
else
    load(sprintf('uc2B_err%i.mat',endtime))
    load(sprintf('tc2B_err%i.mat',endtime))
end
time_array(2,1) = min(diff(tc2B));

% Burgers-style 1+2+3 ROM
if ~(exist(sprintf('uc3B_err%i.mat',endtime),'file') == 2)
    type = 'c3B';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 1;
    [tc3B,uc3B] = PDE_solve(simulation_params);
    save(sprintf('uc3B_err%i.mat',endtime),'uc3B')
    save(sprintf('tc3B_err%i.mat',endtime),'tc3B')
else
    load(sprintf('uc3B_err%i.mat',endtime))
    load(sprintf('tc3B_err%i.mat',endtime))
end
time_array(3,1) = min(diff(tc3B));

% Burgers-style 1+2+3+4 ROM
if ~(exist(sprintf('uc4B_err%i.mat',endtime),'file') == 2)
    type = 'c4B';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 1;
    [tc4B,uc4B] = PDE_solve(simulation_params);
    save(sprintf('uc4B_err%i.mat',endtime),'uc4B')
    save(sprintf('tc4B_err%i.mat',endtime),'tc4B')
else
    load(sprintf('uc4B_err%i.mat',endtime))
    load(sprintf('tc4B_err%i.mat',endtime))
end
time_array(4,1) = min(diff(tc4B));



% KdV-style 1 ROM
if ~(exist(sprintf('uc1KdV_err%i.mat',endtime),'file') == 2)
    type = 'c1KdV';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 0;
    [tc1KdV,uc1KdV] = PDE_solve(simulation_params);
    save(sprintf('uc1KdV_err%i.mat',endtime),'uc1KdV')
    save(sprintf('tc1KdV_err%i.mat',endtime),'tc1KdV')
else
    load(sprintf('uc1KdV_err%i.mat',endtime))
    load(sprintf('tc1KdV_err%i.mat',endtime))
end
time_array(1,2) = min(diff(tc1KdV));

% KdV-style 1+2 ROM
if ~(exist(sprintf('uc2KdV_err%i.mat',endtime),'file') == 2)
    type = 'c2KdV';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 0;
    [tc2KdV,uc2KdV] = PDE_solve(simulation_params);
    save(sprintf('uc2KdV_err%i.mat',endtime),'uc2KdV')
    save(sprintf('tc2KdV_err%i.mat',endtime),'tc2KdV')
else
    load(sprintf('uc2KdV_err%i.mat',endtime))
    load(sprintf('tc2KdV_err%i.mat',endtime))
end
time_array(2,2) = min(diff(tc2KdV));

% KdV-style 1+2+3 ROM
if ~(exist(sprintf('uc3KdV_err%i.mat',endtime),'file') == 2)
    type = 'c3KdV';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 0;
    [tc3KdV,uc3KdV] = PDE_solve(simulation_params);
    save(sprintf('uc3KdV_err%i.mat',endtime),'uc3KdV')
    save(sprintf('tc3KdV_err%i.mat',endtime),'tc3KdV')
else
    load(sprintf('uc3KdV_err%i.mat',endtime))
    load(sprintf('tc3KdV_err%i.mat',endtime))
end
time_array(3,2) = min(diff(tc3KdV));

% KdV-style 1+2+3+4 ROM
if ~(exist(sprintf('uc4KdV_err%i.mat',endtime),'file') == 2)
    type = 'c4KdV';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 0;
    [tc4KdV,uc4KdV] = PDE_solve(simulation_params);
    save(sprintf('uc4KdV_err%i.mat',endtime),'uc4KdV')
    save(sprintf('tc4KdV_err%i.mat',endtime),'tc4KdV')
else
    load(sprintf('uc4KdV_err%i.mat',endtime))
    load(sprintf('tc4KdV_err%i.mat',endtime))
end
time_array(4,2) = min(diff(tc4KdV));






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
errc1B = sum((uc1B - u_exact).*conj(uc1B - u_exact),1);
errc2B = sum((uc2B - u_exact).*conj(uc2B - u_exact),1);
errc3B = sum((uc2B - u_exact).*conj(uc2B - u_exact),1);
errc4B = sum((uc2B - u_exact).*conj(uc2B - u_exact),1);
errc1KdV = sum((uc1KdV - u_exact).*conj(uc1KdV - u_exact),1);
errc2KdV = sum((uc2KdV - u_exact).*conj(uc2KdV - u_exact),1);
errc3KdV = sum((uc3KdV - u_exact).*conj(uc3KdV - u_exact),1);
errc4KdV = sum((uc4KdV - u_exact).*conj(uc4KdV - u_exact),1);

figure(1)
hold off
plot(tc1B,errc1B,'r')
hold on
plot(tc2B,errc2B,'k')
plot(tc3B,errc3B,'c')
plot(tc4B,errc4B,'m')

plot(tc1KdV,errc1KdV,'r--','linewidth',1.2)
plot(tc2KdV,errc2KdV,'k--','linewidth',1.2)
plot(tc3KdV,errc3KdV,'c--','linewidth',1.2)
plot(tc4KdV,errc4KdV,'m--','linewidth',1.2)
legend(leg{:},'location','southwest')

title(sprintf('N = %i',N),'fontsize',16)
xlabel('t')
ylabel('error')
saveas(gcf,sprintf('Burgers_err%i_%i',N,endtime),'png')
close



