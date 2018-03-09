function plot_results(N)

N_full = N;
if exist(sprintf('full_u%i.mat',N_full),'file') == 2
    N_full = N_full+2;
end

load(sprintf('full_u%i.mat',N_full-2));
load(sprintf('full_t%i.mat',N_full-2));