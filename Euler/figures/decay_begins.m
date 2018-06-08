function turn_times = decay_begins(N_list,percentage,end_time)

turn_times = zeros(length(N_list),1);

for i = 1:length(N_list)
    
    N = N_list(i);
    
    % plot the energy in some modes
    load(sprintf('energy_%i_%i.mat',N,end_time))
    load(sprintf('t4_%i_%i',N,end_time))
    
    
    % compute slopes and energy draining times
    turn_percent = percentage;
    
    energy_change = abs(energy-energy(1))/abs(energy(1));
    turn_times(i) = t4(find(energy_change>turn_percent,1));
    
    
end