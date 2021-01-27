function [value,isterminal,direction] = energy_grow(~,u,init_energy,percent_increase,simulation_params)

N = simulation_params.N;
current_energy = get_energy(u,N);

value = double(current_energy < init_energy*(1+percent_increase));
if isnan(current_energy)
    
   value = 0;
   
end

isterminal = 1;
direction = 0;

end