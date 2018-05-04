function deriv = energy_derivative(u,t,params)
%
% deriv = energy_derivative(u,t,params)
%
% Computes the energy derivative of every mode at all times for use in
% renormalization computations
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%       u  =  N x length(t) array of exact solution at all times
%
%       t  =  list of times associated with exact solution
%
%  params  =  parameters needed to compute the full Burgers' RHS
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  deriv  =  the exact energy derivative of every mode at all times


deriv = zeros(length(t),1);

for i = 1:length(t)
    
    u_temp = squeeze(u(:,:,:,:,:,i));
    
    du_dt = RHS(u_temp(:),t(i),params);
    deriv(i) = sum(u_temp(:).*conj(du_dt) + conj(u_temp(:)).*du_dt);
    
end