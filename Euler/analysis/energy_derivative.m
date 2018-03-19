function deriv = energy_derivative(u,t,params)


deriv = zeros(length(t),1);

for i = 1:length(t)
    
    u_temp = squeeze(u(:,:,:,:,:,i));
    
    du_dt = RHS(u_temp(:),t(i),params);
    deriv(i) = sum(u_temp(:).*conj(du_dt) + conj(u_temp(:)).*du_dt);
    
end