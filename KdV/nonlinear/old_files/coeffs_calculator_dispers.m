function coefficients = coeffs_calculator_dispers(RHS_match,nonlinear_terms,N_list)

F = @(x) -RHS_match ...
        + x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).*nonlinear_terms(:,1) ...
        + x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).*(1./(N_list-x(4))).*nonlinear_terms(:,2);

DF = @(x) [(1./(N_list-x(2))).*(1./(N_list-x(3))).*nonlinear_terms(:,1) ...
    + (1./(N_list-x(2))).*(1./(N_list-x(3))).*(1./(N_list-x(4))).*nonlinear_terms(:,2),...
    ...
    x(1)*(1./(N_list-x(2))).^2.*(1./(N_list-x(3))).*nonlinear_terms(:,1) ...
    + x(1)*(1./(N_list-x(2))).^2.*(1./(N_list-x(3))).*(1./(N_list-x(4))).*nonlinear_terms(:,2),...
    ...
    x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).^2.*nonlinear_terms(:,1) ...
    + x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).^2.*(1./(N_list-x(4))).*nonlinear_terms(:,2),...
    ...
    x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).*(1./(N_list-x(4))).^2.*nonlinear_terms(:,2)];


alpha_list = 0.1:0.3:2.8;
beta1_list = 0:0.1:0.9;
beta2_list = 0:0.2:1.8;
beta3_list = 0:0.2:1.8;

[a,b1,b2,b3] = ndgrid(alpha_list,beta1_list,beta2_list,beta3_list);
initial_conditions = [a(:) b1(:) b2(:) b3(:)];
[num_inits,~] = size(initial_conditions);


x = zeros(4,num_inits);
condition = zeros(num_inits,1);

for i = 1:num_inits
    
    [x(:,i),condition(i)] = newton(F,DF,initial_conditions(i,:).',1e-8);
    
end

[~,min_index] = min(condition);
x(:,min_index)

coefficients = [0
                x(1,min_index)*(1/(N_list(1)-x(2,min_index)))*(1/(N_list(1)-x(3,min_index)))
                x(1,min_index)*(1/(N_list(1)-x(2,min_index)))*(1/(N_list(1)-x(3,min_index)))*(1/(N_list(1)-x(4,min_index)))];
