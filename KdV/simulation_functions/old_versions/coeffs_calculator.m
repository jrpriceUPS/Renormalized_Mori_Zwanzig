function coefficients = coeffs_calculator(RHS_match,nonlinear_terms,N_list)
%
%Computes renormalization constants using method parallel to that developed
%in previous works by Panos. Requires "snapshots" of energy derivatives,
%against which we use a nonlinear solve to find the best coefficients to
%minimize the error. Only computes first three coefficients.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  RHS_match        =  vector of energy derivatives of resolved modes for several
%                      subsets of modes (measured when each "full" model becomes
%                      unresolved)
%
%  nonlinear_terms  =  matrix of the energy quantities due to each term in
%                      the memory approximation for different subsets of 
%                      modes (measured when each "full" model becomes unresolved)
%
%  N_list           =  list of resolutions used
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  coefficients  =  3x1 vector of optimal renormalization coefficients

%function handle of quantity we want to minimize
F = @(x) -RHS_match ...
        + x(1)*(1./(N_list-x(2))).*nonlinear_terms(:,1) ...
        + x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).*nonlinear_terms(:,2) ...
        + x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).*(1./(N_list-x(4))).*nonlinear_terms(:,3);

%gradient matrix of that function
DF = @(x) [(1./(N_list-x(2))).*nonlinear_terms(:,1) ...
    + (1./(N_list-x(2))).*(1./(N_list-x(3))).*nonlinear_terms(:,2) ...
    + (1./(N_list-x(2))).*(1./(N_list-x(3))).*(1./(N_list-x(4))).*nonlinear_terms(:,3),...
    ...
    x(1)*(1./(N_list-x(2))).^2.*nonlinear_terms(:,1) ...
    + x(1)*(1./(N_list-x(2))).^2.*(1./(N_list-x(3))).*nonlinear_terms(:,2) ...
    + x(1)*(1./(N_list-x(2))).^2.*(1./(N_list-x(3))).*(1./(N_list-x(4))).*nonlinear_terms(:,3),...
    ...
    x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).^2.*nonlinear_terms(:,2) ...
    + x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).^2.*(1./(N_list-x(4))).*nonlinear_terms(:,3),...
    ...
    x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).*(1./(N_list-x(4))).^2.*nonlinear_terms(:,3)];

%set of initial guesses
alpha_list = 0.1:0.3:2.8;
beta1_list = 0:0.1:0.9;
beta2_list = 0:0.2:1.8;
beta3_list = 0:0.2:1.8;

[a,b1,b2,b3] = ndgrid(alpha_list,beta1_list,beta2_list,beta3_list);
initial_conditions = [a(:) b1(:) b2(:) b3(:)];
[num_inits,~] = size(initial_conditions);


x = zeros(4,num_inits);
condition = zeros(num_inits,1);

%use newton's method to find each solution
for i = 1:num_inits
    
    [x(:,i),condition(i)] = newton(F,DF,initial_conditions(i,:).',1e-8);
    
end

%use the one with the smallest condition number
[~,min_index] = min(condition);
x(:,min_index)

%compute coefficients
coefficients = [x(1,min_index)*(1/(N_list(1)-x(2,min_index)))
                x(1,min_index)*(1/(N_list(1)-x(2,min_index)))*(1/(N_list(1)-x(3,min_index)))
                x(1,min_index)*(1/(N_list(1)-x(2,min_index)))*(1/(N_list(1)-x(3,min_index)))*(1/(N_list(1)-x(4,min_index)))];
