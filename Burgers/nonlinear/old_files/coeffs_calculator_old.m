function [coefficients,resid,t_norm] = coeffs_calculator(u,u_red,M,N,epsilon,alpha,F_modes,G_modes,k,t)
%
%Calculates the coefficients on the t-model, t^2-model, and t^3-model terms
%that will yield the best possible match of the rates of change of the
%resolved modes.

%compute full Markov term
[RHS,~] = markov_term(u,3/2*(M-1),M-1,alpha);

RHS_match = zeros(4,1);
nonlinear_terms = zeros(4,3);

N_list = (N-2*((1:4)-1)).';

for i = 1:4
    
    %compute reduced Markov term
    [nonlin0,u_full] = markov_term(u_red{i},2*N_list(i)+1,N_list(i),alpha);
    
    %compute t-model term
    [nonlin1,uu_star] = tmodel_term(u_full,nonlin0,alpha,F_modes{i});
    
    if i ==1
        t_norm = norm(nonlin1(1:N-2*(i-1)));
    end
    
    %compute t^2-model term
    [nonlin2,uk,uu,uk_uu_u,uk_uu_u_star] = t2model_term(u_full,nonlin0,uu_star,alpha,F_modes{i},G_modes{i},k{i},epsilon);
    
    %compute t^3-model term
    nonlin3 = t3model_term(alpha,F_modes{i},k{i},epsilon,u_full,uu,uu_star,uk,uk_uu_u,uk_uu_u_star);
    
    
    RHS_match(i) = sum(RHS(1:N_list(i)).*conj(u(1:N_list(i)))+conj(RHS(1:N_list(i))).*u(1:N_list(i)));
    RHS_match(i) = RHS_match(i) - sum(nonlin0(1:N_list(i)).*conj(u_red{i}) + conj(nonlin0(1:N_list(i))).*u_red{i});
    nonlinear_terms(i,1) = t*sum(nonlin1(1:N_list(i)).*conj(u_red{i}) + conj(nonlin1(1:N_list(i))).*u_red{i});
    nonlinear_terms(i,2) = -t^2/2*sum(nonlin2(1:N_list(i)).*conj(u_red{i}) + conj(nonlin2(1:N_list(i))).*u_red{i});
    nonlinear_terms(i,3) = t^3/6*sum(nonlin3(1:N_list(i)).*conj(u_red{i}) + conj(nonlin3(1:N_list(i))).*u_red{i}); 
end


F = @(x) -RHS_match + x(1)*(1./(N_list-x(2))).*nonlinear_terms(:,1) ...
                    + x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).*nonlinear_terms(:,2) ...
                    + x(1)*(1./(N_list-x(2))).*(1./(N_list-x(3))).*(1./(N_list-x(4))).*nonlinear_terms(:,3);
                
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
         
         
  if t > 0.25          
x = newton(F,DF,[1.5;0.5;-0.7;0.7],1e-8);
resid = norm(F(x));
coefficients = [x(1)*(1/(N-x(2)))
                x(1)*(1/(N-x(2)))*(1/(N-x(3)))
                x(1)*(1/(N-x(2)))*(1/(N-x(3)))*(1/(N-x(4)))];
            
  else
      resid = 0;
      coefficients = [0;0;0];
  end