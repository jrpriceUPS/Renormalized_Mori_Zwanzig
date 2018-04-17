function answer = find_scaling_law(coeff,epsilon,N_list)
%
%Given a set of coefficients arrayed according to N and epsilon, finds the
%scaling law of the form:
%
%     coeff = a * N^b * epsilon^c
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  coeff    =  length(N_list) x length(epsilon) array of coefficients as a
%              function of N and epsilon
%
%  N_list   =  corresponding list of N values
%
%  epsilon  =  corresponding list of epsilon values
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  answer  =  3 x 1 array containing [a,b,c]

num_eps = length(epsilon);
num_N = length(N_list);

A = zeros(num_eps*num_N,3);
A(:,1) = 1;

for i = 1:length(N_list)
   
    A((i-1)*length(epsilon)+1:i*length(epsilon),2) = log(N_list(i));
    A((i-1)*length(epsilon)+1:i*length(epsilon),3) = log(epsilon);
    
end


b = log(abs(coeff(:)));

answer = A \ b;

answer(1) = exp(answer(1));