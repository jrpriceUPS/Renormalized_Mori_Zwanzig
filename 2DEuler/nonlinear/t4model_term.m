function t4 = t4model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Bhat,Btilde,Ehat,Etilde,Fhat,Ftilde,a,b,k,a_tilde,a_tilde2)
%
% Computes the RHS for every mode in the t^4-model term for 2D Euler
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%    u_full  =  full array of current Fourier state (2Mx2Mx2)
%
%     t0hat  =  full array of current Fourier state of C_hat(u,u)
%
%   t0tilde  =  full array of current Fourier state of C_tilde(u,u)
%
%     t1hat  =  full array of current Fourier state of hat{t1-term}
%
%   t1tilde  =  full array of current Fourier state of tilde{t1-term}
%
%      Ahat  =  resolved part of Dk(hat{u},hat{t0})
%
%    Atilde  =  unresolved part of Dk(hat{u},hat{t0})
%
%      Bhat  =  resolved part of Dk(hat{u},tilde{t1}-Atilde) 
%
%    Btilde  =  unresolved part of Dk(hat{u},tilde{t1}-Atilde) 
%
%      Ehat  =  resolved part of Dk(hat{t0},tilde{t0})
%
%    Etilde  =  unresolved part of Dk(hat{t0},tilde{t0})
%
%      Fhat  =  resolved part of Ck(hat{t0},hat{t0})
%
%    Ftilde  =  unresolved part of Ck(hat{t0},hat{t0})
%
%         a  =  indices of positive resolved modes 1:M
%
%         b  =  indices of negative resolved modes -M:-1
%
%         k  =  array of wavenumbers (2Mx2Mx2)
%
%   a_tilde  =  indices of unresolved modes
%
%  a_tilde2  =  indices corresponding to modes included only for
%               dealiasing
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%       t4  =  t^4-model term of derivative of each resolved mode

[~,term1aa_hat,~] = Dk(u_full,Ahat - 3*t1hat - 5*Atilde + 3*t1tilde,a,b,k,a_tilde,a_tilde2);
[~,~,term1ab_tilde] = Dk(u_full,-3*Ahat + 5*t1hat + 3*Atilde - t1tilde,a,b,k,a_tilde,a_tilde2);


[~,~,term1a_tilde] = Dk(u_full,2*Fhat - 2*Ehat + 6*Bhat + 2*Etilde ...
                                           - 6*Ftilde - 2*Btilde + term1aa_hat ...
                                           + term1ab_tilde,a,b,k,a_tilde,a_tilde2);
                                       
[~,~,term1b_tilde] = Dk(t0hat,3*Ahat - 5*t1hat - 3*Atilde + t1tilde,a,b,k,a_tilde,a_tilde2);
[~,~,term1c_tilde] = Dk(t0tilde,-Ahat + 3*t1hat + 5*Atilde - 3*t1tilde,a,b,k,a_tilde,a_tilde2);

[term1,~,~] = Dk(u_full,term1a_tilde + term1b_tilde + term1c_tilde,a,b,k,a_tilde,a_tilde2);

[~,~,term2a_tilde] = Dk(u_full,Ahat-2*t1hat-2*Atilde+t1tilde,a,b,k,a_tilde,a_tilde2);
[term2,~,~] = Dk(t0tilde,-Etilde + 2*Ftilde + 2*Btilde + term2a_tilde,a,b,k,a_tilde,a_tilde2);

[term3,~,~] = Dk(Atilde,Atilde-2*t1tilde,a,b,k,a_tilde,a_tilde2);

[term4,~,~] = Dk(t1tilde,t1tilde,a,b,k,a_tilde,a_tilde2);


t4 = term1 - 4*term2 - 3*term3 - 3*term4;

