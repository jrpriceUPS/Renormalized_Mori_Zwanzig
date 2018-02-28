function t4 = t4model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Bhat,Btilde,Ehat,Etilde,Fhat,Ftilde,a,b,k,a_tilde)
%
% Computes the RHS for every mode in the t^4-model term for 3D Euler
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%   u_full  =  full array of current Fourier state (2Mx2Mx2Mx3)
%
%    t0hat  =  full array of current Fourier state of C_hat(u,u)
%
%  t0tilde  =  full array of current Fourier state of C_tilde(u,u)
%
%    t1hat  =  full array of current Fourier state of hat{t1-term}
%
%        a  =  indices of positive resolved modes 1:M
%
%        b  =  indices of negative resolved modes -M:-1
%
%        k  =  array of wavenumbers (2Mx2Mx2Mx3)
%
%  a_tilde  =  indices of unresolved modes
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%       t4  =  t^4-model term of derivative of each resolved mode
%
%    t4hat  =  resolved part of the t^4-model term of derivative of each resolved mode
%
%  t4tilde  =  unresolved part of the t^4-model term of derivative of each resolved mode

[term1aa,term1aa_hat,term1aa_tilde] = Dk(u_full,Ahat - 3*t1hat + 3*t1tilde - 5*Atilde,a,b,k,a_tilde);
[term1ab,term1ab_hat,term1ab_tilde] = Dk(u_full,5*t1hat - 3*Ahat + 3*Atilde - t1tilde,a,b,k,a_tilde);


[term1a,term1a_hat,term1a_tilde] = Dk(u_full,2*Fhat - 2*Ehat + 6*Bhat + 2*Etilde ...
                                           - 6*Ftilde - 2*Btilde + term1aa_hat ...
                                           + term1ab_tilde,a,b,k,a_tilde);
                                       
[term1b,term1b_hat,term1b_tilde] = Dk(t0hat,3*Ahat - 5*t1hat + t1tilde - 3*Atilde,a,b,k,a_tilde);
[term1c,term1c_hat,term1c_tilde] = Dk(t0tilde,3*t1hat - Ahat + 5*Atilde - 3*t1tilde,a,b,k,a_tilde);

[term1,term1_hat,term1_tilde] = Dk(u_full,term1a_tilde + term1b_tilde + term1c_tilde,a,b,k,a_tilde);

[term2a,term2a_hat,term2a_tilde] = Dk(u_full,2*t1hat - Ahat + 2*Atilde - t1tilde,a,b,k,a_tilde);
[term2,term2_hat,term2_tilde] = Dk(t0tilde,Etilde - 2*Ftilde - 2*Btilde + term2a_tilde,a,b,k,a_tilde);

[term3,~,~] = Dk(Atilde,Atilde-2*t1tilde,a,b,k,a_tilde);

[term4,~,~] = Dk(t1tilde,t1tilde,a,b,k,a_tilde);


t4 = 1/24*term1 + 1/6*term2 + 1/8*term3 + 1/8*term4;

P1_simpler = 1/8*term3 + 1/8*term4;








% testing other version!

[P11,~,~] = Ck(Ahat,Ahat,a,b,k,a_tilde);
[P12,~,~] = Ck(Ahat+t1hat,Ahat+t1hat,a,b,k,a_tilde);
[P13,~,~] = Ck(Ahat+Atilde,Ahat+Atilde,a,b,k,a_tilde);
[P14,~,~] = Ck(Ahat+t1hat+Atilde+t1tilde,Ahat+t1hat+Atilde+t1tilde,a,b,k,a_tilde);
[P15,~,~] = Dk(Ahat,3*Ahat + t1hat + 2*Atilde,a,b,k,a_tilde);
[P16,~,~] = Dk(Ahat+t1hat+Atilde+t1tilde,-6*Ahat + 6*t1hat,a,b,k,a_tilde);
[P17,~,~] = Dk(Ahat+Atilde,Atilde+t1tilde,a,b,k,a_tilde);

P1 = P11 + 1/4*P12 + P13 + 1/4*P14 - 1/2*P15 -1/24*P16 - 1/2*P17;


max(abs(P1_simpler(:) - P1(:)))
