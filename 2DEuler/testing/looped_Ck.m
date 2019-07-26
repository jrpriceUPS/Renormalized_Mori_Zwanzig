function [C,C_hat,C_tilde] = looped_Ck(N,M,u_full,v_full,a,b,k,a_tilde,a_tilde2)
%
% Computes the convolution of u and v using a nested loop (for testing
% purposes
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%    u_full  =  full Fourier form of first argument of C (2Mx2Mx2)
%
%    v_full  =  full Fourier form of second argument of C (2Mx2Mx2)
%
%         a  =  indices corresponding to positive modes 1:M
%
%         b  =  indices corresponding to negative modes -M:-1
%
%         k  =  array of wavevectors (2Mx2Mx2)
%
%   a_tilde  =  indices corresponding to positive unresolved
%
%  a_tilde2  =  indices corresponding to modes included only for
%               dealiasing
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%        C  =  Ck(v,w) array in compressed form (MxMx2x2)
%
%    C_hat  =  C_hat(v,w) array in compressed form (MxMx2x2)
%
%  C_tilde  =  C_tilde(v,w) array in compressed form (MxMx2x2)

C = zeros(size(u_full));
convo = 0;
for i = 1:(2*N)
    for j = 1:(2*N)
        if i == 1 && j == 1
        else
            k_vec = squeeze(k(i,j,:));
            for m = 1:(2*N)
                for n = 1:(2*N)
                    a_vec = squeeze(k(m,n,:));
                    for p = 1:(2*N)
                        for q = 1:(2*N)
                            b_vec = squeeze(k(p,q,:));
                            if (a_vec(1) + b_vec(1) == k_vec(1)) && (a_vec(2) + b_vec(2) == k_vec(2))
                                u_vec = squeeze(u_full(m,n,:));
                                v_vec = squeeze(v_full(p,q,:));
                                A = eye(2)-(k_vec*k_vec')/(norm(k_vec)^2);
                                d = dot(k_vec,u_vec)*(A*v_vec);
                                convo = convo + -1i*d;
                            end
                        end
                    end
                end
            end
        end
        C(i,j,:) = convo;
        convo = zeros(2,1);
    end
end

C_hat = mode_clearer(C,a_tilde);
C_tilde = C - C_hat;
C_tilde = mode_clearer(C_tilde,a_tilde2);

