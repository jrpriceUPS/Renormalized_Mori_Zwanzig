clear all;close all;

N = 2;
M = 2*N;

% make k array
kvec2 = [0:N-1,-N:1:-1];
[kx,ky] = ndgrid(kvec2,kvec2);
k = zeros(2*N,2*N,2);
k(:,:,1) = kx;
k(:,:,2) = ky;

%u = zeros(2*N,2*N,2);
%u(1,1,1) = 1;
%u(1,2,1) = 2;
%u(2,1,1) = 3;
%u(2,2,1) = rand*10 + 1i*rand*10;
%u(2,4,1) = rand*10 + 1i*rand*10;
% uniform grid
    x = linspace(0,2*pi*(2*M-1)/(2*M),2*M).';
    y = x;
    
    % 2D array of data points
    [X,Y] = ndgrid(x,y);
    
    % create initial condition
    eval = taylor_green(X,Y);
    u_full = fftn_norm(eval);
    newu = u_squishify(u_full,N);
%u(1,1,:) = 0;

%u = [diag(w);diag(w)]
%u(:,:,1) = [0 1i 2i 3i;4i 5i 6i 7i;8i 9i 10i 11i;12i 13i 14i 15i];
%u(:,:,2) = [0 1i 2i 3i;4i 5i 6i 7i;8i 9i 10i 11i;12i 13i 14i 15i];

%v = zeros(2*N,2*N,2);
%v(1,1,1) = 1;
%v(1,2,1) = 4;
%v(2,1,1) = 3;
%v(2,2,1) = 5;
%v(2,4,1) = rand*10 + 1i*rand*10;
%v(:,:,1) = randi(10,2*N,2*N);
%v(:,:,2) = randi(10,2*N,2*N);
%v = fftn_norm(v);
%v(:,:,1) = [0 -1i -2i -3i;-4i -5i -6i -7i;-8i -9i -10i -11i;-12i -13i -14i -15i];
%v(:,:,2) = [0 -1i -2i -3i;-4i -5i -6i -7i;-8i -9i -10i -11i;-12i -13i -14i -15i];
%v(1,1,:) = 0;
    
    % create initial condition
    eval = taylor_green(X,Y);
    v_full = fftn_norm(eval);
    newv = u_squishify(v_full,N);
    
%newu = u_squishify(u,N);
%newv = u_squishify(v,N);

u = u_fullify(newu,N);
v = u_fullify(newv,N);
u_big = u_fullify(newu,M);
v_big = u_fullify(newv,M);

C1 = zeros(size(u));
C2 = zeros(size(u));
convo1 = 0;
convo2 = 0;
count = 1;
count2 = 1;
zero_one = zeros(2,1);
zero_neg_one = zeros(2,1);
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
                                u_vec = squeeze(u(m,n,:));
                                v_vec = squeeze(v(p,q,:));
                                A1 = eye(2);
                                A2 = -(k_vec*k_vec')/(norm(k_vec)^2);
                                d1 = dot(k_vec,u_vec)*(A1*v_vec);
                                d2 = dot(k_vec,u_vec)*(A2*v_vec);
                                
                                if(k_vec(1)==0 && k_vec(2)==1)
                                    k_vec
                                    zero_one(:,count) = -1i*d1
                                    count = count+1;
                                end
                                
                                if(k_vec(1) == 0 && k_vec(2) == -1)
                                    k_vec
                                    zero_neg_one(:,count2) = -1i*d1
                                    count2 = count2+1;
                                end
                                
                                convo1 = convo1 + -1i*d1;
                                convo2 = convo2 + -1i*d2;
                            end
                        end
                    end
                end
            end
        end
        C1(i,j,:) = convo1;
        C2(i,j,:) = convo2;
        convo1 = zeros(2,1);
        convo2 = zeros(2,1);
    end
end
u_squishify(C1,N)
%u_squishify(C2,N)
% make k array
kvec2 = [0:M-1,-M:1:-1];
[kx,ky] = ndgrid(kvec2,kvec2);
k2 = zeros(2*M,2*M,2);
k2(:,:,1) = kx;
k2(:,:,2) = ky;

a = 2:M;
b = 2*M:-1:M+2;
a_tilde = N+1:M;
a_tilde2 = 2*N+1:M;
[C1_new,C2_new] = Ck(u_big,v_big,a,b,k2,a_tilde,a_tilde2);
u_squishify(C1_new,N);
u_squishify(C2,N) - u_squishify(C2_new,N);
max(max(max(max(abs(u_squishify(C1+C2,N)-u_squishify(C1_new+C2_new,N))))))
sum(abs(u(:)));