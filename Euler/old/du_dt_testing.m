% an old script used to debug the full Euler algorithm


clear all;close all;

N = 4; % resolution
end_time = 1;

% uniform grid
x = linspace(0,2*pi*(2*N-1)/(2*N),2*N).';
y = x;
z = x;

% 3D array of data points
[X,Y,Z] = ndgrid(x,y,z);

eval = taylor_green(X,Y,Z);

u_full = fftn_norm(eval);

u = zeros(N,N,N,3,4);

a = 2:N;
b = 2*N:-1:N+2;

u(a,a,a,:,1) = u_full(a,a,a,:);
u(a,a,a,:,2) = u_full(a,a,b,:);
u(a,a,a,:,3) = u_full(a,b,a,:);
u(a,a,a,:,4) = u_full(b,a,a,:);

u(1,a,a,:,1) = u_full(1,a,a,:);
u(a,1,a,:,1) = u_full(a,1,a,:);
u(a,a,1,:,1) = u_full(a,a,1,:);

u(1,a,a,:,2) = u_full(1,a,b,:);
u(a,1,a,:,2) = u_full(a,1,b,:);

u(1,a,a,:,3) = u_full(1,b,a,:);
u(a,a,1,:,3) = u_full(a,b,1,:);

u(a,1,a,:,4) = u_full(b,1,a,:);
u(a,a,1,:,4) = u_full(b,a,1,:);

% make k array
k_vec = [0:N-1,-N:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k = zeros(2*N,2*N,2*N,3);
k(:,:,:,1) = kx;
k(:,:,:,2) = ky;
k(:,:,:,3) = kz;






u_full = u_fullify(u);

convo = convolve(u_full,u_full);
s = size(convo);

term1 = zeros(s(1:end-1));
term2 = zeros(s(1:end-1));

for i = 1:2*N
    for j = 1:2*N
        for l = 1:2*N
            
            if i == 1 && j == 1 && l == 1
                
                term1(i,j,l,:) = zeros(3,1);
                term2(i,j,l,:) = zeros(3,1);
                
            else
                
                wave_vec = squeeze(k(i,j,l,:));
                current_mat = squeeze(convo(i,j,l,:,:));
                
                term1(i,j,l,:) = -1j*current_mat*wave_vec;
                
                term2(i,j,l,:) = wave_vec*wave_vec.'*squeeze(term1(i,j,l,:))/norm(wave_vec)^2;
                
            end
        end
    end
end

m = 2*M-1;

nonzero_array = [3,3,3;
                 m,m,m;
                 3,3,1;
                 m,m,1;
                 3,1,3;
                 m,1,m;
                 1,3,3;
                 1,m,m;
                 3,1,1;
                 m,1,1;
                 1,3,1;
                 1,m,1;
                 1,1,3;
                 1,1,m;
                 1,1,1;
                 3,m,1;
                 m,3,1;
                 3,1,m;
                 m,1,3;
                 1,3,m;
                 1,m,3;
                 3,m,3;
                 m,3,m];
             
index_array = nonzero_array;
index_array(index_array == 3) = 2;
index_array(index_array == 1) = 0;
index_array(index_array == m) = -2;
             
% for i = 1:length(nonzero_array)
%     
%     term1(nonzero_array(i,1),nonzero_array(i,2),nonzero_array(i,3),:) = 0;
%     term2(nonzero_array(i,1),nonzero_array(i,2),nonzero_array(i,3),:) = 0;
%     
% end
% 
% max(max(max(max(abs(term1)))))
% max(max(max(max(abs(term2)))))


du_dt = u_squishify(term1,N) - u_squishify(term2,N);


new_du_dt = u_fullify(du_dt,M);

for i = 1:length(nonzero_array)
    
    index_array(i,:).'
    squeeze(new_du_dt(nonzero_array(i,1),nonzero_array(i,2),nonzero_array(i,3),:))
    pause
    
end
    

