% old testing code used to check if the full Euler method worked correctly

%clear all;close all;
addpath ././simulation_functions
addpath ././nonlinear
addpath ././analysis

N = 4; % resolution
M = 3*N;

% make k array
k_vec = [0:M-1,-M:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k = zeros(2*M,2*M,2*M,3);
k(:,:,:,1) = kx;
k(:,:,:,2) = ky;
k(:,:,:,3) = kz;

u = 1i*u_squishify(k,N);

% % uniform grid
% x = linspace(0,2*pi*(2*M-1)/(2*M),2*M).';
% y = x;
% z = x;
% 
% % 3D array of data points
% [X,Y,Z] = ndgrid(x,y,z);
% 
% eval = taylor_green(X,Y,Z);
% 
% u_full = fftn_norm(eval);
% 
% u = u_squishify(u_full,N);

[du_dt,term1,term2] = full_euler(u,k,N,M);
du_dt = reshape(du_dt,N,N,N,3,4);
term1 = reshape(term1,N,N,N,3,4);
term2 = reshape(term2,N,N,N,3,4);



new_du_dt = u_fullify(du_dt,M);
new_term1 = u_fullify(term1,M);
new_term2 = u_fullify(term2,M);

du_dt_brute = du_dt_slow(u,k,N,M);
save du_dt_brute du_dt_brute
load du_dt_brute
new_du_dt_brute = u_fullify(du_dt_brute,M);

max(abs(new_du_dt_brute(:) - new_du_dt(:)))



% wavevector = squeeze(k(1,1,1,:))
% transform = squeeze(new_du_dt(1,1,1,:))
% brute =squeeze(new_du_dt_brute(1,1,1,:))
% pause

% for i = 2:N
%     for j = 2:N
%         for l = 2:N
%             wavevector = squeeze(k(i,j,l,:))
%             transform = squeeze(new_du_dt(i,j,l,:));
%             brute = squeeze(new_du_dt_brute(i,j,l,:));
%             diff = transform + brute
%             
%             wavevector = squeeze(k(2*M-i+2,2*M-j+2,2*M-l+2,:))
%             transform = squeeze(new_du_dt(2*M-i+2,2*M-j+2,2*M-l+2,:));
%             brute = squeeze(new_du_dt_brute(2*M-i+2,2*M-j+2,2*M-l+2,:));
%             diff = transform + brute
%             
%             pause
%         end
%     end
% end













% m = 2*M-1;
% 
% nonzero_array = [3,3,3;
%                  m,m,m;
%                  3,3,1;
%                  m,m,1;
%                  3,1,3;
%                  m,1,m;
%                  1,3,3;
%                  1,m,m;
%                  3,1,1;
%                  m,1,1;
%                  1,3,1;
%                  1,m,1;
%                  1,1,3;
%                  1,1,m;
%                  1,1,1;
%                  3,m,1;
%                  m,3,1;
%                  3,1,m;
%                  m,1,3;
%                  1,3,m;
%                  1,m,3;
%                  3,m,3;
%                  m,3,m];
%              
% index_array = nonzero_array;
% index_array(index_array == 3) = 2;
% index_array(index_array == 1) = 0;
% index_array(index_array == m) = -2;
% 
% for i = 1:length(nonzero_array)
%     
%     index_array(i,:).'
%     squeeze(new_du_dt_brute(nonzero_array(i,1),nonzero_array(i,2),nonzero_array(i,3),:))
%     pause
%     
% end