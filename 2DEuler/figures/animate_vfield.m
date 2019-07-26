function [	] = animate_vfield(uN,tN)
%UNTITLED 
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%
%  uN  =   array of Fourier modes from time zero to time end_time
%             (N x N x 2 x 2 x length(t))
%
%  tN  =  array of times associated with solution

N = length(uN);

for i = 1:length(tN)
    tN(i);
    u_real = ifftn_norm(u_fullify(uN(:,:,:,:,i),2*N));
    u = u_real(:,:,1);
    v = u_real(:,:,2);
    %plot_vfield(u,v,1000,.5);
    %pause(0.1)
end

