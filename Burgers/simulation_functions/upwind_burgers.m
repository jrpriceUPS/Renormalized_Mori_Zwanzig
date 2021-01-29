function [t_list,u_list,u_list_real] = upwind_burgers(alpha,num_points,endtime,dt,howoften,epsilon)
%
%Computes solution to Burgers' equation using the upwind method
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  alpha       =  coefficient on nonlinear term
%
%  num_points  =  spatial resolution (number)
%
%  endtime     =  final time
%
%  dt          =  time step
%
%  howoften    =  how often to save data (set to zero to save all)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  t_list  =  vector of times
%
%  u_list  =  array of fourier transforms of exact solution
%
%  u_list_real = real space exact solution

%create spatial discretization
x = linspace(0,2*pi*(num_points-1)/num_points,num_points);

dx = x(2)-x(1);

%create time vector
t = 0:dt:endtime;
t_list = 0:dt*howoften:endtime;

%construct output array
u_list = zeros(length(x)/2,length(t_list));
u_list_real = zeros(length(x),length(t_list));

%save initial condition
u = epsilon*sin(x).';
index = 1;

%take FFT of initial condition for output
u_fft = fft_norm(u);
u_list(:,1) = u_fft(1:num_points/2);

%use upwind method to find solution at all later times
for i = 1:length(t)-1
    
    if i == 1
        % Check CFL condition
        dt_dx = abs(dt/dx);
        
        if dt_dx > 1
            sprintf('Failed CFL Test!')
            break;
        end
    end
    
    %find direction of gradient
    a_plus = max([alpha*u zeros(length(x),1)],[],2);
    a_minus = min([alpha*u zeros(length(x),1)],[],2);
    
    u_minus = zeros(length(x),1);
    u_plus = zeros(length(x),1);
    
    u_minus(3:end) = (3*u(3:end)-4*u(2:end-1)+u(1:end-2))/(2*dx);
    u_minus(2) = (3*u(2)-4*u(1)+u(end))/(2*dx);
    u_minus(1) = (3*u(1)-4*u(end)+u(end-1))/(2*dx);
    
    u_plus(1:end-2) = (-u(3:end)+4*u(2:end-1)-3*u(1:end-2))/(2*dx);
    u_plus(end-1) = (-u(1)+4*u(end)-3*u(end-1))/(2*dx);
    u_plus(end) = (-u(2)+4*u(1)-3*u(end))/(2*dx);
    
    %advance each entry
    u = u - dt*(a_plus.*u_minus + a_minus.*u_plus);
    

    %save results periodically
    if mod(i,howoften) == 0
        index = index + 1;
        u_fft = fft_norm(u);
        u_list(:,index) = u_fft(1:num_points/2);
        u_list_real(:,index) = u;
    end
    
    t(i)

end