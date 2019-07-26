function [  ] = plot_vfield(u,v,max_arrows,arrow_len)
% 
% Adapted from quiverC2D
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%       x - 2D matrix, x components of initial points
%       y - 2D matrix, y components of initial points
%       u - 2D matrix, x components of arrows
%       v - 2D matrix, y components of arrows
%       maxNumArrows - a positive integer (non-integer should work as well)
%           number limiting the maximum number of plotted arrows. Since vectors
%           of length zero aren't plotted and the matrices are sampled
%           uniformly, it's possible that fewer arrows will be visible in the
%           end. If maxNumArrows is not given its value is set to 1000.

% uniform grid
N = length(u);
x = linspace(0,2*pi*(2*N-1)/(2*N),2*N).';
y = x;
    
 % 2D array of data points
[X,Y] = ndgrid(x,y);
x = X;
y = Y;

% Line width
lw = 1;
% Maximum of arrow head size
hs = 1;
% Colormap
colormap jet;

% initialization
if numel(u) > max_arrows
    N = ceil(sqrt(numel(u)/max_arrows));
    
    x = x(1:N:end,1:N:end);
    y = y(1:N:end,1:N:end);
    u = u(1:N:end,1:N:end);
    v = v(1:N:end,1:N:end);
end
s = size(u);

% taking care of possible issues
x = issues(x);
y = issues(y);
u = issues(u);
v = issues(v);
 if ~isequal(u,zeros(s)) 
     u = u./norm(u);
 end
 if ~isequal(v,zeros(s))
     v = v./norm(v);
 end
 
% colormap definition
I = sqrt(u.^2 + v.^2);
Ic = round( I/max(max(I))*64);
Ic( Ic == 0) = 1;
C = colormap;

% plotting
hold on;
for n = 1:s(1)
    for m = 1:s(2)
        if u(n,m) == 0 && v(n,m) == 0
        else
            quiver(x(n,m),y(n,m),arrow_len*u(n,m)./sqrt(u(n,m).^2+v(n,m).^2),arrow_len*v(n,m)./sqrt(u(n,m).^2+v(n,m).^2),'Color',C(Ic(n,m),:),'LineWidth',lw,'maxheadsize',hs);
            %quiver(x(n,m),y(n,m),u(n,m),v(n,m),'Color',C(Ic(n,m),:),'LineWidth',lw,'maxheadsize',hs);
        end
    end
end
hold off;

