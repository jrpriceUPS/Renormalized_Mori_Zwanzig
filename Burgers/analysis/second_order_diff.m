function so_diff = second_order_diff(x,u)
% Second order central difference scheme for equally spaced 1-D data
% Uses second order forward and backward difference schemes for end points

so_diff = zeros(length(x),1);
h = x(2) - x(1);

for i = 1:length(x)
    if i == 1
        so_diff(i) = (-u(i+2)+4*u(i+1)-3*u(i))/(2*h);
    elseif i == length(x)
        so_diff(i) = (3*u(i)-4*u(i-1)+u(i-2))/(2*h);
    else
        so_diff(i) = (u(i+1)-u(i-1))/(2*h);
    end 
end