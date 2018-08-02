function err = optimize(x,N_array,t_array,exact,R0,R1,R2,R3,R4)

err = 0;

for i = 1:length(N_array)
    
    resid = R0{i} - exact{i} ... 
            + x(1,1)*N_array{i}.^x(1,2).*t_array{i}.^(1-x(1,3)).*R1{i} ...
            + x(2,1)*N_array{i}.^x(2,2).*t_array{i}.^(2-x(2,3)).*R2{i} ...
            + x(3,1)*N_array{i}.^x(3,2).*t_array{i}.^(3-x(3,3)).*R3{i} ...
            + x(4,1)*N_array{i}.^x(4,2).*t_array{i}.^(4-x(4,3)).*R4{i};
        
    err = err + sum(resid(:).^2);
    
end

err
x(:,3)