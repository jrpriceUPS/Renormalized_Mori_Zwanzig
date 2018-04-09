function [coeffs1,coeffs2,coeffs3,coeffs4] = log_periodic_scaling(N_max)

close all

N_list = 4:2:N_max/2;

for j = 1:4
    load(sprintf('c%i_%i',N_max,j))
end

% c1 results
coeffs1 = zeros(3,1);
for i = 1:1
    coeffs = c1;
    x0 = [-1.100982456346863;exp(0.346140204693076);0];
    coeffs1 = fminsearch(@(x) log_periodic_cost(x,coeffs,N_list),x0);
    coeffs1 = fminsearch(@(x) log_periodic_cost(x,coeffs,N_list),coeffs1);
    
    figure(1)
    subplot(2,2,1)
    plot(log(N_list),-log(coeffs) + log(coeffs1(2)*N_list.^coeffs1(1)),'b.','markersize',20)
    hold on
    plot(log(N_list),-log(coeffs) + log(coeffs1(2)*N_list.^coeffs1(1).*cos(coeffs1(3)*log(N_list))),'r.','markersize',8)
end

% c2 results
coeffs2 = zeros(3,2);
for i = 1:2
    coeffs = c2(i,:);
    
    
    if i == 1
        
        x0 = [ -0.997313455364361; exp(0.688731579483296); 0];
        
    elseif i == 2
        
        x0 = [ -2.126965489206877; exp(1.115048630933493); 0];
        
    end
    
    coeffs2(:,i) = fminsearch(@(x) log_periodic_cost(x,coeffs,N_list),x0);
    coeffs2(:,i) = fminsearch(@(x) log_periodic_cost(x,coeffs,N_list),coeffs2(:,i));
    
    figure(2)
    subplot(2,2,i)
    plot(log(N_list),-log(coeffs) + log(coeffs2(2,i)*N_list.^coeffs2(1,i)),'b.','markersize',20)
    hold on
    plot(log(N_list),-log(coeffs) + log(coeffs2(2,i)*N_list.^coeffs2(1,i).*cos(coeffs2(3,i)*log(N_list))),'r.','markersize',8)
end

% c3 results
coeffs3 = zeros(3,3);
for i = 1:3
    coeffs = c3(i,:);
    
    if i == 1
        
        x0 = [ -0.997205857074403; exp(0.870581358602445); 0];
        
    elseif i == 2
        
        x0 = [ -2.114753283686792; exp(1.638663408564890); 0];
        
    elseif i == 3
        
        x0 = [ -3.289326488342230; exp(1.722885571602311); 0];
        
    end
   
    
    coeffs3(:,i) = fminsearch(@(x) log_periodic_cost(x,coeffs,N_list),x0);
    coeffs3(:,i) = fminsearch(@(x) log_periodic_cost(x,coeffs,N_list),coeffs3(:,i));
    
    figure(3)
    subplot(2,2,i)
    plot(log(N_list),-log(coeffs) + log(coeffs3(2,i)*N_list.^coeffs3(1,i)),'b.','markersize',20)
    hold on
    plot(log(N_list),-log(coeffs) + log(coeffs3(2,i)*N_list.^coeffs3(1,i).*cos(coeffs3(3,i)*log(N_list))),'r.','markersize',8)
end

% c4 results
coeffs4 = zeros(3,4);
for i = 1:4
    coeffs = c4(i,:);
    
    if i == 1
        
        x0 = [ -1.043950293373851; exp(1.052501710233077); 0];
        
    elseif i == 2
        
        x0 = [ -2.226203385845042; exp(2.094977342459796); 0];
        
    elseif i == 3
        
        x0 = [ -3.538251864483979; exp(2.734550827044432); 0];
        
    elseif i == 4
        
        x0 = [ -5.187724919959430; exp(3.058650506815392); 0];
        
    end
   
    coeffs4(:,i) = fminsearch(@(x) log_periodic_cost(x,coeffs,N_list),x0);
    coeffs4(:,i) = fminsearch(@(x) log_periodic_cost(x,coeffs,N_list),coeffs4(:,i));
    
    figure(4)
    subplot(2,2,i)
    plot(log(N_list),-log(coeffs) + log(coeffs4(2,i)*N_list.^coeffs4(1,i)),'b.','markersize',20)
    hold on
    plot(log(N_list),-log(coeffs) + log(coeffs4(2,i)*N_list.^coeffs4(1,i).*cos(coeffs4(3,i)*log(N_list))),'r.','markersize',8)
end
