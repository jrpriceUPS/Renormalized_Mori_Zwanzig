% produces plots demonstrating the locking in of Burgers coefficients


clear all;close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

endtime = 10;
N_list = 10:2:22;

% load the exact data if it exists (or create it)
if ~(exist(sprintf('u_list%i.mat',endtime),'file') == 2)
    
    [t_list,u_list,exact_derivative] = create_data(1,10000,10,1e-4,100);
     
else
     load(sprintf('u_list%i',endtime));
     load(sprintf('t_list%i',endtime));
     load(sprintf('exact_derivative%i',endtime));
end

t_is_2 = find(t_list == 2);
range = t_is_2:100:length(t_list);

coeffs1 = zeros(length(N_list),1,length(range));
coeffs2 = zeros(length(N_list),2,length(range));
coeffs3 = zeros(length(N_list),3,length(range));
coeffs4 = zeros(length(N_list),4,length(range));

coeffs1_c = coeffs1;
coeffs2_c = coeffs2;
coeffs3_c = coeffs3;
coeffs4_c = coeffs4;

for i = 1:length(range)
    
    window = 1:range(i);
    t = t_list(window);
    u = u_list(:,window);
    du = exact_derivative(:,window);
    
    [c1,c2,c3,c4] = renormalize(1,N_list,u,t,du,0);
    [c1_c,c2_c,c3_c,c4_c] = renormalize(1,N_list,u,t,du,1);
    
    coeffs1(:,:,i) = c1.';
    coeffs2(:,:,i) = c2.';
    coeffs3(:,:,i) = c3.';
    coeffs4(:,:,i) = c4.';
    
    coeffs1_c(:,:,i) = c1_c.';
    coeffs2_c(:,:,i) = c2_c.';
    coeffs3_c(:,:,i) = c3_c.';
    coeffs4_c(:,:,i) = c4_c.';
    
end


figure(1)
hold off
plot(t_list(range),squeeze(coeffs1(:,1,:)))
title('Only fit R1 model, decaying','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R1','png')




figure(2)
hold off
subplot(2,1,1)
plot(t_list(range),squeeze(coeffs2(:,1,:)))
title('R1+R2 model, decaying','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,1,2)
plot(t_list(range),squeeze(coeffs2(:,2,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R2','png')



figure(3)
hold off
subplot(3,1,1)
plot(t_list(range),squeeze(coeffs3(:,1,:)))
title('R1+R2+R3 model, decaying','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(3,1,2)
plot(t_list(range),squeeze(coeffs3(:,2,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(3,1,3)
plot(t_list(range),squeeze(coeffs3(:,3,:)))
ylabel('R3','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R3','png')



figure(4)
hold off
subplot(2,2,1)
plot(t_list(range),squeeze(coeffs4(:,1,:)))
title('R1+R2+R3+R4 model, decaying','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,2)
plot(t_list(range),squeeze(coeffs4(:,2,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,3)
plot(t_list(range),squeeze(coeffs4(:,3,:)))
ylabel('R3','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,4)
plot(t_list(range),squeeze(coeffs4(:,4,:)))
ylabel('R4','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R4','png')










figure(5)
hold off
plot(t_list(range),squeeze(coeffs1_c(:,1,:)))
title('Only fit R1 model, constant','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R1_c','png')




figure(6)
hold off
subplot(2,1,1)
plot(t_list(range),squeeze(coeffs2_c(:,1,:)))
title('R1+R2 model, constant','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,1,2)
plot(t_list(range),squeeze(coeffs2_c(:,2,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R2_c','png')



figure(7)
hold off
subplot(3,1,1)
plot(t_list(range),squeeze(coeffs3_c(:,1,:)))
title('R1+R2+R3 model, constant','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(3,1,2)
plot(t_list(range),squeeze(coeffs3_c(:,2,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(3,1,3)
plot(t_list(range),squeeze(coeffs3_c(:,3,:)))
ylabel('R3','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R3_c','png')



figure(8)
hold off
subplot(2,2,1)
plot(t_list(range),squeeze(coeffs4_c(:,1,:)))
title('R1+R2+R3+R4 model, constant','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,2)
plot(t_list(range),squeeze(coeffs4_c(:,2,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,3)
plot(t_list(range),squeeze(coeffs4_c(:,3,:)))
ylabel('R3','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,4)
plot(t_list(range),squeeze(coeffs4_c(:,4,:)))
ylabel('R4','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R4_c','png')