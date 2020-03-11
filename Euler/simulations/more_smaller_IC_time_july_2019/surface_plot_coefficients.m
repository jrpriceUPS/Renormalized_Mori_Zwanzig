clear all;close all;

bigN = 36;
N = 4;

N_list = 4:2:bigN;
alpha_list = [1;0.1;0.01];
scaling = [1;0.75;0.5;0.25;0];

colors = {'b','r'};
whichN = 16;

count = 1;

for l = 1:length(whichN)
    N = whichN(l);
    N_ind = find(N_list == N);
    
    for i = 1:length(scaling)
        
        
        new_list = zeros(3,3);
        
        for j = 1:length(alpha_list)
            alpha = alpha_list(j);
            
            load(['coeffs_eps_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(i)-1) 'k_' num2str(bigN) '.mat'])
            load(['t4_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(i)-1) 'k_' num2str(N) '.mat'])
            
            new_list(j,1) = scaling(i)-1;
            new_list(j,2) = 1/alpha_list(j);
            new_list(j,3) = log(squeeze(coefficients(1,N_ind,1)));
            
            if t4(end) == 1000
                my_color = 'b';
            else
                my_color = 'r';
            end
            
            figure(1)
            h = scatter3(scaling(i)-1,1/alpha_list(j),log(squeeze(coefficients(1,N_ind,1))),my_color);
            hold on
            
            ylabel('1/epsilon','fontsize',16)
            xlabel('Time dependence t^{(alpha*k)}','fontsize',16)
            title('t-model coefficients for size','fontsize',16)
            
        end
        plot3(new_list(:,1),new_list(:,2),new_list(:,3),'k')
    end
end