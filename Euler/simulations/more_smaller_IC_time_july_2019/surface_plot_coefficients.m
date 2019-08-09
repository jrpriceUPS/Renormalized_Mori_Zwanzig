clear all;close all;

bigN = 36;
N = 4;

N_list = 4:2:bigN;
alpha_list = [1;0.1;0.01];
scaling = [1;0.75;0.5;0.25;0];

colors = {'b','r'};
whichN = [4,18];

count = 1;

for l = 1:length(whichN)
    N = whichN(l);
    N_ind = find(N_list == N);
    for j = 1:length(alpha_list)
        
        alpha = alpha_list(j);
        for i = 1:length(scaling)
            
            load(['coeffs_eps_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(i)-1) 'k_' num2str(bigN) '.mat'])
            
            figure(1)
            h = scatter3(scaling(i)-1,1/alpha_list(j),log(squeeze(coefficients(1,N_ind,1))),colors{l});
            hold on
            
            ylabel('1/epsilon','fontsize',16)
            xlabel('Time dependence t^{(alpha*t)}','fontsize',16)
            title('t-model coefficients for size','fontsize',16)
            
            if i == 1 && j == 1
                leg{count} = sprintf('ROM N = %i',N);
                count = count+1;
            else
                set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
        end
    end
end

legend(leg{:})