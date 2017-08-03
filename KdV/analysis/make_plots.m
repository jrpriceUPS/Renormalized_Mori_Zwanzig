function make_plots(t,x,u_real)
%
%Makes real space plot gif of several real space solutions
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  t       =  time vector
%
%  x       =  cell array of x-coordinates for each real space solution
%
%  u_real  =  cell array of each real space solution at each timestep (each
%             entry is length(x) x length(t))

colors = {'b','r','k','g','m'};
figure
for i = 1:length(t)-1
    for j = 1:length(u_real)
        plot(x{j},u_real{j}(:,i),colors{j})
        hold on
        axis([0,2*pi,-3,3])
        title(sprintf('t = %g', t(i+1)))
        drawnow
    end
    hold off
end