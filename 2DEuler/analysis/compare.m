function compare(u,t,N,end_time)
    u_full = u;
    t_full = t;
    fig = 1;
    for i = 4:2:N
        create_Markov_data(u_full,i,end_time);
        load(sprintf('Mu%i.mat',i));
        load(sprintf('Mt%i.mat',i));
        figure(fig);
        E = get_2D_energy(u,i);
        plot(t,E);
        hold on
        Full = get_2D_energy(u_full,i);
        plot(t_full,Full);
        title(sprintf('Energy Evolution in Markov Model vs Exact, N = %i',i))
        ylabel('Energy in resolved modes')
        xlabel('Time')
        legend('Markov model','exact','Location','northwest')
        fig = fig + 1;
    end
end