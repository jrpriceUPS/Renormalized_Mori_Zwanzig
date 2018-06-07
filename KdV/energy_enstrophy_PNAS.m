% creates energy and enstrophy plots for PNAS submission
clear all;close all;

N_list = 4:4:24;
endtime = 1000;

style{1} = '.';
style{2} = 'o';
style{3} = 'x';
style{4} = 's';
style{5} = '+';

[slopes,slopes2,turn_times,ens_max,ens_max_time,vort_max,vort_max_time] = renormalized_multiple_res_colorless(N_list,endtime,style,'eps')