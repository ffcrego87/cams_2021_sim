disp('Plotting')

% actual positions
figure
plot(Nb_min:Nb_max,Err_avg,'ro')
hold on
%er = errorbar(Nb_min:Nb_max,Err_avg,Err_avg-Err_low,Err_high-Err_avg);
er = errorbar(Nb_min:Nb_max,Err_avg,-Err_CI95_low,Err_CI95_high);
er.Color = 'r';                            
er.LineStyle = 'none';
for Nb=Nb_min:Nb_max
    plot(Nb*ones(n_runs,1),err_store(:,Nb-Nb_min+1),'k.')
end
hold off
setniceplot
xlim([Nb_min-1 Nb_max+1])
xticks(Nb_min:Nb_max)
ylabel('Estimation error RMS (m)')
xlabel('Number of bits')
if tr_loc_state
    print -dpdf nb_test_x
else
    print -dpdf nb_test_y
end