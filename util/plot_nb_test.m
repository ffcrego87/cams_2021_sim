disp('Plotting')

% actual positions
figure
bar(Nb_min:Nb_max,Err_avg)
hold on
er = errorbar(Nb_min:Nb_max,Err_avg,Err_avg-Err_low,Err_high-Err_avg);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
setniceplot
ylabel('Estimation error (m)')
xlabel('Number of bits')
if tr_loc_state
    print -dpdf nb_test_x
else
    print -dpdf nb_test_y
end