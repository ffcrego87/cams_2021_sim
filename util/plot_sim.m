disp('Plotting')

% distinguishable colors
colors=distinguishable_colors(NAgents);

% actual positions
figure
for i=1:NAgents
    plot(Agents_log{i}.x(1,:),Agents_log{i}.x(2,:),'Color',colors(i,:))
    hold on
end
for i=1:NAgents
    plot(Output_select{i}(1,:)*Agents_log{i}.x_hat,Output_select{i}(2,:)*Agents_log{i}.x_hat,'Color',colors(i,:),'LineStyle','--')
    hold on
end
for j=1:(floor(T_plot_form)/dt):ifinal
    for i=1:NAgents
        for k=1:NAgents
            plot([Agents_log{i}.x(1,j) Agents_log{k}.x(1,j)],[Agents_log{i}.x(2,j) Agents_log{k}.x(2,j)],'k--')
        end
        plot([Agents_log{i}.x(1,j) Output_select{i}(1,:)*Agents_log{i}.x_hat(:,j)],[Agents_log{i}.x(2,j) Output_select{i}(2,:)*Agents_log{i}.x_hat(:,j)],'k--')
    end
end
hold off
setniceplot
axis equal
xlabel('X (m)')
ylabel('Y (m)')
if tr_loc_state
    print -dpdf pos_x
else
    print -dpdf pos_y
end

%velocity
figure
for i=1:NAgents
    plot(Agents_log{i}.x(3,:),Agents_log{i}.x(4,:),'Color',colors(i,:))
    hold on
end
for i=1:NAgents
    plot(Output_select{i}(3,:)*Agents_log{i}.x_hat,Output_select{i}(4,:)*Agents_log{i}.x_hat,'Color',colors(i,:),'LineStyle','--','Marker','*')
    hold on
end
hold off
setniceplot
axis equal
xlabel('X')
ylabel('Y')
title('Velocity')
%print -dpdf

%Control
figure
for i=1:NAgents
    plot(POutput_select{i}(1,:)*Agents_log{i}.u,POutput_select{i}(2,:)*Agents_log{i}.u,'Color',colors(i,:))
    hold on
end
hold off
setniceplot
axis equal
xlabel('X (m)')
ylabel('Y (m)')
title('Control')
%print -dpdf

%leader speed
figure
plot(dt*(1:ifinal)-dt,Agents_log{1}.y(1,:),'b-');
hold on
plot(dt*(1:ifinal)-dt,Agents_log{1}.y(2,:),'r-');
hold off
setniceplot
xlim([0 dt*(ifinal-1)])
ylabel('velocity (m/s)')
xlabel('t (s)')
title('Leader speed')
%print -dpdf

%inter-vehicle distances
figure
for i=1:NAgents
    for j=(2:NAgents)-1+NBeacons
        plot(dt*(1:ifinal)-dt,Agents_log{i}.y(j+dimp*(i<=N_dvl),:),'Color',colors(i,:))
        hold on
        plot(dt*(1:ifinal)-dt,Agents_log{i}.y_hat(dimp*i*(i<=N_dvl)+dimp*N_dvl*(i>N_dvl)+j+(i-1)*(NAgents+NBeacons-1),:),'Color',colors(i,:),'LineStyle','--')
    end
end
hold off
setniceplot
xlim([0 dt*(ifinal-1)])
ylabel('distance (m)')
xlabel('t (s)')
title('Inter-vehicle range measurements')
if tr_loc_state
    print -dpdf vr_x
else
    print -dpdf vr_y
end

%beacon-vehicle distances
figure
for i=1:NAgents
    for j=1:NBeacons
        plot(dt*(1:ifinal)-dt,Agents_log{i}.y(j+dimp*(i<=N_dvl),:),'Color',colors(i,:))
        hold on
        plot(dt*(1:ifinal)-dt,Agents_log{i}.y_hat(dimp*i*(i<=N_dvl)+dimp*N_dvl*(i>N_dvl)+j+(i-1)*(NAgents+NBeacons-1),:),'Color',colors(i,:),'LineStyle','--')
    end
end
hold off
setniceplot
xlim([0 dt*(ifinal-1)])
xlabel('distance (m)')
xlabel('t (s)')
title('Beacon-vehicle range measurements')
if tr_loc_state
    print -dpdf br_x
else
    print -dpdf br_y
end

%Estimation error norm
figure
for i=1:NAgents
    Eerr = Agents_log{i}.x-Output_select{i}*Agents_log{i}.x_hat;
    plot(dt*(1:ifinal)-dt,sqrt(sum(Eerr.^2,1)),'Color',colors(i,:))
    hold on
end
hold off
setniceplot
xlim([0 dt*(ifinal-1)])
ylabel('estimation error (m)')
xlabel('t (s)')
if tr_loc_state
    print -dpdf ee_x
else
    print -dpdf ee_y
end

%Lambda
figure
if tr_loc_state
    %for i=1:NAgents
    i=1;
    for j=1:dimp
        plot(dt*(1:ifinal)-dt,Agents_log{i}.Lambda(j,:),'b-')
        hold on
        plot(dt*(1:ifinal)-dt,Agents_log{i}.Lambda(j+dimp,:),'r-')
    end
    %end
else
    i=1;
    %for i=1:NAgents
    for j=((i<= N_dvl)*dimp+NBeacons):size(Agents_log{i}.Lambda,1)
       plot(dt*(1:ifinal)-dt,Agents_log{i}.Lambda(j,:),'b-')
       hold on
    end
    for j=(i<= N_dvl)*dimp+(1:NBeacons)
        plot(dt*(1:ifinal)-dt,Agents_log{i}.Lambda(j,:),'g-')
        hold on
    end
    %end
    %for i=1:N_dvl
    for j=1:dimp
        plot(dt*(1:ifinal)-dt,Agents_log{i}.Lambda(j+dimp,:),'r-')
        hold on
    end
    %end
end
hold off
setniceplot
xlim([0 dt*(ifinal-1)])
set(gca, 'YScale', 'log')
ylabel('Quant. interval length')
xlabel('t (s)')
if tr_loc_state
    print -dpdf lambda_x
else
    print -dpdf lambda_y
end

err=0;
for i=1:NAgents
    err=err+sum(sum((Agents_log{i}.x-Output_select{i}*Agents_log{i}.x_hat).^2))/(ifinal*NAgents);
end
fprintf('The error is %f\n',err);


% Position errors
% Compute position error
p_err=cell(NAgents,1);
p_err_norm=cell(NAgents,1);

for i=1:NAgents
    p_err{i}=zeros(dimp,ifinal);
    p_err_norm{i}=zeros(1,ifinal);
end

for j=1:ifinal
    for i=1:NAgents
        p_err{i}(:,j)=Agents_log{i}.x(1:dimp,j)-pd((j-1)*dt)+Dd{1,i};
        p_err_norm{i}(j)=norm(p_err{i}(:,j));
    end
end

figure
for i=1:NAgents
    plot(dt*(1:ifinal)-dt,p_err_norm{i},'Color',colors(i,:))
    hold on
end
hold off
setniceplot
xlim([0 dt*(ifinal-1)])
ylabel('Control error (m)')
xlabel('t (s)')
if tr_loc_state
    print -dpdf ce_x
else
    print -dpdf ce_y
end