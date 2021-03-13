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
xlabel('X')
ylabel('Y')
title('Position')
%print -dpdf

%speed
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
title('Speed')
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
xlabel('X')
ylabel('Y')
title('Control')
%print -dpdf

%leader speed
figure
plot(dt*(1:ifinal)-dt,Agents_log{1}.y(1,:),'b-');
hold on
plot(dt*(1:ifinal)-dt,Agents_log{1}.y(2,:),'r-');
hold off
setniceplot
ylabel('v')
xlabel('t')
title('Leader speed')
%print -dpdf

%inter-vehicle distances
figure
for i=1:NAgents
    for j=(2:NAgents)-1+NBeacons
        plot(dt*(1:ifinal)-dt,Agents_log{i}.y(j+2*(i==1),:),'Color',colors(i,:))
        hold on
        plot(dt*(1:ifinal)-dt,Agents_log{i}.y_hat(2+j+(i-1)*(NAgents+NBeacons-1),:),'Color',colors(i,:),'LineStyle','--')
    end
end
hold off
setniceplot
ylabel('dist')
xlabel('t')
title('Inter-vehicle range measurements')

%beacon-vehicle distances
figure
for i=1:NAgents
    for j=1:NBeacons
        plot(dt*(1:ifinal)-dt,Agents_log{i}.y(j+2*(i==1),:),'Color',colors(i,:))
        hold on
        plot(dt*(1:ifinal)-dt,Agents_log{i}.y_hat(2+j+(i-1)*(NAgents+NBeacons-1),:),'Color',colors(i,:),'LineStyle','--')
    end
end
hold off
setniceplot
xlabel('dist')
xlabel('t')
title('Beacon-vehicle range measurements')