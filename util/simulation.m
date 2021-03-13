%% Initialization
disp('Initialization')
% Initialize agents
Agents = cell(NAgents,1);
Agents_broadcast = cell(NAgents,1);
for i=1:NAgents
    Agents{i}.x=[p_nominal(:,i);speed]+diag([xstd*ones(1,dimp) sstd*ones(1,dimp)])*randn(dim,1);
    %Observer
    Agents{i}.Observer = EKF_range(A,B,p_beacon,N_dvl,wstd_o,vstd_o,dvlstd_o,xstd,sstd,p_nominal,speed);
    %Control
    Agents{i}.Controller = Controller(speed,Freq,Amp_pert,p_nominal,Ks,Kp);
end

% Initalize log
Agents_log = cell(NAgents,1);
for i=1:NAgents
    Agents_log{i}.x=zeros(dim,ifinal);
    Agents_log{i}.x_hat=zeros(NAgents*dim,ifinal);
    Agents_log{i}.u=zeros(NAgents*dimp,ifinal);
    Agents_log{i}.y=zeros(NAgents+NBeacons-1+2*(i==1),ifinal);
    Agents_log{i}.y_hat=zeros(2+NAgents*(NAgents+NBeacons-1),ifinal);
end

% Messages between agents
Messages = cell(NAgents,1);

%% Simulation
disp('Simulation')

for j=1:ifinal

    %Measure and Broadcast
    for i=1:NAgents
        Agents{i}.y = measurement(i,Agents,vstd,dvlstd,p_beacon,N_dvl);
        Messages{i} = Agents{i}.y;
    end
    
    %Update and Control
    for i=1:NAgents
        % update
        Agents{i}.Observer.Update(Messages);
        % control
        if debug_observer
            x_cont=zeros(NAgents*dim,1);
            for k=1:NAgents
                x_cont=x_cont+Output_select{k}'*Agents{k}.x;
            end
            Agents{i}.Controller.update_control(x_cont,(j-1)*dt);
        else
            Agents{i}.Controller.update_control(Agents{i}.Observer.x_hat,(j-1)*dt);
        end
    end
    
    %Log
    for i=1:NAgents
        Agents_log{i}.x(:,j)=Agents{i}.x;
        Agents_log{i}.x_hat(:,j)=Agents{i}.Observer.x_hat;
        Agents_log{i}.u(:,j)=Agents{i}.Controller.u;
        Agents_log{i}.y(:,j)=Agents{i}.y;
        Agents_log{i}.y_hat(:,j)=Agents{i}.Observer.y_hat;
    end
    
    %Propagate and Predict
    for i=1:NAgents
        % Propagate
        Agents{i}.x=A*Agents{i}.x+B*POutput_select{i}*Agents{i}.Controller.u+wstd*randn(dim,1);
        % Predict
        Agents{i}.Observer.Predict(Agents{i}.Controller.u);
    end
end