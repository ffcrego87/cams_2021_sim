%% Initialization
disp('Initialization')
% Initialize agents
Agents = cell(NAgents,1);
Agents_broadcast = cell(NAgents,1);
for i=1:NAgents
    Agents{i}.x=[p_nominal(:,i);zeros(dimp,1)]+diag([xstd*ones(1,dimp) sstd*ones(1,dimp)])*randn(dim,1);
    %Observer
    if tr_loc_state
        Agents{i}.Observer = EKF_loc_state(i,A,B,p_beacon,N_dvl,wstd_o,vstd_o,v_distance_flag_o,dvlstd_o,xstd_o,sstd_o,p_nominal,speed,alpha,beta,Lambda_zero,Nb);
    else
        Agents{i}.Observer = EKF_range(A,B,p_beacon,N_dvl,wstd_o,vstd_o,v_distance_flag_o,dvlstd_o,xstd_o,sstd_o,p_nominal,speed,alpha,beta,Lambda_zero,Nb);
    end
    %Control
    Agents{i}.Controller = Controller(speed,Freq,Amp_pert,p_nominal,Ks,Kp);
end

% Initalize log
Agents_log = cell(NAgents,1);
for i=1:NAgents
    Agents_log{i}.x=zeros(dim,ifinal);
    Agents_log{i}.x_hat=zeros(NAgents*dim,ifinal);
    Agents_log{i}.u=zeros(NAgents*dimp,ifinal);
    Agents_log{i}.y=zeros(NAgents+NBeacons-1+dimp*(i<=N_dvl),ifinal);
    Agents_log{i}.y_hat=zeros(dimp*N_dvl+NAgents*(NAgents+NBeacons-1),ifinal);
    if tr_loc_state
        Agents_log{i}.Lambda=zeros(dim,ifinal);
    else
        Agents_log{i}.Lambda=zeros(NAgents+NBeacons-1+dimp*(i<=N_dvl),ifinal);
    end
end

% Messages between agents
Messages = cell(NAgents,1);

%% Simulation
disp('Simulation')

for j=1:ifinal

    %Measure and Broadcast / Measure, Update and Broadcast
    for i=1:NAgents
        %Measure
        Agents{i}.y = measurement(i,Agents,vstd,v_distance_flag,dvlstd,p_beacon,N_dvl);
        %Update
        if tr_loc_state
            Agents{i}.Observer.Update(Agents{i}.y)
        end
        %Broadcast
        if tr_loc_state
            Messages{i} = Agents{i}.Observer.Broadcast();
        else
            Messages{i} = Agents{i}.Observer.Broadcast(Agents{i}.y,i);
        end
    end
    
    %Update and Control / Fuse and Control
    for i=1:NAgents
        % update / Fuse
        if tr_loc_state
            Agents{i}.Observer.Fuse(Messages);
        else
            Agents{i}.Observer.Update(Messages);
        end
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
        Agents_log{i}.Lambda(:,j)=Agents{i}.Observer.Lambda{i};
    end
    
    %Propagate and Predict
    for i=1:NAgents
        % Propagate
        Agents{i}.x=A*Agents{i}.x+B*POutput_select{i}*Agents{i}.Controller.u+wstd*randn(dim,1);
        % Predict
        Agents{i}.Observer.Predict(Agents{i}.Controller.u);
    end
end