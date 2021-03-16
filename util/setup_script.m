%% System
disp('Setup script')

% Number of agents
NAgents = size(p_nominal,2);
NBeacons = size(p_beacon,2);

% Position dimension
dimp = dim/2;

%Auxiliar selection matrices
Output_select = cell(NAgents,1);
for i=1:NAgents
    Output_select{i}=zeros(dim,NAgents*dim);
    Output_select{i}(:,(dim*(i-1)+1):(dim*i)) = eye(dim);
end

POutput_select = cell(NAgents,1);
for i=1:NAgents
    POutput_select{i}=zeros(dimp,NAgents*dimp);
    POutput_select{i}(:,(dimp*(i-1)+1):(dimp*i)) = eye(dimp);
end

Dd = cell(NAgents,NAgents);
for i=1:NAgents
    for j=1:NAgents
        Dd{i,j}=p_nominal(:,i)-p_nominal(:,j);
    end
end