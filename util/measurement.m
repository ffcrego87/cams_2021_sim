function y = measurement(id,Agents,vstd_d,dvlstd,p_beacon,N_dvl)
%auxiliary variables
dimp=size(p_beacon,1);
xi=Agents{id}.x;
pos_o=xi(1:dimp);

%initialization
y=zeros((id==1)+size(Agents,1),1);

%dvl measurements
if id==1
    y(1:dimp) = xi((dimp+1):end,1)+dvlstd*randn(dimp,1);
end

%leader range
for i = 1:size(p_beacon,2)
    y((id==1)*dimp+i) = abs(norm(pos_o-p_beacon(:,i)))+vstd_d*randn;%norm(pos_o-p_beacon)*vstd_d*randn);
end

%ranges
idx = (id==1)*dimp+size(p_beacon,2)+1;
for j=1:size(Agents,1)
    if j~=id
        pos_j=Agents{j}.x(1:dimp);
        y(idx) = abs(norm(pos_o-pos_j))+vstd_d*randn;%norm(pos_o-pos_j)*vstd_d*randn);
        idx=idx+1;
    end
end
end

