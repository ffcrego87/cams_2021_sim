function y = measurement(id,Agents,vstd_d,v_distance_flag,dvlstd,p_beacon,N_dvl)
%auxiliary variables
dimp=size(p_beacon,1);
xi=Agents{id}.x;
pos_o=xi(1:dimp);

%initialization
y=zeros(dimp*(id<=N_dvl)+size(Agents,1)+size(p_beacon,2)-1,1);

%dvl measurements
if id<=N_dvl
    y(1:dimp) = xi((dimp+1):end,1)+dvlstd*randn(dimp,1);
end

%leader range
for i = 1:size(p_beacon,2)
    if v_distance_flag
        y((id<=N_dvl)*dimp+i) = abs(norm(pos_o-p_beacon(:,i))+norm(pos_o-p_beacon)*vstd_d*randn);
    else
        y((id<=N_dvl)*dimp+i) = abs(norm(pos_o-p_beacon(:,i))+vstd_d*randn);
    end
end

%ranges
idx = (id<=N_dvl)*dimp+size(p_beacon,2)+1;
for j=1:size(Agents,1)
    if j~=id
        pos_j=Agents{j}.x(1:dimp);
        if v_distance_flag
            y(idx) = abs(norm(pos_o-pos_j)+norm(pos_o-pos_j)*vstd_d*randn);
        else
            y(idx) = abs(norm(pos_o-pos_j)+vstd_d*randn);
        end
        idx=idx+1;
    end
end
end

