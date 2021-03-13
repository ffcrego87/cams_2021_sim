classdef EKF_range < handle
   properties
       dim
       NAg
       NB
       dimp
       Output_select
       A
       B
       p_beacon
       N_dvl
       Q
       vstd_d
       dvlstd
       x_hat
       P
       y
   end
   properties (Dependent)
       C
       R
       y_hat
   end
   methods
       function obj=EKF_range(A,B,p_beacon,N_dvl,wstd,vstd_d,dvlstd,xstd,sstd,p_nominal,speed)
           obj.dimp = size(p_nominal,1);
           obj.dim = obj.dimp*2;
           obj.NAg = size(p_nominal,2);
           obj.NB = size(p_beacon,2);
           obj.Output_select = cell(obj.NAg,1);
           for i=1:obj.NAg
               obj.Output_select{i}=zeros(obj.dim,obj.NAg*obj.dim);
               obj.Output_select{i}(:,(obj.dim*(i-1)+1):(obj.dim*i)) = eye(obj.dim);
           end
           obj.A = kron(eye(obj.NAg),A);
           obj.B = kron(eye(obj.NAg),B);
           obj.p_beacon = p_beacon;
           obj.N_dvl = N_dvl;
           obj.Q = kron(eye(obj.NAg),wstd^2*eye(obj.dim));
           obj.dvlstd = dvlstd;
           obj.vstd_d = vstd_d;
           obj.x_hat = zeros(obj.NAg*obj.dim,1);
           obj.y = zeros(2+obj.NAg*(obj.NAg-1),1);
           for j=1:obj.NAg
               obj.x_hat=obj.x_hat+obj.Output_select{j}'*[p_nominal(:,j);speed];
           end
           obj.P = kron(eye(obj.NAg),diag([xstd*ones(1,obj.dimp) sstd*ones(1,obj.dimp)]));
       end
       function Update(obj,Messages)
           obj.Compute_Meas(Messages)
           C_aux = obj.C;
           y_hat_aux = obj.y_hat;
           S=C_aux*obj.P*C_aux'+obj.R;
           K=obj.P*C_aux'*S^-1;
           obj.x_hat = obj.x_hat + K*(obj.y-y_hat_aux);
           obj.P = (eye(size(obj.P))-K*obj.C)*obj.P;
       end
       function Predict(obj,u)
           obj.x_hat = obj.A*obj.x_hat+obj.B*u;
           obj.P = obj.A*obj.P*obj.A'+obj.Q;
       end
       function Compute_Meas(obj,Messages)
           idx=0;
           for i=1:obj.NAg
               obj.y(idx+(1:size(Messages{i},1)))=Messages{i};
               idx = idx+size(Messages{i},1);
           end
       end
       function y_hat=get.y_hat(obj)
           idx=0;
           y_hat=zeros(2+obj.NAg*(obj.NAg+obj.NB-1),1);
           for i=1:obj.NAg
               %auxiliary variables
               xi=obj.Output_select{i}*obj.x_hat;
               pos_o=xi(1:obj.dimp);
               
               %initialization
               y_loc=zeros((i==1)+obj.NAg,1);
               
               %dvl measurements
               if i==1
                   y_loc(1:obj.dimp) = xi((obj.dimp+1):end,1);
                   %disp(y_loc(1:obj.dimp))
               end
               
               %leader range
               for j = 1:obj.NB
                   y_loc((i==1)*obj.dimp+j) = abs(norm(pos_o-obj.p_beacon(:,j)));%+norm(pos_o-p_beacon)*vstd_d*randn);
               end
               
               %ranges
               idxr = (i==1)*obj.dimp+obj.NB+1;
               for j=1:obj.NAg
                   if j~=i
                       pos_j=obj.Output_select{j}(1:obj.dimp,:)*obj.x_hat;
                       y_loc(idxr) = norm(pos_o-pos_j);
                       idxr=idxr+1;
                   end
               end
               
               y_hat(idx+(1:size(y_loc,1)))=y_loc;
               idx = idx+size(y_loc,1);
           end
       end
       function C=get.C(obj)
           idx=0;
           C=zeros(2+obj.NAg*(obj.NAg+obj.NB-1),obj.NAg*obj.dim);
           for i=1:obj.NAg
               %auxiliary variables
               xi=obj.Output_select{i}*obj.x_hat;
               pos_o=xi(1:obj.dimp);
               
               %initialization
               C_loc=zeros((i==1)+obj.NAg,obj.NAg*obj.dim);
               
               %dvl measurements
               if i==1
                   C_loc(1:obj.dimp,:) = obj.Output_select{i}((obj.dimp+1):end,:);
               end
               
               %leader range
               for j = 1:obj.NB
                   C_loc((i==1)*obj.dimp+j,:) = (pos_o-obj.p_beacon(:,j))'*obj.Output_select{i}(1:obj.dimp,:)*norm(pos_o-obj.p_beacon(:,j))^-1;
               end
               
               %ranges
               idxr = (i==1)*obj.dimp+obj.NB+1;
               for j=1:obj.NAg
                   if j~=i
                       pos_j=obj.Output_select{j}(1:obj.dimp,:)*obj.x_hat;
                       C_loc(idxr,:) = (pos_o-pos_j)'*(obj.Output_select{i}(1:obj.dimp,:)-obj.Output_select{j}(1:obj.dimp,:))*norm(pos_o-pos_j)^-1;
                       idxr=idxr+1;
                   end
               end
               
               C(idx+(1:size(C_loc,1)),:)=C_loc;
               idx = idx+size(C_loc,1);
           end
       end
       function R=get.R(obj)
           idx=0;
           R=zeros(2+obj.NAg*(obj.NAg+obj.NB-1),2+obj.NAg*(obj.NAg+obj.NB-1));
           for i=1:obj.NAg
               %auxiliary variables
               xi=obj.Output_select{i}*obj.x_hat;
               pos_o=xi(1:obj.dimp);
               
               %initialization
               R_loc=zeros((i==1)*obj.dimp+obj.NAg,(i==1)*obj.dimp+obj.NAg);
               
               %dvl measurements
               if i==1
                   R_loc(1:obj.dimp,1:obj.dimp) = obj.dvlstd^2*eye(obj.dimp);
               end
               
               %leader range
               for j = 1:obj.NB
                   R_loc((i==1)*obj.dimp+j,(i==1)*obj.dimp+j) = obj.vstd_d^2;%(norm(pos_o-obj.p_beacon)*obj.vstd_d)^2;
               end
               
               %ranges
               idxr = (i==1)*obj.dimp+obj.NB+1;
               for j=1:obj.NAg
                   if j~=i
                       pos_j=obj.Output_select{j}(1:obj.dimp,:)*obj.x_hat;
                       R_loc(idxr,idxr) = obj.vstd_d^2;%(norm(pos_o-pos_j)*obj.vstd_d)^2;
                       idxr=idxr+1;
                   end
               end
               
               R(idx+(1:size(R_loc,1)),idx+(1:size(R_loc,1)))=R_loc;
               idx = idx+size(R_loc,1);
           end
       end
       
   end
end

