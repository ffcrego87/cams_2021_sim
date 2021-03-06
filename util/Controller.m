classdef Controller < handle
   properties
       dim
       dimp
       NAg
       Ks
       Kp
       Dd
       pos_cell
       v_cell
       Output_select
       POutput_select
       u
       pd
       dpd
   end
   methods
       function obj=Controller(p_nominal,Ks,Kp,pd,dpd)
           obj.dimp = size(p_nominal,1);
           obj.dim = obj.dimp*2;
           obj.NAg = size(p_nominal,2);
           obj.Ks = Ks;
           obj.Kp = Kp;
           obj.pd = pd;
           obj.dpd = dpd;
           obj.Dd = cell(obj.NAg,obj.NAg);
           for i=1:obj.NAg
               for j=1:obj.NAg
                   obj.Dd{i,j}=p_nominal(:,i)-p_nominal(:,j);
               end
           end
           obj.pos_cell = cell(obj.NAg,1);
           obj.v_cell = cell(obj.NAg,1);
           obj.Output_select = cell(obj.NAg,1);
           for i=1:obj.NAg
               obj.Output_select{i}=zeros(obj.dim,obj.NAg*obj.dim);
               obj.Output_select{i}(:,(obj.dim*(i-1)+1):(obj.dim*i)) = eye(obj.dim);
           end
           obj.POutput_select = cell(obj.NAg,1);
           for i=1:obj.NAg
               obj.POutput_select{i}=zeros(obj.dimp,obj.NAg*obj.dimp);
               obj.POutput_select{i}(:,(obj.dimp*(i-1)+1):(obj.dimp*i)) = eye(obj.dimp);
           end
           obj.u = zeros(obj.NAg*obj.dimp,1);
       end
       function update_control(obj,x_hat,t)
           % Retrieve positions from observer
           for i=1:obj.NAg
               xi=obj.Output_select{i}*x_hat;
               obj.pos_cell{i}=xi(1:obj.dimp,1);
               obj.v_cell{i}=xi(obj.dimp+(1:obj.dimp),1);
           end
           
           obj.u = zeros(obj.NAg*obj.dimp,1);
           
           % Leader controller
           i=1;
           %pd=obj.p0+obj.vd*t+obj.Pert_basis*sin(2*pi*obj.Freq*t);
           ep=obj.pos_cell{i}-obj.pd(t);
           p_dot_d=obj.dpd(t)-obj.Kp*ep;
           es=obj.v_cell{i}-p_dot_d;
           ui=-obj.Ks*es;
           obj.u=obj.u+obj.POutput_select{i}'*ui;
           
           % Formation controllers
           for i=2:obj.NAg
               ep=zeros(obj.dimp,1);
               for j=2:obj.NAg
                   ep=ep+obj.pos_cell{i}-obj.pos_cell{j}-obj.Dd{i,j};
               end
               if obj.NAg>2
                   ep=ep./(obj.NAg-2);
               end
               ep=0.3*ep+0.7*(obj.pos_cell{i}-obj.pos_cell{1}-obj.Dd{i,1});
               p_dot_d=obj.v_cell{1}-obj.Kp*ep;
               es=obj.v_cell{i}-p_dot_d;
               ui=-obj.Ks*es;
               obj.u=obj.u+obj.POutput_select{i}'*ui;
           end
       end
   end
end

