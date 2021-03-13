function sim_output = sim_cons(sim_input)

alpha = sim_input.alpha ;
nb = sim_input.nb;
lf = sim_input.lf;
r = sim_input.r;
ifinal = sim_input.ifinal;
NAg = sim_input.NAg;
nx = sim_input.nx;
ny = sim_input.ny;
nu = sim_input.nu;
A = sim_input.A;
Pi = sim_input.Pi;
Bglob = sim_input.Bglob;
Cglob = sim_input.Cglob;
Kglobal = sim_input.Kglobal;
Lglobal = sim_input.Lglobal;
Cdist = sim_input.Cdist;
Kdist = sim_input.Kdist;
Ldist = sim_input.Ldist;
wbloc = sim_input.wbloc;
vbloc = sim_input.vbloc;

% initialization
xbar=zeros(NAg*nx,1);
xdhat=zeros(NAg*nx,1);
xhat=zeros(nx,1);
x=zeros(nx,1);
w=gen_noise(nx,wbloc);
v=gen_noise((2*NAg-1)*ny,vbloc);
u=zeros(nu,1);
y=v;

% log variables
xbarlog=zeros(NAg*nx,ifinal);
xdhatlog=zeros(NAg*nx,ifinal);
xhatlog=zeros(nx,ifinal);
xlog=zeros(nx,ifinal);
wlog=zeros(nx,ifinal);
vlog=zeros((2*NAg-1)*ny,ifinal);
ulog=zeros(nu,ifinal);
ylog=zeros((2*NAg-1)*ny,ifinal);

xbarlog(:,1)=xbar;
xdhatlog(:,1)=xdhat;
xhatlog(:,1)=xhat;
xlog(:,1)=x;
wlog(:,1)=w;
vlog(:,1)=v;
ulog(:,1)=u;
ylog(:,1)=y;

%cicle
for i=2:ifinal
    %disp(i)
    xdhat=kron(eye(NAg),A+Bglob*Kglobal)*xdhat+Ldist*(y-Cdist*xdhat);
    %xdhat=NAg^-1*kron(ones(NAg),eye(nx))*xdhat;
    %xdhat=consensus_alg(xdhat,Pi,nx,lf);
    [xdhat,xbar]=q_consensus_alg(xdhat,xbar,Pi,nx,lf,nb,alpha,r);
    xbar=kron(eye(NAg),A+Bglob*Kglobal)*xbar;
    xhat=A*xhat+Bglob*u+Lglobal*(y-Cglob*xhat);
    x=A*x+Bglob*u+w;
    w=gen_noise(nx,wbloc);
    v=gen_noise((2*NAg-1)*ny,vbloc);
    u=Kdist*xdhat;
    y=Cglob*x+v;
    
    %log variables
    xbarlog(:,i)=xbar;
    xdhatlog(:,i)=xdhat;
    xhatlog(:,i)=xhat;
    xlog(:,i)=x;
    wlog(:,i)=w;
    vlog(:,i)=v;
    ulog(:,i)=u;
    ylog(:,i)=y;
end

sim_output = struct;

sim_output.xbarlog=xbarlog;
sim_output.xdhatlog=xdhatlog;
sim_output.xhatlog=xhatlog;
sim_output.xlog=xlog;
sim_output.wlog=wlog;
sim_output.vlog=vlog;
sim_output.ulog=ulog;
sim_output.ylog=ylog;

end