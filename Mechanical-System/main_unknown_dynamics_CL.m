clear all
close all
clc

n=2; % number of states
nv=2; % number of exostates 
m=1; % number of controls
%% Control Gains
% Problem Definition
T = 100;
tspan=0:0.01:T;

Env.c1=[0;0]; Env.r1=3;
Env.c2=[0;0]; Env.r2=sqrt(0.5);

R=1;
Q=diag([1 1]);

params.etac1 = 0.1;
params.etac2 = 0.75;
params.etaa1 = 0.75;
params.etaa2 = 0.1;
params.beta = 0.01;
params.nu = 1;
% Safeguarding
params.cb = 0.01;
% System ID
params.k = 0.01;
params.l = 0.1;
params.GammaTheta = 10;
%% Initial Conditions stack
[~,L]=Basis(zeros(n+nv,1));
v0 = [2 2]';
x0 = [0 1.3]';
WcH0 = 5*ones(L,1);
WaH0 = 3*ones(L,1);
Gamma0 = 10*eye(L);

[~,p]=SIDBasis(x0);
thetaH0 = 0*ones(p,n);
xf0 = x0;
phif0 = zeros(p,1);
varphif0 = zeros(n,1);
z0 = [x0;v0;WcH0;WaH0;Gamma0(:);thetaH0(:);xf0;phif0;varphif0;zeros(p*p,1);zeros(p*n,1)];
%% Concurrent learning - Create database
cldata=create_CL_term(n,nv,L,p,Q,R);
%% Simulation
options = odeset('OutputFcn',@(t,y,flag)phaseplot(t,y,flag,Env),'OutputSel',[1 2 3 4]);
% options = odeset('OutputFcn',@odeplot,'OutputSel',[1 2 3 4]);
tic
[t,y] = ode45(@(t,y) closedLoopDynamics(t,y,n,nv,L,p,Q,R,params,cldata,Env),tspan,z0,options);
toc

[~,RA,U,Mu,Udhat,Usafe,MEP,Bcell,hcell] = cellfun(@(t,y)closedLoopDynamics(t,y.', ...
    n,nv,L,p,Q,R,params,cldata,Env), num2cell(t), num2cell(y,2), 'uni',0);
%% Plot
ADPTr_plot;
% Compute RMSE and IAU
e = y(:,1:n) - y(:,n+1:n+nv);          % tracking error e = x - xd
e_norm_sq = sum(e.^2, 2);              % ||e||^2
integral_e = trapz(t, e_norm_sq);      % ∫||e||^2 dτ
RMSE = sqrt(integral_e / T);           % RMSE

u_vec = cell2mat(U);                   % Convert cell to vector
u_abs = abs(u_vec);                    % |u| (m=1)
IAU = trapz(t, u_abs);                 % ∫|u| dτ

fprintf('RMSE = %.6f\n', RMSE);
fprintf('IAU = %.6f\n', IAU);

%--------------------------------------------------------------------------
function [ydot,ra,u,mu,udhat,usafe,meP,B,h1]=closedLoopDynamics(t,y, ...
    n,nv,L,p,QT,R,params,cldata,Env)
% parameters
etac1=params.etac1;
etac2=params.etac2;
etaa1=params.etaa1;
etaa2=params.etaa2;
beta =params.beta;
nu   =params.nu;
k    =params.k;
l    =params.l;
GammaTheta=params.GammaTheta;
cb   =params.cb;
% CL data
QQ    =cldata.QQ;
SIGPF1=cldata.SIGPF1;  % ∇σF
SIGP  =cldata.SIGP;
SIGPGD=cldata.SIGPGD;
PHI   =cldata.PHI;
PHID  =cldata.PHID;
GSIGMA=cldata.GSIGMA;  % ∇σGR^(-1) G^T ∇σ^T

x=y(1:n); v=y(n+1:n+nv);
WcH = y(n+nv+1:n+nv+L);
WaH = y(n+nv+L+1:n+nv+2*L);
Gamma = reshape(y(n+nv+2*L+1:n+nv+2*L+L*L),L,L);
thetaH=reshape(y(n+nv+2*L+L*L+1:n+nv+2*L+L*L+p*n),p,n);
xf=y(n+nv+2*L+L*L+p*n+1:2*n+nv+2*L+L*L+p*n); phif=y(2*n+nv+2*L+L*L+p*n+1:2*n+nv+2*L+L*L+p*n+p); 
varphif=y(2*n+nv+2*L+L*L+p*n+p+1:3*n+nv+2*L+L*L+p*n+p);
P=reshape(y(3*n+nv+2*L+L*L+p*n+p+1:3*n+nv+2*L+L*L+p*n+p+p*p),p,p); 
Q=reshape(y(3*n+nv+2*L+L*L+p*n+p+p*p+1:3*n+nv+2*L+L*L+2*p*n+p+p*p),p,n); 

vdot=a_fun(v);  xd=v;  % exosystem: \dot{v}=a(v), xd=hd(v)
xd_dot=a_fun(v);  % ∂hd(v)/∂v a(v)

e = x-xd;
Zeta = [e;v];

% Dynamics
g=[0; 1];
gplusd=[0 1];

phi = SIDBasis(x);
phid= SIDBasis(xd);
FTH=[thetaH'*phi-g*gplusd*thetaH'*phid;zeros(nv,1)];
F1 = [-xd_dot+g*gplusd*xd_dot;vdot];  % known part defined as F1

udhat = gplusd*(xd_dot-thetaH'*phid);  % estimated steady-state controller
G = [g;zeros(nv,size(g,2))];

sig_p = Basis(Zeta);
mu = -0.5*(R\G')*sig_p'*WaH;

% Concurrent Learning part
N = length(QQ);
thetaH1 = [thetaH zeros(p,nv)];

GSIGMAWaH = reshape(GSIGMA*WaH,L,N);

SIGPTH=reshape(permute(reshape(reshape((SIGP*thetaH1')',L,p,N),p,L,N),[2,1,3]),L,p*N);

SIGPGDTH=reshape(permute(reshape(reshape((SIGPGD*thetaH')',L,p,N),p,L,N),[2,1,3]),L,p*N);

rc = QQ + (1/4)*(GSIGMAWaH'*WaH);  % Cost function in the meshgrid (column stack: N×1)
Omegac = (SIGPF1 + SIGPTH*PHI - SIGPGDTH*PHID - (1/2)*(GSIGMAWaH));

deltac = Omegac'*WcH + rc;  % (N×1)
normm = 1./(1+nu*sum(Omegac.^2,1));  % (1×N)

clWc = -etac2*Gamma*(Omegac.*normm)*deltac/N;
clWa = etac2*GSIGMAWaH*(Omegac.*normm)'*WcH/(4*N);
clGamma = -etac2*Gamma*(Omegac.*normm.^2)*Omegac'*Gamma/N;
ra=real(min(eig( (Omegac.*normm.^2)*Omegac' )));

if mod(t,5)<=0.0001
    fprintf('t=%.4f, ADPRank=%.0f \n',t,rank(Omegac))
end

% ADP
omega = sig_p*(FTH+F1+G*mu);
Gsigma = sig_p*G*(R\G')*sig_p';
r = e'*QT*e + mu'*R*mu;  % Cost Function
delta = WcH'*omega + r;
rho = (1+nu*(omega'*omega));
WcHD = -etac1*Gamma*omega*delta/rho +clWc;
GammaD = reshape(beta*Gamma-etac1*Gamma*(omega*omega'/rho^2)*Gamma +clGamma,L*L,1);
WaHD = etaa1*(WcH-WaH)-etaa2*WaH ...
    +etac1*Gsigma'*WaH*omega'*WcH/(4*rho) +clWa;
WaHD = proj_rectangle(WaH, WaHD, -100, 100, 1);

% Add safeguarding controller==============================================
c1=Env.c1; r1=Env.r1; c2=Env.c2; r2=Env.r2; 
h1 = r1^2-norm(x-c1)^2;
nabla_h1 = -2*(x-c1);

h2 = norm(x-c2)^2-r2^2;
nabla_h2 = 2*(x-c2);

cbf1=1/h1;  cbf2=1/h2;
B = norm(e)^2*(cbf1 + cbf2);
gradxB = 2*e*(cbf1 + cbf2) + norm(e)^2*(-nabla_h1/h1^2 -nabla_h2/h2^2);

usafe = -cb*g'*gradxB;

% Controller
u = udhat+mu+usafe;

% Identifier
[thetaHD,dxf,dphif,dvarphif,dP,dQ,meP] = identifier(thetaH,u,t,x, ...
    xf,phif,varphif,P,Q, k,l,GammaTheta);

xD = f_fun(x)+g*(u);

ydot=[xD;vdot;WcHD;WaHD;GammaD;
    thetaHD;dxf;dphif;dvarphif;reshape(dP,p*p,1);reshape(dQ,p*n,1)];
end


function a = f_fun(x)
a = [x(2); -((0.8+0.2*exp(-100*abs(x(2))))*tanh(10*x(2))+x(2))-x(1)];
end

function dv=a_fun(v)
    dv=[v(2); -v(1)+(1-v(1)^2)*v(2)];
end


function cldata=create_CL_term(n,nv,L,p,Q,R)
temp=linspace(-2,2,3);  % range of the extrapolated points
[temp1, temp2, temp3, temp4] = ndgrid(temp,temp,temp,temp);
XC=[temp1(:) temp2(:) temp3(:) temp4(:)]';
N = length(temp)^(n+nv);  % number of the extrapolated points (N=na^{n+nv})

index=1;
cldata.QQ=zeros(N,1);
cldata.SIGPF1=zeros(L,N);
% cldata.SIGPF1m=zeros(L,N);
cldata.SIGP=zeros(N*L,n+nv);
cldata.SIGPGD=zeros(N*L,n);
cldata.PHI=zeros(N*p,N);
cldata.PHID=zeros(N*p,N);
cldata.GSIGMA=zeros(N*L,L);
for i=1:N
    ec=[XC(1,i) XC(2,i)]';
    vc = [XC(3,i) XC(4,i)]';
    
    vdotc = a_fun(vc);  xdc=vc; % exosystem: \dot{v}=a(v), xd=hd(v)
    xd_dotc= a_fun(vc);  % ∂hd(v)/∂v a(v)
    xc=ec+xdc;   % state
    Zc=[ec;vc];
    
    gc = [0; 1];  % Modified
    gplusdc = [0 1];  % pseudo-inverse
    Gc = [gc;zeros(nv,size(gc,2))];  %Gc=[gc;0nv×m]

    % Fc=[f_fun(xc)-gc*gplusdc*f_fun(xdc);zeros(nv,1)];
    F1c = [-xd_dotc+gc*gplusdc*xd_dotc;vdotc];  % known part defined as F1

    % F1cm=F1c+Fc;
    
    phic=SIDBasis(xc);
    phidc=SIDBasis(xdc);
    sig_pc=Basis(Zc);  % Actor-critic NN Basis ∈R^{L×(n+nv)}

    cldata.QQ(index)=ec'*Q*ec;  % QQ∈R^{N}
    cldata.SIGPF1(:,index)=sig_pc*F1c;
    % cldata.SIGPF1m(:,index)=sig_pc*F1cm;  % SIGPF1m∈R^{L×N}
    cldata.SIGP(index*L-(L-1):index*L,:)=sig_pc;  % SIGP∈R^{NL×(n+nv)}
    cldata.SIGPGD(index*L-(L-1):index*L,:)=sig_pc*[gc*gplusdc;zeros(nv,n)];
    cldata.GSIGMA(index*L-(L-1):index*L,:)=sig_pc*Gc*(R\Gc')*sig_pc';  % GSIGMA∈R^{NL×L}
    cldata.PHI(index*p-(p-1):index*p,index)=phic;
    cldata.PHID(index*p-(p-1):index*p,index)=phidc;
    index = index+1;
end
end