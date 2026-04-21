clear all
close all
clc

n=6; % number of states
nv=2; % number of exostates 
m=3; % number of controls
%% Time step
T = 20;
tspan = 0:0.01:T;
%% Control Gains
% Problem Definition
Cons = deg2rad(45);
Consdeg = 45;

R = eye(m);
Q = diag([10 10 10 2 2 2]);

params.etac1 = 0.1;
params.etac2 = 1;
params.etaa1 = 1;
params.etaa2 = 0.1;
params.beta = 0.01;
params.nu = 1;
params.cb = 0.01;
% System ID
params.k = 0.01;
params.l = 0.1;
params.GammaTheta = 10;
%% Initial Conditions stack
[~,L]=Basis(zeros(n+nv,1));
v0 = deg2rad([0 40]');
x0 = deg2rad([21, 28, 28, 0, 0, 0])';
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
options = odeset('OutputFcn',@odeplot,'OutputSel',[1 2 3]);
tic
[t,y] = ode45(@(t,y) closedLoopDynamics(t,y,n,nv,L,p,Q,R,params,cldata,Cons),tspan,z0,options);
toc

[~,ADPRank,U,Mu,Udhat,MEP,Bcell,Usafe] = cellfun(@(t,y)closedLoopDynamics(t,y.', ...
    n,nv,L,p,Q,R,params,cldata,Cons), num2cell(t), num2cell(y,2), 'uni',0);
%% Plot
ADPTr_Plot
% % 保存动画所需数据
% % t 和 y 已经是 ode45 的输出
% % 提取外系统状态 v（y 的第 n+1 到 n+nv 列）
% v_sim = y(:, n+1:n+nv);
% % 实际欧拉角：x 的前3个状态 (phi, theta, psi)
% phi_act   = y(:, 1);
% theta_act = y(:, 2);
% psi_act   = y(:, 3);
% % 参考欧拉角：根据外系统定义，xd = [v1; v1; v1; v2; v2; v2]，所以参考角度为 v_sim(:,1)
% phi_ref   = v_sim(:, 1);
% theta_ref = v_sim(:, 1);   % 与 phi_ref 相同，您可根据需要修改
% psi_ref   = v_sim(:, 1);   % 与 phi_ref 相同
% 
% % 保存到文件
% save('sim_data.mat', 't', 'phi_act', 'theta_act', 'psi_act', ...
%                      'phi_ref', 'theta_ref', 'psi_ref');
%--------------------------------------------------------------------------
function [ydot,ra,u,mu,udhat,meP,B,usafe]=closedLoopDynamics(t,y, ...
    n,nv,L,p,QT,R,params,cldata,Cons)
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

vdot=a_fun(v); 
xd=[v(1); v(1); v(1); v(2); v(2); v(2)];  % exosystem: \dot{v}=a(v), xd=hd(v)
xd_dot = [v(2); v(2); v(2); -v(1); -v(1); -v(1)];  % ∂hd(v)/∂v a(v)

e = x-xd;
Zeta = [e;v];

% Dynamics
Jx = 0.0211;
Jy = 0.0219;
Jz = 0.0366;
J = diag([Jx,Jy,Jz]);
g = [zeros(3,3); J^(-1)];
gplusd = [zeros(3,3) J];

phi = SIDBasis(x);
phid= SIDBasis(xd);
FTH=[thetaH'*phi-g*gplusd*thetaH'*phid;zeros(nv,1)];
F1 = [f0_fun(x)-xd_dot+g*gplusd*xd_dot;vdot];  % known part defined as F1

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
beta_cbf=5;
h01 = (Cons+x(1))*(Cons-x(1));
h1 = -2*x(1)*x(4) + beta_cbf*h01;
nabla_h1 = [-2*x(4)-2*beta_cbf*x(1); 0; 0; -2*x(1); 0; 0];

h02 = (Cons+x(2))*(Cons-x(2));
h2 = -2*x(2)*x(5) + beta_cbf*h02;
nabla_h2 = [0; -2*x(5)-2*beta_cbf*x(2); 0; 0; -2*x(2); 0];

h03 = (Cons+x(3))*(Cons-x(3));
h3 = -2*x(3)*x(6) + beta_cbf*h03;
nabla_h3 = [0; 0; -2*x(6)-2*beta_cbf*x(3); 0; 0; -2*x(3)];

cbf1=1/h1;  cbf2=1/h2;  cbf3=1/h3;
B = norm(e)^2*(cbf1+cbf2+cbf3);
gradxB = 2*e*(cbf1+cbf2+cbf3) + norm(e)^2*(-nabla_h1/h1^2 -nabla_h2/h2^2 -nabla_h3/h3^2);

usafe = -cb*g'*gradxB;

% Controller
u = udhat+mu+usafe;

% Identifier
f0 = f0_fun(x);
[thetaHD,dxf,dphif,dvarphif,dP,dQ,meP] = identifier(thetaH,u,t,x, ...
    xf,phif,varphif,f0,P,Q, k,l,GammaTheta);

xD = f0+g*(u)+g*W_fun(x);

ydot=[xD;vdot;WcHD;WaHD;GammaD;
    thetaHD;dxf;dphif;dvarphif;reshape(dP,p*p,1);reshape(dQ,p*n,1)];
end


function a = f0_fun(x)
Jx = 0.0211;
Jy = 0.0219;
Jz = 0.0366;
J = diag([Jx,Jy,Jz]);

phi = x(1);
theta = x(2);
psi = x(3);
dphi = x(4);
dtheta = x(5);
dpsi = x(6);

c11 = 0;
c12 = (Jy-Jz)*(dtheta*cos(phi)*sin(phi)+dpsi*sin(phi)^2*cos(theta)) + (Jz-Jy)*dpsi*cos(phi)^2*cos(theta)...
    -Jx*dpsi*cos(theta);
c13 = (Jz-Jy)*dpsi*cos(phi)*sin(phi)*cos(theta)^2;
c21 = (Jz-Jy)*(dtheta*cos(phi)*sin(phi)+dpsi*sin(phi)^2*cos(theta)) + (Jy-Jz)*dpsi*cos(phi)^2*cos(theta)...
    +Jx*dpsi*cos(theta);
c22 = (Jz-Jy)*dphi*cos(phi)*sin(phi);
c23 = -Jx*dpsi*sin(theta)*cos(theta) +Jy*dpsi*sin(phi)^2*cos(theta)*sin(theta) ...
    +Jz*dpsi*cos(phi)^2*sin(theta)*cos(theta);
c31 = (Jy-Jz)*dpsi*cos(theta)^2*sin(phi)*cos(phi) -Jx*dtheta*cos(theta);
c32 = (Jz-Jy)*(dtheta*cos(phi)*sin(phi)*sin(theta)+dphi*sin(phi)^2*cos(theta))...
    +(Jy-Jz)*dphi*cos(phi)^2*cos(theta)...
    +Jx*dpsi*sin(theta)*cos(theta) -Jy*dpsi*sin(phi)^2*sin(theta)*cos(theta)...
    -Jz*dpsi*cos(phi)^2*sin(theta)*cos(theta);
c33 = (Jy-Jz)*dphi*cos(phi)*sin(phi)*cos(theta)^2 -Jy*dtheta*sin(phi)^2*cos(theta)*sin(theta)...
    -Jz*dtheta*cos(phi)^2*cos(theta)*sin(theta) +Jx*dtheta*cos(theta)*sin(theta);
C = [c11 c12 c13;
    c21 c22 c23;
    c31 c32 c33];

a = [x(4:6); 
    J\(-C*x(4:6))];
end

function W = W_fun(x)
rho = [0.6294 0.8116 -0.7460 0.8268 0.2647 -0.8049];

phi = x(1);
theta = x(2);
psi = x(3);
dphi = x(4);
dtheta = x(5);
dpsi = x(6);

W = [rho(1)*phi*cos(dtheta)+rho(2)*theta*sin(phi)
    rho(3)*psi*sin(phi)+rho(4)*dphi*cos(theta)
    rho(5)*dtheta*cos(dpsi)+rho(6)*dpsi*sin(theta)];
end

function vdot = a_fun(v)
vdot = [0 1
       -1 0]*v;
end

function cldata=create_CL_term(n,nv,L,p,Q,R)
rng(0) % To reproduce
N=100;
XC=unifrnd(-0.5,0.5,n+nv,N);

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
    Zc=XC(:,i);  % extrapolated points
    ec=Zc(1:n,1);  % error
    vc=Zc(n+1:n+nv,1); % exosystem
    
    vdotc = a_fun(vc);  
    xdc=[vc(1); vc(1); vc(1); vc(2); vc(2); vc(2)];  % exosystem: \dot{v}=a(v), xd=hd(v)
    xd_dotc = [vc(2); vc(2); vc(2); -vc(1); -vc(1); -vc(1)];  % ∂hd(v)/∂v a(v)
    xc=ec+xdc;   % state
    Zc=[ec;vc];
    
    Jx = 0.0211;
    Jy = 0.0219;
    Jz = 0.0366;
    J = diag([Jx,Jy,Jz]);
    gc = [zeros(3,3); J^(-1)];
    gplusdc = [zeros(3,3) J];
    Gc = [gc;zeros(nv,size(gc,2))];  %Gc=[gc;0nv×m]

    % Fc=[f_fun(xc)-gc*gplusdc*f_fun(xdc);zeros(nv,1)];
    F1c = [f0_fun(xc)-xd_dotc+gc*gplusdc*xd_dotc;vdotc];  % known part defined as F1

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