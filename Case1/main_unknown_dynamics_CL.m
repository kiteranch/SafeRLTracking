clear all
close all
clc

n=2; % number of states
m=1; % number of controls
%% Control Gains
% Problem Definition
T_end = 100;
tspan=0:0.01:T_end;

Env.c1=[0;0]; Env.r1=3;
Env.c2=[0;0]; Env.r2=sqrt(0.5);
R=1;
Q=diag([10 10]);

etac1 = 0.01;
etac2 = 1;
etaa1 = 1;
etaa2 = 0.01;
beta = 0.1;
nu = 1;
% GammaBar = 100000;

cb=0.01;
%% System ID
k=0.01; % 滤波器常数
l=0.1;
GammaTheta=10;
% kp=0.001;
% kTheta=1;
% M = 20; % Number of recorded points
%% Initial Conditions stack
[~,L]=Basis(zeros(2*n,1));
xd0 = [2 2]';
x0 = [0 1.3]';
WcH0 = 5*ones(L,1);
WaH0 = 3*ones(L,1);
Gamma0 = reshape(10*eye(L),L*L,1);

[~,p]=SIDBasis(x0);
xf0 = x0;
phif0 = zeros(p,1);
varphif0 = zeros(n,1);
thetaH0 = 0*ones(p*n,1);
z0 = [x0;xd0;WcH0;WaH0;Gamma0;thetaH0;xf0;phif0;varphif0;zeros(p*p,1);zeros(p*n,1);0];
%% Concurrent learning - Create database
[QQ,SIGPF1,SIGP,SIGPGD,PHI,PHID,GSIGMA,index]=create_CL_term(n,L,p,Q,R);
%% Simulation
options = odeset('OutputFcn',@(t,y,flag)phaseplot(t,y,flag,Env),'OutputSel',[1 2 3 4]);
[t,y] = ode45(@(t,y) closedLoopDynamics(t,y,n,L,p,Q,R,etac1,etac2,etaa1,etaa2,beta,nu, ...
    QQ,SIGPF1,SIGP,SIGPGD,PHI,PHID,GSIGMA, k,l,GammaTheta, cb, Env),tspan,z0,options);

% for k=1:length(t)
%     [~,ra(k)] = closedLoopDynamics(t(k),y(k,:)',n,L,p,Q,R,etac1,etac2,etaa1,etaa2,beta,nu, ...
%     QQ,SIGPF1,SIGP,SIGPGD,PHI,PHID,GSIGMA, k,l,GammaTheta,M,kp,kTheta, cb);
% end
[~,RA,U,Mu,Udhat,Usafe,MEP,Bcell,hcell,Norme] = cellfun(@(t,y)closedLoopDynamics(t,y.', ...
    n,L,p,Q,R,etac1,etac2,etaa1,etaa2,beta,nu, QQ,SIGPF1,SIGP,SIGPGD,PHI,PHID,GSIGMA, ...
    k,l,GammaTheta, cb,Env), num2cell(t), num2cell(y,2), 'uni',0);
%% Plot
ADPTr_plot;

%--------------------------------------------------------------------------
function [ydot,ra,u,mu,udhat,usafe,meP,B,h1,norme]=closedLoopDynamics(t,y, ...
    n,L,p,QT,R,etac1,etac2,etaa1,etaa2,beta,nu, ...
    QQ,SIGPF1,SIGP,SIGPGD,PHI,PHID,GSIGMA, k,l,GammaTheta, cb,Env)
x=y(1:n); xd=y(n+1:2*n);
WcH = y(2*n+1:2*n+L);
WaH = y(2*n+L+1:2*n+2*L);
Gamma = reshape(y(2*n+2*L+1:2*n+2*L+L*L),L,L);
thetaH=reshape(y(2*n+2*L+L*L+1:2*n+2*L+L*L+p*n),p,n);
xf=y(2*n+2*L+L*L+p*n+1:3*n+2*L+L*L+p*n); phif=y(3*n+2*L+L*L+p*n+1:3*n+2*L+L*L+p*n+p); 
varphif=y(3*n+2*L+L*L+p*n+p+1:4*n+2*L+L*L+p*n+p);
P=reshape(y(4*n+2*L+L*L+p*n+p+1:4*n+2*L+L*L+p*n+p+p*p),p,p); 
Q=reshape(y(4*n+2*L+L*L+p*n+p+p*p+1:4*n+2*L+L*L+2*p*n+p+p*p),p,n); 

hd=hd_fun(xd);

e = x-xd;
Zeta = [e;xd];
norme=norm(e);

% dynamics
g=[0; 1];
gplusd=[0 1];

phi = SIDBasis(x);
phid= SIDBasis(xd);
FTH=[thetaH'*phi-g*gplusd*thetaH'*phid;zeros(n,1)];
F1 = [-hd+g*gplusd*hd;hd];  %已知部分定义为F1

udhat = gplusd*(hd-thetaH'*phid);  %期望稳态控制
G = [g;zeros(size(g))];

sig_p = Basis(Zeta);
mu = -0.5*(R\G')*sig_p'*WaH;

% add safeguarding controller==============================================
c1=Env.c1; r1=Env.r1; c2=Env.c2; r2=Env.r2; 
h1 = r1^2-norm(x-c1)^2;
nabla_h1 = -2*(x-c1);
% h1d= 2-2*xd(1)^2-xd(2)^2;

h2 = norm(x-c2)^2-r2^2;
nabla_h2 = 2*(x-c2);
% h2d= 2*xd(1)^2+xd(2)^2-0.32;

% B = (1/h1-1/h1d)^2 + (1/h2-1/h2d)^2;
% gradxB = 2*(1/h1-1/h1d)*(-1/h1^2)*nabla_h1 + 2*(1/h2-1/h2d)*(-1/h2^2)*nabla_h2;

cbf1=1/h1;  cbf2=1/h2;
B = norm(e)^2*(cbf1 + cbf2);
gradxB = 2*e*(cbf1 + cbf2) + norm(e)^2*(-nabla_h1/h1^2 -nabla_h2/h2^2);
%==========================================================================
usafe = -cb*g'*gradxB;

% Controller
u = udhat+mu+usafe;

% identifier
[thetaHD,dxf,dphif,dvarphif,dP,dQ,meP] = identifier(thetaH,u,t,x, ...
    xf,phif,varphif,P,Q, k,l,GammaTheta);

% Concurrent Learning part
N = length(QQ);
thetaH1 = [thetaH zeros(p,n)];

GSIGMAWaH = reshape(GSIGMA*WaH,L,N);

SIGPTH=reshape(permute(reshape(reshape((SIGP*thetaH1')',L,p,N),p,L,N),[2,1,3]),L,p*N);

SIGPGDTH=reshape(permute(reshape(reshape((SIGPGD*thetaH')',L,p,N),p,L,N),[2,1,3]),L,p*N);

rc = QQ + (1/4)*(GSIGMAWaH'*WaH);  % Cost function in the meshgrid （列堆叠：N×1）
Omegac = (SIGPF1 + SIGPTH*PHI - SIGPGDTH*PHID - (1/2)*(GSIGMAWaH));

deltac = Omegac'*WcH + rc;  %（N×1）
normm = 1./(1+nu*sum(Omegac.^2,1));  %（1×N）

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
r = e'*QT*e + mu'*R*mu;  %Cost Function
delta = WcH'*omega + r;
rho = (1+nu*(omega'*omega));
WcHD = -etac1*Gamma*omega*delta/rho +clWc;
GammaD = reshape(beta*Gamma-etac1*Gamma*(omega*omega'/rho^2)*Gamma +clGamma,L*L,1);
WaHD = etaa1*(WcH-WaH)-etaa2*WaH ...
    +etac1*Gsigma'*WaH*omega'*WcH/(4*rho) +clWa;
WaHD = proj_rectangle(WaH, WaHD, -30, 30, 1);

xD = f_fun(x)+g*(u);

ydot=[xD;hd;WcHD;WaHD;GammaD;
    thetaHD;dxf;dphif;dvarphif;reshape(dP,p*p,1);reshape(dQ,p*n,1);r];
end


function a = f_fun(x)
a = [x(2); -((0.8+0.2*exp(-100*abs(x(2))))*tanh(10*x(2))+x(2))-x(1)];
end

function dv=hd_fun(v)
    dv=[v(2); -v(1)+(1-v(1)^2)*v(2)];
end


function [QQ,SIGPF1,SIGP,SIGPGD,PHI,PHID,GSIGMA,index]=create_CL_term(n,L,p,Q,R)
Ec = linspace(-2,2,3);  %离线轨迹点范围
N = length(Ec)^(2*n);   %离线轨迹点数量（考虑跟踪信号^(2*n)）
index=1;
QQ=zeros(N,1);
SIGPF1=zeros(L,N);
% SIGPF1m=zeros(L,N);
SIGP=zeros(N*L,2*n);
SIGPGD=zeros(N*L,n);
PHI=zeros(N*p,N);
PHID=zeros(N*p,N);
GSIGMA=zeros(N*L,L);
for i=1:length(Ec)
    for ii=1:length(Ec)
        for iii=1:length(Ec)
            for iiii=1:length(Ec)
                xdc = [Ec(iii) Ec(iiii)]';  %期望轨迹离线点
                hdc = hd_fun(xdc);

                ec=[Ec(i) Ec(ii)]';
                xc=ec+xdc;
                Zc=[ec;xdc];
                
                gc = [0; 1];  %修改时注意
                gdc= [0; 1];
                gplusdc = [0 1];  %伪逆
                Gc = [gc;zeros(size(gc))];  %Gc=[gc;0n×m]

                % Fc=[f_fun(xc)-gc*gplusdc*f_fun(xdc);zeros(n,1)];
                F1c = [-hdc+gc*gplusdc*hdc;hdc];  %已知部分定义为F1

                % F1cm=F1c+Fc;
                
                phic=SIDBasis(xc);
                phidc=SIDBasis(xdc);
                sig_pc=Basis(Zc);  % Actor-critic NN Basis ∈R^{L×2n}

                QQ(index)=ec'*Q*ec;  % QQ∈R^{N}
                SIGPF1(:,index)=sig_pc*F1c;
                % SIGPF1m(:,index)=sig_pc*F1cm;  % SIGPF1m∈R^{L×N}
                SIGP(index*L-(L-1):index*L,:)=sig_pc;  % SIGP∈R^{NL×2n}
                SIGPGD(index*L-(L-1):index*L,:)=sig_pc*[gc*gplusdc;zeros(n)];
                GSIGMA(index*L-(L-1):index*L,:)=sig_pc*Gc*(R\Gc')*sig_pc';  % GSIGMA∈R^{NL×L}
                PHI(index*p-(p-1):index*p,index)=phic;
                PHID(index*p-(p-1):index*p,index)=phidc;
                index = index+1;
            end
        end
    end
end
end