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

kappa=0.1;
Gamma=10;
%% System ID
k=0.01; % 滤波器常数
l=0.1;
GammaTheta=10;
%% Initial Conditions stack
[~,L]=Basis(zeros(2*n,1));
xd0 = [2 2]';
x0 = [0 1.3]';
WcH0 = 1*ones(L,1);

[~,p]=SIDBasis(x0);
xf0 = x0;
phif0 = zeros(p,1);
varphif0 = zeros(n,1);
thetaH0 = 0*ones(p*n,1);
z0 = [x0;xd0;WcH0;zeros(L*L,1);zeros(L,1);thetaH0;xf0;phif0;varphif0;zeros(p*p,1);zeros(p*n,1);0];
%% Simulation
options = odeset('OutputFcn',@(t,y,flag)phaseplot(t,y,flag,Env),'OutputSel',[1 2 3 4]);
[t,y] = ode45(@(t,y) closedLoopDynamics(t,y,n,L,p,Q,R, ...
    kappa,Gamma,k,l,GammaTheta),tspan,z0,options);

% for k=1:length(t)
%     [~,ra(k)] = closedLoopDynamics(t(k),y(k,:)',n,L,p,Q,R,etac1,etac2,etaa1,etaa2,beta,nu, ...
%     QQ,SIGPF1,SIGP,SIGPGD,PHI,PHID,GSIGMA, k,l,GammaTheta,M,kp,kTheta, cb);
% end
[~,RA,MEP] = cellfun(@(t,y) closedLoopDynamics(t,y.',n,L,p,Q,R, ...
    kappa,Gamma,k,l,GammaTheta), num2cell(t), num2cell(y,2), 'uni',0);
%% Plot
x_plot = y(1:2001,1:2);
xd_plot = y(1:2001,3:4);
figure(1);clf
% subplot(211)
% plot(t,xd_plot(:,1),t,x_plot(:,1), 'LineWidth',1);
% subplot(212)
% plot(t,xd_plot(:,2),t,x_plot(:,2), 'LineWidth',1);
c1=Env.c1; r1=Env.r1; c2=Env.c2; r2=Env.r2; 
theta=0:0.01:2*pi;
circle1 = [c1(1)+r1*cos(theta); c1(2)+r1*sin(theta)];
circle2 = [c2(1)+r2*cos(theta); c2(2)+r2*sin(theta)];

pg1 = polyshape(circle1(1,:),circle1(2,:));
pg2 = polyshape(circle2(1,:),circle2(2,:));
pg_final = subtract(pg1,pg2);

plot(pg_final,'FaceColor',[0.7 0.9 0.7],'FaceAlpha',0.3, 'EdgeColor','none','HandleVisibility','off')

% axis([-5 4 -3.5 4.5])
% 获取坐标轴范围和刻度
ax = gca;
xlims = ax.XLim;
ylims = ax.YLim;
xticks = ax.XTick;
yticks = ax.YTick;

% 手动绘制网格线（置于顶层）
hold on;  
dis=0.1; %边界距离 % 2 end-1不画四周
for i = 2:length(xticks)-1
    plot([xticks(i), xticks(i)], [ylims(1)+dis, ylims(2)-dis], 'color','[0.15 0.15 0.15 0.15]','HandleVisibility','off'); % 垂直线
end
for i = 2:length(yticks)-1
    plot([xlims(1)+dis, xlims(2)-dis], [yticks(i), yticks(i)], 'color','[0.15 0.15 0.15 0.15]','HandleVisibility','off'); % 水平线
end

plot(circle1(1,:),circle1(2,:), 'k-.', 'LineWidth',1,'HandleVisibility','off')  % 边缘
plot(circle2(1,:),circle2(2,:), 'k-.', 'LineWidth',1,'HandleVisibility','off')  % 边缘
% hold on
plot(x_plot(:,1),x_plot(:,2),'r-','LineWidth',1)
plot(xd_plot(:,1),xd_plot(:,2),'b--','LineWidth',1)
xlabel('$x_1$','Interpreter','latex','FontSize',12);
ylabel('$x_2$','Interpreter','latex','FontSize',12);
% grid on
box on
legend('$x$','$x_d$','Interpreter','latex','fontsize',11)

figure(2); clf
WcH = y(:,2*n+1:2*n+L);
plot(t,WcH,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\hat{W}_c$','Interpreter','latex','FontSize',12);
title('Estimation of Critic NN');
grid on;

figure(3); clf
ra = cell2mat(RA);
plot(t,ra,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
title('Satisfaction of Assumption');
grid on;

figure(4); clf
theta = y(:,2*n+2*L+L*L+1:2*n+2*L+L*L+p*n);
plot(t,theta,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\hat{\theta}$','Interpreter','latex','FontSize',12);
title('Estimation of uncertain parameters');
grid on;

figure(5); clf
meP = cell2mat(MEP);
plot(t,meP,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
title('Satisfaction of PE Assumption');
grid on;

figure(6); clf
cost=y(:,end);
plot(t,cost,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\int_{0}^{t}r(\zeta(\tau),\mu(\tau)){\rm d}\tau$','Interpreter','latex','FontSize',12);
title('Instantaneous cost');
grid on;

% fig=figure(1);
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% exportgraphics(fig, 'fig_compare1.eps','ContentType','vector')

%--------------------------------------------------------------------------
function [ydot,ADPEig,meP]=closedLoopDynamics(t,y,n,L,p,QT,R, ...
    kappa,Gamma,k,l,GammaTheta)
x=y(1:n); xd=y(n+1:2*n);
WcH = y(2*n+1:2*n+L);
M = reshape(y(2*n+L+1:2*n+L+L*L),L,L);
V = y(2*n+L+L*L+1:2*n+2*L+L*L);
thetaH=reshape(y(2*n+2*L+L*L+1:2*n+2*L+L*L+p*n),p,n);
xf=y(2*n+2*L+L*L+p*n+1:3*n+2*L+L*L+p*n); phif=y(3*n+2*L+L*L+p*n+1:3*n+2*L+L*L+p*n+p); 
varphif=y(3*n+2*L+L*L+p*n+p+1:4*n+2*L+L*L+p*n+p);
P=reshape(y(4*n+2*L+L*L+p*n+p+1:4*n+2*L+L*L+p*n+p+p*p),p,p); 
Q=reshape(y(4*n+2*L+L*L+p*n+p+p*p+1:4*n+2*L+L*L+2*p*n+p+p*p),p,n); 

hd=hd_fun(xd);

e = x-xd;
Zeta = [e;xd];

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
mu = -0.5*(R\G')*sig_p'*WcH;

% Controller
u = udhat+mu;

% if t<=10
%     noise=sum(0.2*sin([1 3 7 11 13 15]*t));
% else
%     noise=0;
% end
% u = u+noise;

% identifier
[thetaHD,dxf,dphif,dvarphif,dP,dQ,meP] = identifier(thetaH,u,t,x, ...
    xf,phif,varphif,P,Q, k,l,GammaTheta);

% ADP
omega = sig_p*(FTH+F1+G*mu);
r = e'*QT*e + mu'*R*mu;  %Cost Function

varpi=omega/(omega'*omega+1);  %归一化
% 辅助矩阵
dM=-kappa*M+varpi*varpi';
dV=-kappa*V+varpi*r/(omega'*omega+1);

WcHD=-Gamma*(M*WcH+V);

ADPEig = min(eig(M));

xD = f_fun(x)+g*(u);

ydot=[xD;hd;WcHD;reshape(dM,L*L,1);dV;
    thetaHD;dxf;dphif;dvarphif;reshape(dP,p*p,1);reshape(dQ,p*n,1);r];
end


function a = f_fun(x)
a = [x(2); -((0.8+0.2*exp(-100*abs(x(2))))*tanh(10*x(2))+x(2))-x(1)];
end

function dv=hd_fun(v)
    dv=[v(2); -v(1)+(1-v(1)^2)*v(2)];
end