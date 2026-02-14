% Wang XY Automatic mechanical system
clear all
close all
clc

tspan=0:0.01:100;

Env.c1=[0;0]; Env.r1=3;
Env.c2=[0;0]; Env.r2=sqrt(0.5);

x0=[0 1.3]';
v0=[2 2]';

options=odeset('OutputFcn',@(t,y,flag) phaseplot(t,y,flag,Env),'OutputSel',1:4);
[t,y] = ode23(@(t,y) control(t,y,Env), tspan, [x0;v0], options);

[~,Wcell] = cellfun(@(t,y)control(t,y.',Env), num2cell(t), num2cell(y,2), 'uni',0);
%% plot
figure(1); clf
c1=Env.c1; r1=Env.r1; c2=Env.c2; r2=Env.r2; 
theta=0:0.01:2*pi;
circle1 = [c1(1)+r1*cos(theta); c1(2)+r1*sin(theta)];
circle2 = [c2(1)+r2*cos(theta); c2(2)+r2*sin(theta)];

pg1 = polyshape(circle1(1,:),circle1(2,:));
pg2 = polyshape(circle2(1,:),circle2(2,:));
pg_final = subtract(pg1,pg2);

plot(pg_final,'FaceColor',[0.7 0.9 0.7],'FaceAlpha',0.3, 'EdgeColor','none','HandleVisibility','off')

% ylim([-3.5 4.5])
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
for i = 2:length(yticks)-1 %注意变更
    plot([xlims(1)+dis, xlims(2)-dis], [yticks(i), yticks(i)], 'color','[0.15 0.15 0.15 0.15]','HandleVisibility','off'); % 水平线
end

plot(circle1(1,:),circle1(2,:), 'k-.', 'LineWidth',1,'HandleVisibility','off')  % 边缘
plot(circle2(1,:),circle2(2,:), 'k-.', 'LineWidth',1,'HandleVisibility','off')  % 边缘

plot(y(1:2001,1), y(1:2001,2), 'r-', 'LineWidth',1)
plot(y(1:2001,3), y(1:2001,4), 'b--', 'LineWidth',1)
% grid on
xlabel('$x_1$','Interpreter','latex','FontSize',12);
ylabel('$x_2$','Interpreter','latex','FontSize',12);
% grid on
box on
legend('$x$','$x_d$','Interpreter','latex','Fontsize',11)

figure(2); clf
W=cell2mat(Wcell);
e=y(:,1:2)-y(:,3:4);
subplot(211)
plot(t,W,'LineWidth',1)
grid on
subplot(212)
plot(t,sum(e.^2,2),'LineWidth',1)
grid on

% fig=figure(1);
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% exportgraphics(fig, 'fig_compare2.eps','ContentType','vector')
%--------------------------------------------------------------------------
function [dydt,W]=control(t,y,Env)
    x=y(1:2); v=y(3:4);

    dv=reference(v);
    
    fv=-((0.8+0.2*exp(-100*abs(x(2))))*tanh(10*x(2))+x(2))-x(1);
    
    e=x-v;
    
    % nominal parameters
    beta1=1; beta2=1;
    u0=(-fv+dv(2)-beta1*e(1)-beta2*e(2));
    
    % safety critical
    c1=Env.c1; r1=Env.r1; c2=Env.c2; r2=Env.r2; 
    h1 = r1^2-norm(x-c1)^2;
    h1_partial = -2*(x-c1)';
    
    h2 = norm(x-c2)^2-r2^2;
    h2_partial = 2*(x-c2)';

    P= [1.5 0.5;
        0.5 1];

    W=(1/2)*e'*P*e*(1/h1+1/h2);
    dWdx=e'*P*(1/h1+1/h2) + (1/2)*e'*P*e*(-h1_partial/h1^2 -h2_partial/h2^2);
    dWdv=-e'*P*(1/h1+1/h2);

    g=[0; 1];

    p=dWdx*g;
    d=dWdx*([x(2);fv]+g*u0)+dWdv*dv;
    if p==0
        delta=0
    else
        delta=-(d+sqrt(norm(d)^2+0.1*norm(p)^4))*p'/norm(p)^2;
    end

    % u=u0-5*p';
    u=u0+delta;

    dx=dyn_fcn(x,u);
    dydt=[dx;dv];
end

function dv=reference(v)
    dv=[v(2); -v(1)+(1-v(1)^2)*v(2)];
end

function dx=dyn_fcn(x,u)
    f=[x(2); -((0.8+0.2*exp(-100*abs(x(2))))*tanh(10*x(2))+x(2))-x(1)];
    g=[0; 1];
    dx=f+g*u;
end