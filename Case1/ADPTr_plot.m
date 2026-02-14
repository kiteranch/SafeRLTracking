
%% WcH
figure(1); clf
WcH = y(:,2*n+1:2*n+L);
plot(t,WcH,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\hat{W}_c(t)$','Interpreter','latex','FontSize',12);
% title('Estimation of Critic NN');
grid on;
legend('$\hat{W}_{c1}$','$\hat{W}_{c2}$','$\hat{W}_{c3}$','$\hat{W}_{c4}$', ...
    '$\hat{W}_{c5}$','$\hat{W}_{c6}$','$\hat{W}_{c7}$','interpreter','latex','fontsize',11, ...
    'numcolumns',2,'Position',[0.617280487787156,0.704308243403359-0.13,0.26843379792713,0.195691756596641]);
% axis([0 40 -4 14])

%% WaH
figure(2); clf
WaH = y(:,2*n+L+1:2*n+2*L);
plot(t,WaH,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\hat{W}_a(t)$','Interpreter','latex','FontSize',12);
% title('Estimation of Actor NN');
grid on;
legend('$\hat{W}_{a1}$','$\hat{W}_{a2}$','$\hat{W}_{a3}$','$\hat{W}_{a4}$', ...
    '$\hat{W}_{a5}$','$\hat{W}_{a6}$','$\hat{W}_{a7}$','interpreter','latex','fontsize',11, ...
    'numcolumns',2,'Position',[0.613766579400925,0.704308243403359-0.13,0.27194770631336,0.195691756596641]);
% axis([0 40 -4 14])

%% Theta
figure(3); clf
theta = y(:,2*n+2*L+L*L+1:2*n+2*L+L*L+p*n);
% plot(t,theta,'LineWidth',1)
styles={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],...
    [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};
hold on
for j=1:p
    plot(t,theta(:,j),'color',styles{j},'LineStyle','--','LineWidth',1)
end
for j=p+1:2*p
    plot(t,theta(:,j),'color',styles{j-p},'LineStyle','-','LineWidth',1)
end
% end
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\hat{\theta}(t)$','Interpreter','latex','FontSize',12);
% title('Estimation of uncertain parameters');
grid on;
box on
legend('$\hat{\theta}_{11}$','$\hat{\theta}_{21}$','$\hat{\theta}_{31}$','$\hat{\theta}_{41}$','$\hat{\theta}_{51}$', ...
    '$\hat{\theta}_{12}$','$\hat{\theta}_{22}$','$\hat{\theta}_{32}$','$\hat{\theta}_{42}$','$\hat{\theta}_{52}$', ...
    'interpreter','latex', 'numcolumns',2);
% xlim([0 40])

%%
figure(6); clf
subplot(311)
WcH = y(:,2*n+1:2*n+L);
plot(t,WcH,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex');
ylabel('$\hat{W}_c(t)$','Interpreter','latex','FontSize',12);
% title('Estimation of Critic NN');
grid on;
legend('$\hat{W}_{c1}$','$\hat{W}_{c2}$','$\hat{W}_{c3}$','$\hat{W}_{c4}$', ...
    '$\hat{W}_{c5}$','$\hat{W}_{c6}$','$\hat{W}_{c7}$','interpreter','latex',...
    'numcolumns',4,'Orientation', 'horizontal','fontsize',9.5, ...
    'Position',[0.372380592709496+0.015,0.820751610253447-0.005,0.51333369300479,0.089562777867393]);
xlim([0 40])

subplot(312)
WaH = y(:,2*n+L+1:2*n+2*L);
plot(t,WaH,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex');
ylabel('$\hat{W}_a(t)$','Interpreter','latex','FontSize',12);
% title('Estimation of Actor NN');
grid on;
legend('$\hat{W}_{a1}$','$\hat{W}_{a2}$','$\hat{W}_{a3}$','$\hat{W}_{a4}$', ...
    '$\hat{W}_{a5}$','$\hat{W}_{a6}$','$\hat{W}_{a7}$','interpreter','latex',...
    'numcolumns',4,'Orientation', 'horizontal','fontsize',9.5, ...
    'Position',[0.366311118716285+0.015,0.520751610253447-0.005,0.519403166998,0.089562777867393]);
xlim([0 40])

subplot(313);
theta = y(:,2*n+2*L+L*L+1:2*n+2*L+L*L+p*n);
% plot(t,theta,'LineWidth',1)
styles={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],...
    [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};
hold on
for j=1:p
    plot(t,theta(:,j),'color',styles{j},'LineStyle','--','LineWidth',1)
end
for j=p+1:2*p
    plot(t,theta(:,j),'color',styles{j-p},'LineStyle','-','LineWidth',1)
end
% end
xlabel('$t~({\rm s})$','Interpreter','latex');
ylabel('$\hat{W}_{\theta}(t)$','Interpreter','latex','FontSize',12);
% title('Estimation of uncertain parameters');
grid on;
box on
legend('$\hat{W}_{\theta 1}$','$\hat{W}_{\theta 2}$','$\hat{W}_{\theta 3}$','$\hat{W}_{\theta 4}$','$\hat{W}_{\theta 5}$', ...
    '$\hat{W}_{\theta 6}$','$\hat{W}_{\theta 7}$','$\hat{W}_{\theta 8}$','$\hat{W}_{\theta 9}$','$\hat{W}_{\theta 10}$', ...
    'interpreter','latex', 'numcolumns',5,'Orientation', 'horizontal','fontsize',6, ...
    'Position',[0.304464367457799+0.015,0.249801050667633-0.08,0.581249918256487,0.060513337453206]);
xlim([0 40])

%% State
x_plot = y(1:2001,1:2);
xd_plot = y(1:2001,3:4);
figure(4);clf
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
% hold on
plot(x_plot(:,1),x_plot(:,2),'r-','LineWidth',1)
plot(xd_plot(:,1),xd_plot(:,2),'b--','LineWidth',1)
xlabel('$x_1$','Interpreter','latex','FontSize',12);
ylabel('$x_2$','Interpreter','latex','FontSize',12);
% grid on
box on
legend('$x$','$x_d$','Interpreter','latex','FontSize',11)
%% Error
figure(5);clf
e = y(:,1:2) - y(:,3:4);
plot(t,e,'LineWidth',1);
% plot(t,e(:,1),t,e(:,2),'LineWidth',1);
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$e$','Interpreter','latex','FontSize',12);
grid on
title('Tracking Error')
%% Control
% u = cell2mat(U); mu = cell2mat(Mu); udhat = cell2mat(Udhat); usafe = cell2mat(Usafe);
% figure(6); clf
% % plot(t,u,'LineWidth',1)
% plot(t,udhat,'r',t,mu,'g',t,usafe,'b--',t,u,'k--','LineWidth',1)
% xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
% ylabel('$u(t)$','Interpreter','latex','FontSize',12);
% grid on
% % title('Control')
% xlim([0 20])
% legend('$\hat{u}_d$','$\hat{u}_e$','$u_s$','$u=\hat{u}_d+\hat{u}_e+u_s$','interpreter','latex','fontsize',12)
%% ADPRank
figure(7); clf
ra = cell2mat(RA);
plot(t,ra,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
title('Satisfaction of Assumption');
grid on;
%% SIDRank
figure(8); clf
meP = cell2mat(MEP);
plot(t,meP,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
title('Satisfaction of PE Assumption');
grid on;
%% CBF
B=cell2mat(Bcell); h=cell2mat(hcell);
figure(9); clf
subplot(211)
plot(t,B,'r','LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$B(x,v)$','Interpreter','latex','FontSize',12);
grid on
xlim([0 20])
subplot(212)
plot(t,sqrt(sum(e.^2,2)),'b','LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\Vert e(t)\Vert$','Interpreter','latex','FontSize',12);
xlim([0 20])
grid on;
%% Cost
figure(11); clf
cost=y(:,end);
plot(t,cost,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\int_{0}^{t}r(\zeta(\tau),\mu(\tau)){\rm d}\tau$','Interpreter','latex','FontSize',12);
title('Instantaneous cost');
grid on;
%% savefigures
% fig=figure(6);
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% exportgraphics(fig, 'fig_weight.eps','ContentType','vector')