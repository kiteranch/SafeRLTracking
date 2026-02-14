
%% WcH
figure(1); clf
WcH = y(:,2*n+1:2*n+L);
plot(t,WcH,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\hat{W}_c$','Interpreter','latex','FontSize',12);
title('Estimation of Critic NN');
grid on;

%% WaH
figure(2); clf
WaH = y(:,2*n+L+1:2*n+2*L);
plot(t,WaH,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\hat{W}_a$','Interpreter','latex','FontSize',12);
title('Estimation of Actor NN');
grid on;

%% Theta
figure(3); clf
theta = y(:,2*n+2*L+L*L+1:2*n+2*L+L*L+p*n);
% theta=theta(:,12);
plot(t,theta,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\hat{\theta}$','Interpreter','latex','FontSize',12);
title('Estimation of uncertain parameters');
grid on;
legend;

%% State
x_plot = y(1:2001,1:2);
xd_plot = y(1:2001,3:4);
figure(4);clf
% subplot(211)
% plot(t,xd_plot(:,1),t,x_plot(:,1), 'LineWidth',1);
% subplot(212)
% plot(t,xd_plot(:,2),t,x_plot(:,2), 'LineWidth',1);

c1=Env.c1; r1a=Env.r1a; r1b=Env.r1b;
theta=0:0.01:2*pi;
circle1 = [c1(1)+r1a*cos(theta); c1(2)+r1b*sin(theta)];

pg1 = polyshape(circle1(1,:),circle1(2,:));

plot(pg1,'FaceColor',[243, 135, 183] / 255,'FaceAlpha',0.5,'EdgeColor','k','EdgeAlpha',0.5, ...
    'LineWidth',1,'LineStyle','-.')

axis([-4 4 -4 4])
% % 获取坐标轴范围和刻度
% ax = gca;
% xlims = ax.XLim;
% ylims = ax.YLim;
% xticks = ax.XTick;
% yticks = ax.YTick;
% 
% % 手动绘制网格线（置于顶层）
% hold on;  
% dis=0.1; %边界距离 % 2 end-1不画四周
% for i = 2:length(xticks)-1
%     plot([xticks(i), xticks(i)], [ylims(1)+dis, ylims(2)-dis], 'color','[0.15 0.15 0.15 0.15]','HandleVisibility','off'); % 垂直线
% end
% for i = 1:length(yticks) %注意变更
%     plot([xlims(1)+dis, xlims(2)-dis], [yticks(i), yticks(i)], 'color','[0.15 0.15 0.15 0.15]','HandleVisibility','off'); % 水平线
% end

% plot(circle1(1,:),circle1(2,:), 'k-.', 'LineWidth',1,'DisplayName','Safety boundary')  % 边缘
hold on
plot(x_plot(:,1),x_plot(:,2),'r-','LineWidth',1)
plot(xd_plot(:,1),xd_plot(:,2),'b--','LineWidth',1)
xlabel('$x_1$','Interpreter','latex','FontSize',12);
ylabel('$x_2$','Interpreter','latex','FontSize',12);
grid on
box on
legend('Obstacle','$x$','$x_d$','Interpreter','latex','fontsize',11)
%% Error
figure(5);clf
e = y(:,1:2)-y(:,3:4);
plot(t,e(:,1),t,e(:,2),'LineWidth',1);
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$e$','Interpreter','latex','FontSize',12);
grid on
title('Tracking Error')
%% Control
% u = cell2mat(U); mu = cell2mat(Mu); udhat = cell2mat(Udhat); usafe = cell2mat(Usafe);
% figure(6); clf
% % plot(t,u,'LineWidth',1)
% plot(t,mu,'r',t,udhat,'g',t,usafe,'b--',t,u,'k--','LineWidth',1)
% xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
% ylabel('$u(t)$','Interpreter','latex','FontSize',12);
% grid on
% title('Control')
% xlim([0 20])
% legend('$\hat{\mu}$','$\hat{u}_d$','$u_s$','$u=\hat{\mu}+\hat{u}_d+u_s$','interpreter','latex','fontsize',12)
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
fig2=figure(91); clf
plot(t,B,'r','LineWidth',1)
xlabel('$t$ (s)','Interpreter','latex','FontSize',12);
ylabel('$B(x,x_d)$','Interpreter','latex','FontSize',12);
grid on
xlim([0 20])
set(fig2,'Position', [573 444 560 210])

fig1=figure(9); clf
norme=sqrt(sum(e.^2,2));
plot(t,norme,'b','LineWidth',1);
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\Vert e(t)\Vert$','Interpreter','latex','FontSize',12);
xlim([0 20])
grid on;
set(fig1,'Position', [573 444 560 210])
%% Cost
figure(11); clf
cost=y(:,end);
plot(t,cost,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$\int_{0}^{t}r(\zeta(\tau),\mu(\tau)){\rm d}\tau$','Interpreter','latex','FontSize',12);
title('Instantaneous cost');
grid on;
%% savefigures
% fig=figure(9);
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% exportgraphics(fig, 'obst_error.eps')