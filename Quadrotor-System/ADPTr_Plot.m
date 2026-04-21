
x_plot = rad2deg(y(:,1:6));
v = rad2deg(y(:,7:8));
xd_plot = [v(:,1) v(:,1) v(:,1) v(:,2) v(:,2) v(:,2)];
figure(1); clf
subplot(311)
plot(t,x_plot(:,1),'r-','LineWidth',1)
hold on
plot(t,xd_plot(:,1),'b--','LineWidth',1)
yline(-Consdeg,'k-.','LineWidth',1)
yline(Consdeg,'k-.','LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',11);
ylabel('$\phi~({\rm ^\circ})$','Interpreter','latex','FontSize',11);
grid on
% legend('$x$','$x_d$','Interpreter','latex','Location','best')
ylim([-50 55])

subplot(312)
plot(t,x_plot(:,2),'r-','LineWidth',1)
hold on
plot(t,xd_plot(:,2),'b--','LineWidth',1)
yline(-Consdeg,'k-.','LineWidth',1)
yline(Consdeg,'k-.','LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',11);
ylabel('$\theta~({\rm ^\circ})$','Interpreter','latex','FontSize',11);
grid on
ylim([-50 55])

subplot(313)
plot(t,x_plot(:,3),'r-','LineWidth',1)
hold on
plot(t,xd_plot(:,3),'b--','LineWidth',1)
yline(-Consdeg,'k-.','LineWidth',1)
yline(Consdeg,'k-.','LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',11);
ylabel('$\psi~({\rm ^\circ})$','Interpreter','latex','FontSize',11);
grid on
ylim([-50 55])

figure(2); clf
subplot(311)
plot(t,x_plot(:,4),'r-','LineWidth',1)
hold on
plot(t,xd_plot(:,4),'b--','LineWidth',1)
yline(-Consdeg,'k-.','LineWidth',1)
yline(Consdeg,'k-.','LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',11);
ylabel('$\dot{\phi}~({\rm ^\circ/s})$','Interpreter','latex','FontSize',11);
grid on
y_range = ylim;             % Get the current y-axis range
% disp(y_range);
subplot(312)
plot(t,x_plot(:,5),'r-','LineWidth',1)
hold on
plot(t,xd_plot(:,5),'b--','LineWidth',1)
yline(-Consdeg,'k-.','LineWidth',1)
yline(Consdeg,'k-.','LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',11);
ylabel('$\dot{\theta}~({\rm ^\circ/s})$','Interpreter','latex','FontSize',11);
grid on
subplot(313)
plot(t,x_plot(:,6),'r-','LineWidth',1)
hold on
plot(t,xd_plot(:,6),'b--','LineWidth',1)
yline(-Consdeg,'k-.','LineWidth',1)
yline(Consdeg,'k-.','LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',11);
ylabel('$\dot{\psi}~({\rm ^\circ/s})$','Interpreter','latex','FontSize',11);
grid on

% figure(3); clf
% subplot(211)
% plot(t,y(:,n+nv+1:n+nv+L))
% subplot(212)
% plot(t,y(:,n+nv+L+1:n+nv+2*L))
% 
% figure(4); clf
% plot(t,y(:,n+nv+2*L+L*L+1:n+nv+2*L+L*L+p*n))
% % legend
% 
% u = reshape(cell2mat(U),3,[]);
% figure(5); clf
% plot(t,u)
% 
% figure(6); clf
% plot(t,cell2mat(ADPRank))
% figure(7); clf
% plot(t,cell2mat(MEP))


%% CBF
u = cell2mat(U); mu = cell2mat(Mu); udhat = cell2mat(Udhat); usafe = reshape(cell2mat(Usafe),3,[]);
B=cell2mat(Bcell);
figure(9); clf
subplot(211)
plot(t,B,'r','LineWidth',1)

idx_interval = find(t >= 3 & t <= 4);  % Find the index within the interval [t1, t2]
[B_max, local_idx] = max(B(idx_interval));  % Find the maximum value and local index
t_max = t(idx_interval(local_idx));  % The global timestamp

% xline(1.26, '--', 'LineWidth',1,'Color',[0.3010 0.7450 0.9330],'Alpha',1);
% xline(3.42, '--', 'LineWidth',1,'Color',[0.4660 0.6740 0.1880],'Alpha',1);
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$B(x,v)$','Interpreter','latex','FontSize',12);
grid on
xlim([0 4])
% set(gca,'xtick',[0 1 1.26 2 3 3.42 4]);
subplot(212)
% plot(t,sqrt(sum(e.^2,2)),'b','LineWidth',1)
plot(t,usafe,'b','LineWidth',1)
% xline(1.26, '--', 'LineWidth',1,'Color',[0.3010 0.7450 0.9330],'Alpha',1);
% xline(3.42, '--', 'LineWidth',1,'Color',[0.4660 0.6740 0.1880],'Alpha',1);
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
% ylabel('$\Vert e(t)\Vert$','Interpreter','latex','FontSize',12);
ylabel('$u_s(x,v)$','Interpreter','latex','FontSize',12);
xlim([0 4])
% set(gca,'xtick',[0 1 1.26 2 3 3.42 4]);
grid on;
%%
% fig=figure(1);
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% exportgraphics(fig, 'quad_state_uncons.eps','ContentType','vector')