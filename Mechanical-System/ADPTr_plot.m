
%% Weights
figure(3); clf
subplot(311)
WcH = y(:,2*n+1:2*n+L);
plot(t,WcH,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex');
ylabel('$\hat{W}_c(t)$','Interpreter','latex','FontSize',12);
% title('Estimation of Critic NN');
grid on;
legend('$\hat{W}_{c1}$','$\hat{W}_{c2}$','$\hat{W}_{c3}$','$\hat{W}_{c4}$', ...
    '$\hat{W}_{c5}$','$\hat{W}_{c6}$','$\hat{W}_{c7}$','interpreter','latex',...
    'numcolumns',2);
xlim([0 40])
ylim([-0.2 5])

subplot(312)
WaH = y(:,2*n+L+1:2*n+2*L);
plot(t,WaH,'LineWidth',1)
xlabel('$t~({\rm s})$','Interpreter','latex');
ylabel('$\hat{W}_a(t)$','Interpreter','latex','FontSize',12);
% title('Estimation of Actor NN');
grid on;
legend('$\hat{W}_{a1}$','$\hat{W}_{a2}$','$\hat{W}_{a3}$','$\hat{W}_{a4}$', ...
    '$\hat{W}_{a5}$','$\hat{W}_{a6}$','$\hat{W}_{a7}$','interpreter','latex',...
    'numcolumns',2);
xlim([0 40])
ylim([-0.2 5])

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
    'interpreter','latex', 'numcolumns',2);
xlim([0 40])
y_range = ylim;             % Get the current y-axis range
ylim([y_range(1) 1.0001])
% set(gca,'ytick',[-1 0 1]);

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
% Get axes limits and ticks
ax = gca;
xlims = ax.XLim;
ylims = ax.YLim;
xticks = ax.XTick;
yticks = ax.YTick;

% Manually draw grid lines (placed on top)
hold on;  
dis=0.1; % boundary distance % 2 end-1 do not draw around the four edges
for i = 2:length(xticks)-1
    plot([xticks(i), xticks(i)], [ylims(1)+dis, ylims(2)-dis], 'color','[0.15 0.15 0.15 0.15]','HandleVisibility','off'); % vertical line
end
for i = 2:length(yticks)-1 % notice
    plot([xlims(1)+dis, xlims(2)-dis], [yticks(i), yticks(i)], 'color','[0.15 0.15 0.15 0.15]','HandleVisibility','off'); % horizontal line
end

plot(circle1(1,:),circle1(2,:), 'k-.', 'LineWidth',1,'HandleVisibility','off')  % boundary
plot(circle2(1,:),circle2(2,:), 'k-.', 'LineWidth',1,'HandleVisibility','off')  % boundary
% hold on
plot(x_plot(:,1),x_plot(:,2),'r-','LineWidth',1)
plot(xd_plot(:,1),xd_plot(:,2),'b--','LineWidth',1)

% Find the state positions corresponding to time points
[~, idx1] = min(abs(t - 1.26));
[~, idx2] = min(abs(t - 3.42));
% plot(x_plot(idx1,1), x_plot(idx1,2), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Intervention 1');
plot(x_plot(idx1,1), x_plot(idx1,2), 'x', 'MarkerSize',8, 'LineWidth', 1.5,'Color',[0.3010 0.7450 0.9330]);
plot(x_plot(idx2,1), x_plot(idx2,2), 'x', 'MarkerSize',8, 'LineWidth', 1.5,'Color',[0.4660 0.6740 0.1880]);
% Add annotation near the point (adjust offset)
text(x_plot(idx1,1)+0.1, x_plot(idx1,2)+0.1, '$t=1.26$s', 'Interpreter', 'latex', 'Color', 'k');
text(x_plot(idx2,1)-0.9, x_plot(idx2,2)-0.2, '$t=3.42$s', 'Interpreter', 'latex', 'Color', 'k');
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
B=cell2mat(Bcell);  usafe = cell2mat(Usafe);
figure(9); clf
subplot(211)
plot(t,B,'r','LineWidth',1)

idx_interval = find(t >= 3 & t <= 4);  % Find the index within the interval [t1, t2]
[B_max, local_idx] = max(B(idx_interval));  % Find the maximum value and local index
t_max = t(idx_interval(local_idx));  % The global timestamp

xline(1.26, '--', 'LineWidth',1,'Color',[0.3010 0.7450 0.9330],'Alpha',1);
xline(3.42, '--', 'LineWidth',1,'Color',[0.4660 0.6740 0.1880],'Alpha',1);
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
ylabel('$B(x,v)$','Interpreter','latex','FontSize',12);
grid on
xlim([0 4])
set(gca,'xtick',[0 1 1.26 2 3 3.42 4]);
subplot(212)
% plot(t,sqrt(sum(e.^2,2)),'b','LineWidth',1)
plot(t,usafe,'b','LineWidth',1)
xline(1.26, '--', 'LineWidth',1,'Color',[0.3010 0.7450 0.9330],'Alpha',1);
xline(3.42, '--', 'LineWidth',1,'Color',[0.4660 0.6740 0.1880],'Alpha',1);
xlabel('$t~({\rm s})$','Interpreter','latex','FontSize',12);
% ylabel('$\Vert e(t)\Vert$','Interpreter','latex','FontSize',12);
ylabel('$u_s(x,v)$','Interpreter','latex','FontSize',12);
xlim([0 4])
set(gca,'xtick',[0 1 1.26 2 3 3.42 4]);
grid on;
%% savefigures
% fig=figure(4);
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% exportgraphics(fig, 'fig_state.eps','ContentType','vector')