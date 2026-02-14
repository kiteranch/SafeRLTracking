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
    'numcolumns',2);
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
    'numcolumns',2);
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
    'interpreter','latex', 'numcolumns',2);
xlim([0 40])

%% savefigures
% fig=figure(6);
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% exportgraphics(fig, 'fig_weight.eps','ContentType','vector')



T_steady = 5; % 稳态开始时间（根据实际仿真调整）
steady_state_indices = t >= T_steady; % 稳态时间索引

% 计算稳态期间的最大跟踪误差范数
e_norm_steady = sqrt(sum(e(steady_state_indices,:).^2,2)); % e 为跟踪误差矩阵 [e1; e2; ...]
epsilon = max(e_norm_steady)