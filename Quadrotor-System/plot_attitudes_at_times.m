%% 绘制 t = 0,1,2,3 秒处的四旋翼姿态（2x2 子图）
clear;clc;

%% 1. 加载仿真数据（或直接使用工作区变量）
load('sim_data.mat');   % 假设已保存

% 如果需要的时间点
target_times = [0, 1, 2, 3];
num_times = length(target_times);

% 找到每个目标时间在 t 向量中的最近索引
indices = zeros(1, num_times);
for i = 1:num_times
    [~, idx] = min(abs(t - target_times(i)));
    indices(i) = idx;
end

%% 2. 四旋翼模型定义（机体坐标系下）
arm_length = 0.4;   % 臂长
% 四个臂端点的局部坐标：前(y+), 右(x+), 后(y-), 左(x-)
arm_ends = [0,  arm_length, 0;    % 前
            arm_length, 0, 0;     % 右
            0, -arm_length, 0;    % 后
            -arm_length, 0, 0];   % 左

%% 3. 创建 2x2 图形窗口
% fig = figure('Name', 'Quadrotor Attitude at t = 0,1,2,3 s', ...
       % 'NumberTitle', 'off', 'Color', 'w', 'Position', [100, 100, 1200, 900]);
% set(gcf, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]);
figure(2); clf

for i = 1:num_times
    idx = indices(i);
    t_now = t(idx);
    
    % 获取该时刻的实际与参考欧拉角
    phi_a = phi_act(idx);
    theta_a = theta_act(idx);
    psi_a = psi_act(idx);
    
    phi_r = phi_ref(idx);
    theta_r = theta_ref(idx);
    psi_r = psi_ref(idx);
    
    % 计算旋转矩阵 (ZYX 顺序，与仿真一致)
    R_act = eul2rotm([psi_a, theta_a, phi_a], 'ZYX');
    R_ref = eul2rotm([psi_r, theta_r, phi_r], 'ZYX');
    
    % 计算臂端点的世界坐标
    act_endpoints = (R_act * arm_ends')';
    ref_endpoints = (R_ref * arm_ends')';
    
    % 创建子图
    subplot(2,2,i);
    % 手动指定每个子图的位置 [left, bottom, width, height]（归一化坐标）
    % if i == 1
    %     set(gca, 'Position', [0.10, 0.55, 0.42, 0.38]);
    % elseif i == 2
    %     set(gca, 'Position', [0.55, 0.55, 0.42, 0.38]);
    % elseif i == 3
    %     set(gca, 'Position', [0.10, 0.08, 0.42, 0.38]);
    % else % i == 4
    %     set(gca, 'Position', [0.55, 0.08, 0.42, 0.38]);
    % end
    hold on; grid on; 
    % axis equal;
    view(45, 30);   % 视角与动画一致
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title(sprintf('t = %.0f s', t_now));
    axis([-0.8 0.8 -0.8 0.8 -0.5 0.5]);
    
    % 绘制参考模型（蓝色半透明）
    for j = 1:4
        p = ref_endpoints(j, :);
        plot3([0, p(1)], [0, p(2)], [0, p(3)], 'b-', 'LineWidth', 2, 'Color', [0, 0, 1, 0.3]);
    end
    
    % 绘制实际模型（红色实体）
    for j = 1:4
        p = act_endpoints(j, :);
        plot3([0, p(1)], [0, p(2)], [0, p(3)], 'r-', 'LineWidth', 2);
        % 在臂端点绘制旋翼小球
        scatter3(p(1), p(2), p(3), 40, 'r', 'filled');
    end
    % 绘制实际中心点
    scatter3(0, 0, 0, 53, 'k', 'filled');
    
    % 添加图例（仅在第一个子图显示，避免重复）
    if i == 1
        % 创建不可见对象用于图例
        h_act_line = plot3(NaN, NaN, NaN, 'r-', 'LineWidth', 2);
        h_ref_line = plot3(NaN, NaN, NaN, 'b-', 'LineWidth', 2, 'Color', [0,0,1,0.3]);
        legend([h_act_line, h_ref_line], {'Actual', 'Reference'}, 'Location', 'best');
    end
    
    hold off;
end

% sgtitle('Quadrotor Attitude at Selected Times');  % 总标题

%%
fig=figure(1);
set(fig,'Position', [100 100 800 600])

fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
exportgraphics(fig, 'quad_1.eps','ContentType','vector')