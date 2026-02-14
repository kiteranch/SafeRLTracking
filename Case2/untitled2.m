function precise_crossing_detection()
    % 定义椭圆
    ellipse_fun = @(x1, x2) ((x1 + 2).^2)/(0.7^2) + x2.^2 - 1;
    
    % 范德波尔系统
    function dx = vdp_sys(t, x)
        dx = [x(2);
              -x(1) + (1 - x(1)^2)*x(2)];
    end
    
    % 改进的事件函数：检测符号变化
    function [value, isterminal, direction] = crossing_events(t, x)
        persistent prev_value;
        
        if isempty(prev_value)
            prev_value = ellipse_fun(x(1), x(2));
        end
        
        current_value = ellipse_fun(x(1), x(2));
        
        % 检测符号变化
        value = prev_value * current_value;
        isterminal = 0;  % 不停止
        direction = -1;  % 只检测负向穿越（符号变化）
        
        prev_value = current_value;
    end
    
    % 初始条件
    x0 = [2; 2];
    tspan = [0 25];
    
    options = odeset('Events', @crossing_events, 'RelTol', 1e-10);
    [t, x, te, xe, ie] = ode45(@vdp_sys, tspan, x0, options);
    
    % 分析结果
    analyze_crossings(t, x, te, xe, ellipse_fun);
end

function analyze_crossings(t, x, te, xe, ellipse_fun)
    fprintf('\n=== 穿越事件分析 ===\n');
    
    if isempty(te)
        fprintf('未检测到穿越事件\n');
        return;
    end
    
    % 计算每个时间点的椭圆函数值
    ellipse_vals = ellipse_fun(x(:,1), x(:,2));
    
    entry_times = [];
    exit_times = [];
    
    for i = 1:length(te)
        % 找到事件前后的索引
        idx_before = find(t < te(i), 1, 'last');
        idx_after = find(t > te(i), 1, 'first');
        
        if isempty(idx_before) || isempty(idx_after)
            continue;
        end
        
        val_before = ellipse_vals(idx_before);
        val_after = ellipse_vals(idx_after);
        
        % 判断穿越方向
        if val_before > 0 && val_after < 0
            entry_times = [entry_times; te(i)];
            fprintf('进入: t = %.6f, 位置: (%.6f, %.6f)\n', te(i), xe(i,1), xe(i,2));
        elseif val_before < 0 && val_after > 0
            exit_times = [exit_times; te(i)];
            fprintf('离开: t = %.6f, 位置: (%.6f, %.6f)\n', te(i), xe(i,2), xe(i,2));
        end
    end
    
    % 绘制详细分析图
    figure;
    
    % 轨迹和椭圆
    subplot(2,1,1);
    plot(x(:,1), x(:,2), 'b-');
    hold on;
    
    % 椭圆边界
    theta = linspace(0, 2*pi, 200);
    x1_e = -2 + 0.7*cos(theta);
    x2_e = sin(theta);
    plot(x1_e, x2_e, 'r--', 'LineWidth', 2);
    
    if ~isempty(entry_times)
        entry_points = interp1(t, x, entry_times);
        plot(entry_points(:,1), entry_points(:,2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    end
    if ~isempty(exit_times)
        exit_points = interp1(t, x, exit_times);
        plot(exit_points(:,1), exit_points(:,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end
    
    legend('轨迹', '椭圆', '进入点', '离开点');
    xlabel('x_1'); ylabel('x_2');
    title('轨迹穿越分析');
    axis equal;
    grid on;
    
    % 椭圆函数值随时间变化
    subplot(2,1,2);
    plot(t, ellipse_vals, 'b-');
    hold on;
    plot([t(1) t(end)], [0 0], 'k--');
    if ~isempty(te)
        plot(te, zeros(size(te)), 'ro', 'MarkerSize', 6);
    end
    xlabel('时间 t');
    ylabel('椭圆函数值');
    title('轨迹相对椭圆的位置');
    grid on;
end