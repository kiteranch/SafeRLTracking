function calculate_crossing_times()
    % 定义椭圆方程
    ellipse = @(x1, x2) ((x1 + 2).^2)/(0.7^2) + x2.^2 - 1;
    
    % 定义范德波尔系统
    function dx = vanderpol(t, x)
        dx = [x(2);
              -x(1) + (1 - x(1)^2)*x(2)];
    end
    
    % 定义事件函数：当轨迹穿过椭圆边界时触发
    function [value, isterminal, direction] = ellipse_events(t, x)
        value = ellipse(x(1), x(2));  % 椭圆边界方程
        isterminal = 0;  % 不停止积分，只记录事件
        direction = 0;   % 检测所有方向的穿越
    end
    
    % 初始条件（根据实际情况调整）
    x0 = [2; 2];  % 初始状态
    
    % 设置积分时间
    tspan = [0 20];
    
    % 设置ODE选项，启用事件检测
    options = odeset('Events', @ellipse_events, 'RelTol', 1e-8, 'AbsTol', 1e-10);
    
    % 求解ODE并检测事件
    [t, x, te, xe, ie] = ode45(@vanderpol, tspan, x0, options);
    
    % 计算每个事件点与椭圆边界的距离符号变化
    ellipse_values = ellipse(x(:,1), x(:,2));
    
    % 分析穿越事件
    if ~isempty(te)
        fprintf('检测到 %d 次穿越事件\n', length(te));
        
        % 区分进入和离开事件
        entry_times = [];
        exit_times = [];
        
        for i = 1:length(te)
            % 检查事件前后的符号变化
            if i == 1
                prev_sign = sign(ellipse(x(1,1), x(1,2)));
            else
                prev_time_idx = find(t < te(i), 1, 'last');
                prev_sign = sign(ellipse(x(prev_time_idx,1), x(prev_time_idx,2)));
            end
            
            current_sign = sign(ellipse(xe(i,1), xe(i,2)));
            
            % 从正到负：进入椭圆；从负到正：离开椭圆
            if prev_sign > 0 && current_sign <= 0
                entry_times = [entry_times; te(i)];
                fprintf('进入椭圆: t = %.4f, 位置: (%.4f, %.4f)\n', te(i), xe(i,1), xe(i,2));
            elseif prev_sign < 0 && current_sign >= 0
                exit_times = [exit_times; te(i)];
                fprintf('离开椭圆: t = %.4f, 位置: (%.4f, %.4f)\n', te(i), xe(i,1), xe(i,2));
            end
        end
    else
        fprintf('在时间区间内未检测到穿越事件\n');
    end
    
    % 绘制结果
    plot_results(t, x, te, xe);
end

function plot_results(t, x, te, xe)
    figure;
    
    % 子图1: 轨迹和椭圆
    subplot(2,2,[1,2]);
    plot(x(:,1), x(:,2), 'b-', 'LineWidth', 1.5);
    hold on;
    
    % 绘制椭圆
    theta = linspace(0, 2*pi, 100);
    x1_ellipse = -2 + 0.7*cos(theta);
    x2_ellipse = sin(theta);
    plot(x1_ellipse, x2_ellipse, 'r--', 'LineWidth', 2);
    
    % 标记穿越点
    if ~isempty(te)
        plot(xe(:,1), xe(:,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
    end
    
    grid on;
    xlabel('x_1');
    ylabel('x_2');
    title('范德波尔轨迹与椭圆');
    legend('轨迹', '椭圆边界', '穿越点', 'Location', 'best');
    axis equal;
    
    % 子图2: x1随时间变化
    subplot(2,2,3);
    plot(t, x(:,1), 'b-', 'LineWidth', 1.5);
    hold on;
    if ~isempty(te)
        plot(te, xe(:,1), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'red');
    end
    grid on;
    xlabel('时间 t');
    ylabel('x_1');
    title('x_1 随时间变化');
    
    % 子图3: x2随时间变化
    subplot(2,2,4);
    plot(t, x(:,2), 'b-', 'LineWidth', 1.5);
    hold on;
    if ~isempty(te)
        plot(te, xe(:,2), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'red');
    end
    grid on;
    xlabel('时间 t');
    ylabel('x_2');
    title('x_2 随时间变化');
end