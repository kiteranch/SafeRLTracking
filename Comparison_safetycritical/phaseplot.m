function status = phaseplot(t, y, flag, Env)
    persistent h1 h2; % 持久变量保存图形句柄
    status = 0;   % 返回0表示继续解算
    
    switch flag
        case 'init'
            % 初始化图形
            figure(10);
            h1 = plot(y(3), y(4), 'r-','LineWidth',1); % 初始线段
            hold on;
            h2 = plot(y(1), y(2), 'b-','LineWidth',1); % 初始线段

            % plot(y(1), y(2), 'ko');     % 标记初始点
            % plot(0, 0, 'go', 'MarkerFaceColor','g');     % 标记目标点
            
            % rectangle('Position', [-3 -3 6 6], 'EdgeColor','k','LineStyle','-.','LineWidth',1)
            c1=Env.c1; r1=Env.r1; c2=Env.c2; r2=Env.r2; 
            theta=0:0.01:2*pi;
            circle1 = [c1(1)+r1*cos(theta); c1(2)+r1*sin(theta)];
            circle2 = [c2(1)+r2*cos(theta); c2(2)+r2*sin(theta)];
            
            plot(circle1(1,:),circle1(2,:),'k-.')
            plot(circle2(1,:),circle2(2,:),'k-.')
            
            xlabel('x_1');
            ylabel('x_2');
            title('二维相图');
            grid on;
            % axis([-0.6 2.5 -0.6 3.5]);
        case ''
            % 更新数据点
            if ~isempty(h1)
                newX = [h1.XData, y(3)];
                newY = [h1.YData, y(4)];
                set(h1, 'XData', newX, 'YData', newY);
            end
            if ~isempty(h2)
                newX2 = [h2.XData, y(1)];
                newY2 = [h2.YData, y(2)];
                set(h2, 'XData', newX2, 'YData', newY2);
            end
            drawnow; % 实时更新图形
        case 'done'
            % 清理持久变量
            clear h1 h2;
        otherwise
    end
end