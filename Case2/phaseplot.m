function status = phaseplot(t, y, flag, Env)
    persistent h1 h2; % 持久变量保存图形句柄
    status = 0;   % 返回0表示继续解算
    
    switch flag
        case 'init'
            % 初始化图形
            figure(10);
            h1 = plot(y(3), y(4), 'b-'); % 初始线段：期望轨迹
            hold on;
            h2 = plot(y(1), y(2), 'r-'); % 初始线段：状态

            % plot(y(1), y(2), 'ro');     % 标记初始点
            % rectangle('Position',[-0.8 -1.25 1.4 2.5],'LineStyle','-.')
            c1=Env.c1; r1a=Env.r1a; r1b=Env.r1b;
            theta=0:0.01:2*pi;
            circle1 = [c1(1)+r1a*cos(theta); c1(2)+r1b*sin(theta)];
            
            plot(circle1(1,:),circle1(2,:),'k-.')

            xlabel('x_1');
            ylabel('x_2');
            title('二维相图');
            % legend
            grid on;
            % axis([-1.5 1.5 -1.5 1.5]);
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
            drawnow limitrate; % 实时更新图形
        case 'done'
            % 清理持久变量
            clear h1 h2;
        otherwise
    end
end