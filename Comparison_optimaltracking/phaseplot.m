function status = phaseplot(t, y, flag, Env)
    persistent h1 h2; % persistent variables to store graphics handles
    status = 0;
    
    switch flag
        case 'init'
            % Initialize figure
            figure(10);
            h1 = plot(y(3), y(4), 'b-'); % desired trajectory
            hold on;
            h2 = plot(y(1), y(2), 'r-'); % state

            % plot(y(1), y(2), 'ro');     % initial point
            c1=Env.c1; r1=Env.r1; c2=Env.c2; r2=Env.r2; 
            theta=0:0.01:2*pi;
            circle1 = [c1(1)+r1*cos(theta); c1(2)+r1*sin(theta)];
            circle2 = [c2(1)+r2*cos(theta); c2(2)+r2*sin(theta)];
            
            plot(circle1(1,:),circle1(2,:),'k-.')
            plot(circle2(1,:),circle2(2,:),'k-.')

            xlabel('x_1');
            ylabel('x_2');
            title('Phase Portrait');
            % legend
            % axis([-1.5 1.5 -1.5 1.5]);
        case ''
            % Update data points
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
            drawnow limitrate;
        case 'done'
            % Clear persistent variables
            clear h1 h2;
        otherwise
    end
end