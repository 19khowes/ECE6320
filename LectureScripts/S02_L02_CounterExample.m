function S02_L02_CounterExample
close all;
    % Setup simulation
    t = [0:.01:10];
    x0 = [1; 0];
    
    % Simulate system
    %[tvec xvec] = ode45(@f_A1, t, x0);
    %[tvec xvec] = ode45(@f_A2, t, x0);
    %[tvec xvec] = ode45(@f_unstable, t, x0);
    [tvec xvec] = ode45(@f_stable, t, x0);
    
    % Plot the system
    plot(xvec(:,1), xvec(:,2), 'linewidth', 3); hold on;
    plot(x0(1), x0(2), 'go', 'linewidth', 3);
    set(gca, 'fontsize', 20)
    axis equal
end

function xdot = f_A1(t, x)
    A1 = [-0.1 1; -10 -.1];
    xdot = A1*x;
end

function xdot = f_A2(t, x)
    A2 = [-0.1 10; -1 -.1];
    xdot = A2*x;
end

function xdot = f_unstable(t, x)
    if x(1) >= 0 && x(2) >= 0 % first quadrant
        xdot = f_A2(t,x);
    elseif x(1) >= 0 && x(2) <= 0 % second quadrant
        xdot = f_A1(t,x);
    elseif x(1) <= 0 && x(2) <= 0 % third quadrant
        xdot = f_A2(t,x);
    else                        % Fourth quadrant
        xdot = f_A1(t,x);
    end
end

function xdot = f_stable(t, x)
    if x(1) >= 0 && x(2) >= 0 % first quadrant
        xdot = f_A1(t,x);
    elseif x(1) >= 0 && x(2) <= 0 % second quadrant
        xdot = f_A2(t,x);
    elseif x(1) <= 0 && x(2) <= 0 % third quadrant
        xdot = f_A1(t,x);
    else                        % Fourth quadrant
        xdot = f_A2(t,x);
    end
end