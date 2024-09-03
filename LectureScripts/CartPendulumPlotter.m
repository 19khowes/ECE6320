classdef CartPendulumPlotter < handle
    %CartPendulum Plots the cart pendulum
    
    % Cart-Pendulum properties
    properties
        l = 0.3 % length of the pendulum
        h = 0.1 % Height of the cart
        w = 0.2 % Width of the cart
    end
    
    % Plotting variables
    properties
        % Handles for plotting
        ax = [] % Handle to the plot axes
        h_cart = [] % Handle to patch object for the body of the cart
        h_pend = [] % Handle to the bar of the pendulum
        h_mass = [] % Handle to the pendulum mass
        
        % Other plotting variables
        v_x % Nominal x coordinates of the cart
        v_y % Nominal y coordinates of the cart
        xlim_nom % Nominal limit for x axis
        xlim_act % Actual limit for the x axis
        xlim_diff % Length of plotting on x axis
    end
    
    
    methods
        function obj = CartPendulumPlotter()
            %CartPendulum Constructor which initializes the plots
            obj.initializePlot();
        end
        
        function initializePlot(obj)
            %initializePlot Creates the initial cart, shaft, and pendulum
            %mass
            obj.v_x = [-obj.w/2; obj.w/2; obj.w/2; -obj.w/2];
            obj.v_y = [obj.h; obj.h; 0; 0];
            
            % create the cart in the zero position
            c = [3, 63, 99]./255; % Metallic blue
            obj.h_cart = patch(obj.v_x, obj.v_y, c); hold on;
            
            % Create the pendulum
            p1 = [0;0]; %[0; obj.h];
            p2 = [0; -obj.l]; %[0; obj.h+obj.l];
            obj.h_pend = plot([p1(1) p2(1)], [p1(2) p2(2)], 'k', 'linewidth', 2);
            
            % Plot the mass
            obj.h_mass = plot(p2(1), p2(2), 'ro', 'linewidth', 10);
            
            % Adjust the plot
            obj.ax = gca;
            obj.xlim_nom = [-obj.w - obj.l, obj.w+obj.l];
            obj.xlim_act = obj.xlim_nom;
            obj.xlim_diff = obj.xlim_nom(2) - obj.xlim_nom(1);
            ylim = [-obj.h - obj.l, obj.h + obj.l];
            set(obj.ax, 'xlim', obj.xlim_nom);
            set(obj.ax, 'ylim', ylim);
        end
        
        function plot(obj, x, theta)
           % Update the cart position
           v_x_adjusted = obj.v_x+x;
           set(obj.h_cart, 'xdata', v_x_adjusted);
           
           % Plot the pendulum
           p1 = [x; 0]; %[x; obj.h];
           c = cos(theta);
           s = sin(theta);
           p2 = [s; -c] * obj.l + p1;
           set(obj.h_pend, 'xdata', [p1(1) p2(1)], 'ydata', [p1(2) p2(2)]);
           
           % Plot the mass
           set(obj.h_mass, 'xdata', p2(1), 'ydata', p2(2));
           
           % Set the x-axis limit
           x_max = max([v_x_adjusted; p2(1)]);
           x_min = min([v_x_adjusted; p2(1)]);
           if x_max > obj.xlim_act(2)
                obj.xlim_act = [x_max - obj.xlim_diff, x_max];
                set(obj.ax, 'xlim', obj.xlim_act);
           elseif x_min < obj.xlim_act(1)
                obj.xlim_act = [x_min, x_min+obj.xlim_diff];
                set(obj.ax, 'xlim', obj.xlim_act);
           end
           
        end
        
    end
end

