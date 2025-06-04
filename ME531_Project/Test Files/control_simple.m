function [next_state, next_input] = control_simple(init_positions, initial_vel, ref, t_step, t_max, verbose)

    clc, close all

    %example commands to test 
    %in = [1 2 3 4 5 6 7 8 9 10 11 12]
    %r = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2]
    %control_simple(in, in_v, r, 0.001, 10, true)

    %in = [4 4 4 4 4 4 4 4 4 4 4 4]
    %r = [4.01 4.01 4.01 4.01 4.01 4.01 4.01 4.01 4.01 4.01 4.01 4.01]
    %in_v = [0 0 0 0 0 0 0 0 0 0 0 0]
    
    %known values:
    m = 1; b =1; k =1; 

    if verbose == true
        %tspan = 0:t_step:t_max; %calculates a pretty graph over a selected 
        tspan = [0 t_max];
        % amount of seconds, but actually only sends back the next steps
    else
        tspan = [0 t_step t_step*2];
    end
    Q = [1 0
         0 1]; %Q was chosen to be very aggressive
    R = 1;
    C = [1 0];
    B = [0; 1/m];

    init = [init_positions; initial_vel];
    ref = [ref; zeros(1,length(init))];

    A = [0, 1; -k/m, -b/m];
        if verbose == true
            figure
            hold on
        end
        %for each entry in the initial values we run an LQR controller
        for n = 1:length(init)
            s = ss(A,B,C,0);
            [K_lqr,~,~] = lqr(A, B, Q,R);
            u_lqr = @(x) - K_lqr*(x-ref(:,n)); % control law
            [t,x] = ode45(@(t,x)springmass(x,m, b, k, u_lqr(x)),tspan,init(:,n));
 
            if verbose == true
                plot(t,x(:,1))
            end
            
            for q = 1:length(x(:,1))
                input_line = u_lqr(x(q,:));
                input(q,1) = input_line(1);
                input(q,2) = input_line(2);
            end
            all_inputs{n} = input(:,2);
            all_t{n} = t;
            all_positions{n} = x;
        end
        if verbose == true
            title('Position LQR')
            xlabel('Time (sec)')
            ylabel('Position (m)')
            figure
            hold on
            for n = 1:length(init)
                plot(all_t{n},all_positions{n}(:,2))
            end
            title('Velocity LQR')
            xlabel('Time (sec)')
            ylabel('Velocity (m/s)')
            figure 
            hold on
            input = [];
            for j = 1:length(init)
                plot(all_t{j}, all_inputs{j})
            end
            title("Control Input")
            ylabel("Control Input")
            xlabel("Time (sec)")
            % figure
            % hold on
        end


    next_state = [];
    next_input = [];
    for l =1:length(all_positions)
        next_state(1,l) = all_positions{l}(2,1);%gets the 2nd entry in the states, so the one right after the initial one, or the next step/state
        next_state(2,l) = all_positions{l}(2,2);
        next_input = [next_input, all_inputs{l}(1,:)]; %gets the first input applied, which should be the one that is connected to the 2nd state/next state
        %should we actually be getting the second entry? The first are all
        %the same?
    end

    function dx = springmass(x,m, b,k,u) %actual function
        dx(1,1) = x(2);
        dx(2,1) = u-b/m*x(2) - k/m*x(1);
    end
end