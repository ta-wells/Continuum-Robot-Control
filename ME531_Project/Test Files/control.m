function [next_state, next_input] = control(init_positions, initial_vel, ref, t_step, t_max, random, verbose)
    %init and ref are 1xn arrays of position information
    %random is a true or false value that allows the noise aspect to be
    %turned off or on
    %verbose plots the graphs

    %instead of a tspan, I have the code request a step size.
    %the code will plot a tspan of ten seconds But the output is the very
    %next out put based on tstep. So the plots are really just for show if
    %you gave the controller more time

    clc, close all

    %example commands to test 
    %in = [1 2 3 4 5 6 7 8 9 10 11 12]
    %r = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2]
    %control(in, in_v, r, 0.001, 10, false, true)

    %in = [4 4 4 4 4 4 4 4 4 4 4 4]
    %r = [4.01 4.01 4.01 4.01 4.01 4.01 4.01 4.01 4.01 4.01 4.01 4.01]
    %in_v = [0 0 0 0 0 0 0 0 0 0 0 0]
    
    %known values:
    m = 1; b =1; k =1; 

    %maximum 5% change for k and b, can be tuned as neeeded
    k_percent = 0.05;
    b_percent =0.05;
    
    %finds zeta and omega for pole placing later on
    zeta = b/(2*sqrt(k*m));
    omega_n = k/m;
    
    %sets some of the variables and aspects for the LQR and matrices
    
    if verbose == true
        tspan = 0:t_step:t_max; %calculates a pretty graph over a selected 
        %tspan = [0 t_max];
        % amount of seconds, but actually only sends back the next steps
    else
        tspan = [0 t_step t_step*2];
    end
    %tspan = [0 t_step t_step*2];%if you don't want the graphs just
    %uncomment this and it will make only find a tspan of a couple of
    %entries
    %the lqr controller does weird stuff if it's only 2 entries

    Q = [5000 0 ; 0 5000]; %Q was chosen to be very aggressive
    R = 0.01;
    %Q = [100 0 ; 0 100]; %Q was chosen to be very aggressive
    %R = 0.001;
    C = [1 0];
    B = [0; 1/m];

    %the function is past an array 1x12 array for both reference and
    %initial positions, we need to make this a 2x12 array to include the
    %velocities which will be 0 for the reference array but up to user
    %input for the intial values
    init = [init_positions; initial_vel];
    ref = [ref; zeros(1,length(init))];
    
    %if random is true, then we don't just graph the values of b and k but
    %also we randomly sample between min and max k and be (set by the
    %percentage) and then also graph -3, -2, -1 sigma and +3, +2, +1 sigma
    if random == true
        %uses the function "actuator noise" to get a list of values that
        %have a range of different sigma values
        ks = actuator_noise(k, k_percent);
        bs = actuator_noise(b, b_percent);
        if verbose ==true
            figure
            hold on
        end
        for a =1 :length(ks)
            A = [0, 1; -ks(a)/m, -bs(a)/m]; %A matrix is changing
            
            for c = 1:length(init)
                s = ss(A,B,C,0);
                [K_lqr,S,e] = lqr(s,Q,R); %uses LQR controller
                u_lqr = @(x) - K_lqr*(x-ref(:,c)); % control law
                
                [all_t{c},all_positions{c}] = ode45(@(t,x)springmass(x,m, b, k, u_lqr(x)),tspan,init(:,c));
                if a == 4 %save the array for when it is just the mean
                    mean_all_t = all_t;
                    mean_all_positions = all_positions;
                end
                if verbose == true
                    plot(all_t{c},all_positions{c}(:,1))
                end
            end
       
        end
        if verbose == true
            title('Position LQR')
            xlabel('Time (sec)')
            ylabel('Position (m)')
            figure
            hold on
        end
        input = [];
        for j = 1:length(init) %use the mean value solution to 
                % calculate inputs, we could technically do it for all of 
                % the sifferent z-scores but I don't really see why
            for q = 1:length(mean_all_positions{j}(:,1))
                input_line = u_lqr(mean_all_positions{j}(q,:));
                input(q,1) = input_line(1);
                input(q,2) = input_line(2);
            end
            mean_all_inputs{j} = input(:,2);
            if verbose == true
                plot(mean_all_t{j}, input(:,2));
            end
        end
        all_positions = mean_all_positions; %since there are multiple 
            % arrays with the name all positions we have to clarify that we
            % want to save the mean or 0 zscore one
        all_inputs = mean_all_inputs;
        if verbose == true
            title("Control Input")
            ylabel("Control Input")
            xlabel("Time (sec)")
            figure
            hold on
        end
        %we also can compare the solution that does pole placement rather
        %than an LQR controller
        for c =1:length(ks) 
            A = [0, 1; -ks(c)/m, -bs(c)/m];
            for m = 1:length(init)
                %for each entry in the list of initial coord we find then place
            %our poles
                p1 = -zeta*omega_n + 1i*omega_n*sqrt(1-zeta^2);
                p2 = -zeta*omega_n - 1i*omega_n*sqrt(1-zeta^2);
                K_pole = place(A,B,[p1 p2]); %finds the gains
                u_pole = @(x)-K_pole*(x-ref(:,c)); % control l [t, positions] = control(in, r)aw
                [t,x] = ode45(@(t,x)springmass(x,m, b, k, u_pole(x)),tspan,init(:,c));
                if verbose == true
                    plot(t,x(:,1))
                end
            end
        end

        if verbose == true
            title('Position Poles')
            xlabel('Time (sec)')
            ylabel('Position (m)')
        end
    %if random is false then we don't bother taking a sample, we just graph
    %with the values of k and b without any noise
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %The A matrix only needs to be defined once
        A = [0, 1; -k/m, -b/m];
        if verbose == true
            figure
            hold on
        end
        %for each entry in the initial values we run an LQR controller
        for n = 1:length(init)
            clear u_lqr
            clear x
            s = ss(A,B,C,0);
            [K_lqr,~,~] = lqr(s,Q,R);
            u_lqr = @(x) - K_lqr*(x-ref(:,n)); % control law
            [all_t{n},all_positions{n}] = ode45(@(t,x)springmass(x,m, b, k, u_lqr(x)),tspan,init(:,n));
   
            if verbose == true
                plot(all_t{n},all_positions{n}(:,1))
            end
            for q = 1:length(all_positions{n}(:,1))
                input_line = u_lqr(all_positions{n}(q,:));
                input(q,1) = input_line(1);
                input(q,2) = input_line(2);
            end
            all_inputs{n} = input(:,2);

        end
        if verbose == true
            title('Position LQR')
            xlabel('Time (sec)')
            ylabel('Position (m)')
            figure
            hold on
            input = [];
            for j = 1:length(init)
                plot(all_t{j}, all_inputs{j})
            end
            title("Control Input")
            ylabel("Control Input")
            xlabel("Time (sec)")
            figure
            hold on
        end
        %we also can compare the solution that does pole placement rather
        %than an LQR controller
        % 
        for v = 1:length(init)
            %for each entry in the list of initial coord we find then place
            %our poles
            p1 = -zeta*omega_n + 1i*omega_n*sqrt(1-zeta^2);
            p2 = -zeta*omega_n - 1i*omega_n*sqrt(1-zeta^2);
            K_pole = place(A,B,[p1 p2]); %finds the gains
            u_pole = @(x)-K_pole*(x-ref(:,v));
            [t,x] = ode45(@(t,x)springmass(x,m, b, k, u_pole(x)),tspan,init(:,v));
            if verbose == true
                plot(t,x(:,1))
            end
        end
        if verbose == true
            title('Position Poles')
            xlabel('Time (sec)')
            ylabel('Position (m)')
        end
    end

    %writes out the values
    next_state = [];
    next_input = [];
    for l =1:length(all_positions)
        next_state(1,l) = all_positions{l}(2,1);%gets the 2nd entry in the states, so the one right after the initial one, or the next step/state
        next_state(2,l) = all_positions{l}(2,2);
        next_input = [next_input, all_inputs{l}(1,:)]; %gets the first input applied, which should be the one that is connected to the 2nd state/next state
        %should we actually be getting the second entry? The first are all
        %the same?
    end

    %this function aims to simulate the possible range of vlaues we could
    %get from sensor noise.
    function [ks, mu] = actuator_noise(k, percent)
        %creates an upper and lower bound using the percentage change
        lower_bound_k = k*(1-percent);
        upper_bound_k = k*(1+percent);
        %r is data set from random samples bound to the range we set
        r = (upper_bound_k-lower_bound_k).*rand(5000,1) + lower_bound_k;
        r_range = [min(r) max(r)]; %checks the range, technically superflous but good to double check everything
        sigma = std(r); %finds sigma and the mean from the data set
        mu = mean(r);
        %since data set was random we can assume gaussian distribution
        sigmas = [-3, -2, -1, 0, 1, 2, 3]; %sigma values/z-scores of interest
        %based on the sigma values of interest, we find the actual values
        %of k or b that match the sigma values/z-scores
        for i = 1:length(sigmas)
            ks(i) = sigmas(i)*sigma+k;
        end
    end

    function dx = springmass(x,m, b,k,u) %actual function
        dx(1,1) = x(2);

        dx(2,1) = u-b/m*x(2) - k/m*x(1);
    end
end