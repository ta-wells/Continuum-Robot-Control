function [all_t, all_positions] = control(init, ref, random)
    %init and ref need to be a cell of 2x1 values ie, in = {[1;1] [2;2]}
    %random is a true or false value that allows the noise aspect to be
    %turned off or on

    clc, close all
    %known values:
    m = 1; b =1; k =1; 

    %maximum 5% change for k and b, can be tuned as neeeded
    k_percent = 0.05;
    b_percent =0.05;
    
    %finds zeta and omega for pole placing later on
    zeta = b/(2*sqrt(k*m));
    omega_n = k/m;
    
    %sets some of the variables and aspects for the LQR and matrices
    tspan = 0:.1:10;
    Q = [10 0 ; 0 100];
    R = 2000;
    C = [1 0];
    B = [0; 1/m];

    %if random is true, then we don't just graph the values of b and k but
    %also we randomly sample between min and max k and be (set by the
    %percentage) and then also graph -3, -2, -1 sigma and +3, +2, +1 sigma
    if random == true
        %uses the function "actuator noise" to get a list of values that
        %have a range of different sigma values
        ks = actuator_noise(k, k_percent);
        bs = actuator_noise(b, b_percent);
        figure
        hold on
        for a =1 :length(ks)
            A = [0, 1; -ks(a)/m, -bs(a)/m]; %A matrix is changing
            
            for k = 1:length(init)
                s = ss(A,B,C,0);
                [K_lqr,S,e] = lqr(s,Q,R); %uses LQR controller
                u_lqr = @(x) - K_lqr*(x-ref{k}); % control law
                [all_t{k},all_positions{k}] = ode45(@(t,x)springmass(x,m, b, k, u_lqr(x)),tspan,init{k});
                plot(all_t{k},all_positions{k}(:,1))
            end
       
        end
        title('Position LQR')
        xlabel('Time (sec)')
        ylabel('Position (m)')
        figure
        hold on
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
                u_pole = @(x)-K_pole*(x-ref{m}); % control l [t, positions] = control(in, r)aw
                [t,x] = ode45(@(t,x)springmass(x,m, b, k, u_pole(x)),tspan,init{m});
                plot(t,x(:,1))
            end
        end
        title('Position Poles')
        xlabel('Time (sec)')
        ylabel('Position (m)')
    
    %if random is false then we don't bother taking a sample, we just graph
    %with the values of k and b without any noise
    else
        %The A matrix only needs to be defined once
        A = [0, 1; -k/m, -b/m];

        figure
        hold on
        %for each entry in the intial values we run an LQR controller
        for k = 1:length(init)
            s = ss(A,B,C,0);
            [K_lqr,S,e] = lqr(s,Q,R);
            u_lqr = @(x) - K_lqr*(x-ref{k}); % control law
            [all_t{k},all_positions{k}] = ode45(@(t,x)springmass(x,m, b, k, u_lqr(x)),tspan,init{k});
            plot(all_t{k},all_positions{k}(:,1))
            % u = @(x)-K*(x-r); % control law
            % [t,x] = ode45(@(t,x)springmass(x,m, b, k, u(x)),tspan,x0);
        end
        title('Position LQR')
        xlabel('Time (sec)')
        ylabel('Position (m)')
        
        %we also can compare the solution that does pole placement rather
        %than an LQR controller
        figure
        hold on
        for m = 1:length(init)
            %for each entry in the list of initial coord we find then place
            %our poles
            p1 = -zeta*omega_n + 1i*omega_n*sqrt(1-zeta^2);
            p2 = -zeta*omega_n - 1i*omega_n*sqrt(1-zeta^2);
            K_pole = place(A,B,[p1 p2]); %finds the gains
            u_pole = @(x)-K_pole*(x-ref{m});
            [t,x] = ode45(@(t,x)springmass(x,m, b, k, u_pole(x)),tspan,init{m});
            plot(t,x(:,1))
        end
        title('Position Poles')
        xlabel('Time (sec)')
        ylabel('Position (m)')
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