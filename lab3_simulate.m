function [t, o, theta, odes,K] = lab3_simulate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   MODIFY
    %
    %   - Change parameter values and sample time to match your quadrotor (see
    %     your results from Labs #1 and #2)
    %   - Change initial time, final time, and initial conditions as you like
    %

    % Parameters
    g = 9.81;                 % acceleration of gravity m/s^2
    m = (19.365+5.845)*.0283;                  % mass kg
    J = diag([0.0023,0.0023,0.008]);    % moment of inertia matrix in body frame
    l = .384/2;                  % spar length, length of rotor arm (m)
    kF = 7.4e-6;              % aerodynamic force coefficient
    kM = 1.1e-7;              % aerodynamic torque coefficient
    sigmamax = 1e3;         % maximum spin rate

    % Initial time
    t0 = 0;

    % Initial state
    o0 = [0; 0; -1];
    theta0 = [0; 0; 0];
    v0 = [0; 0; 0];
    w0 = [0; 0; 0];

    % Sample time
    % from lab2 simulation
    dt = (1/50);

    % Final time
    t1 = 10;

    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   MODIFY
    %
    %   - Design a control policy (see HW2.5.2):
    %
    %       (xe, ue)    equilibrium point
    %       K           gain matrix
    %

    % equilibrium point
    % where we want the quadorotor to end up
    xe = [0;0;-2;
          0;0;0;
          0;0;0;
          0;0;0];
    ue = [0;0;0; m*g];
    
    %%%%% Compute K %%%%%
    syms o1 o2 o3 t1 t2 t3 w1 w2 w3 v1 v2 v3 u1 u2 u3 u4 
    
    pos   = [o1; o2; o3];
    theta = [t1; t2; t3];
    v     = [v1; v2; v3];
    w     = [w1; w2; w3];
    
    % symbolic state
    x = [pos;theta;v;w];
    % symbolic input
    u = [u1; u2; u3; u4]; 
    
    odot = v;
    %--Define R--
    rz = [cos(theta(1)) -sin(theta(1)) 0; sin(theta(1)) cos(theta(1)) 0; 0 0 1];
    ry = [cos(theta(2)) 0 sin(theta(2)); 0 1 0; -sin(theta(2)) 0 cos(theta(2))];
    rx = [1 0 0; 0 cos(theta(3)) -sin(theta(3)); 0 sin(theta(3)) cos(theta(3))];
    R = rz*ry*rx;
    %--Defining N--
    % first axis of rotation:
    R_3in1 = ry*rx;
    R_1in3 = transpose(R_3in1);
    % unit vector
    e1 = [0;0;1];
    column1 = R_1in3*e1;
    % second axis of rotation
    R_3in2 = rx;
    R_2in3 = transpose(R_3in2);
    % unit vector
    e2 = [0;1;0];
    column2 = R_2in3*e2;
    % third axis
    % unit vector
    e3 = [1;0;0];
    column3 = e3;
    N = ([column1 column2 column3])^-1;
    % -------------
    % angular velocity
    thetadot = N*w;
    % net force (positive z is downward)
    f = [0; 0; m*g]+R*[0; 0; -u(4)];
    % linear acceleration
    vdot = 1/m*f;
    % angular acceleration
    tau = [u(1); u(2); u(3)];
    w_wedge = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    wdot = inv(J)*(tau-w_wedge*J*w);
    
    xdot = pack_state(odot, thetadot, vdot, wdot);
    hsym = xdot;
    
    % compute Ac and Bc
    Ac = double((subs(jacobian(hsym,x),[x;u],[xe;ue])));
    Bc = double((subs(jacobian(hsym,u),[x;u],[xe;ue])));

    % dicretize (compute Ad and Bd)
    [Ad, Bd] = ssdata(c2d(ss(Ac,Bc,[],[]),dt));

    Q = eye(12);
    R = eye(4);
    % Compute the optimal gain matrix, (a row vector)
    K = dlqr(Ad, Bd, Q, R);
    %%%%%%%%%%%%%%%%%%%%%%

    
    % Create variables to keep track of time, state, input, and desired position
    t = [t0];
    x = [o0; theta0; v0; w0];
    u = [];
    odes = [];
    % Iterate over t1/dt sample intervals.
    t1 = 10;
    for i = 1:(t1/dt)

        % Get time and state at start of i'th sample interval
        ti = t(:, i);
        xi = x(:, i);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %   MODIFY
        %
        %   - Get desired position at start of i'th sample interval
        %
        %     (You may also need to redefine your equilibrium point!)
        %
        
        % desired position (can be defined as a function of time (ti))
        % desired position moves in a circle at a constant height
        odesi = [cos(ti); sin(ti); -2];     
        xe(1:3) = odesi;
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        odes(:, i) = odesi;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %   MODIFY
        %
        %   - Get input that will be applied throughout the i'th sample
        %     interval (in other words, implement your control policy)
        %   - Don't forget to make sure that this input can be realized by
        %     spin rates that are between 0 and sigmamax
        %
        % calculate desired input
        u_desired = ue - K * (xi-xe);
        % check if u_desired is achievable and return best option
        ui = GetBoundedInputs(u_desired, kF, kM, l, sigmamax);

        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u(:, i) = ui;

        % Get time and state at start of (i+1)'th sample interval
        [tsol, xsol] = ode45(@(t, x) h(t, x, ui, g, m, J), [ti ti+dt], xi);
        t(:, i+1) = tsol(end, :)';
        x(:, i+1) = xsol(end, :)';
    end

    % Get position and orientation
    o = x(1:3, :);
    theta = x(4:6, :);

end

function xdot = h(t, x, u, g, m, J)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   MODIFY
    %
    %   - Compute xdot given x, u, and other parameters (see HW2.2)

    [o, theta, v, w] = unpack_state(x);
    odot = v;
    %--Define R--
    rz = [cos(theta(1)) -sin(theta(1)) 0; sin(theta(1)) cos(theta(1)) 0; 0 0 1];
    ry = [cos(theta(2)) 0 sin(theta(2)); 0 1 0; -sin(theta(2)) 0 cos(theta(2))];
    rx = [1 0 0; 0 cos(theta(3)) -sin(theta(3)); 0 sin(theta(3)) cos(theta(3))];
    R = rz*ry*rx;
    %---Defining N---
    R_3in1 = ry*rx;
    R_1in3 = transpose(R_3in1);
    e1 = [0;0;1];
    R_3in2 = rx;
    R_2in3 = transpose(R_3in2);
    e2 = [0;1;0];
    e3 = [1;0;0];
    N = ([R_1in3*e1, R_2in3*e2, e3]);
    % --------------
    % angular velocity
    thetadot = inv(N)*w;
    % net force (positive z is downward)
    f = [0; 0; m*g]+R*[0; 0; -u(4)];
    % linear acceleration
    vdot = 1/m*f;
    % angular acceleration
    tau = [u(1); u(2); u(3)];
    w_wedge = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    wdot = inv(J)*(tau-w_wedge*J*w);
    
    xdot = pack_state(odot, thetadot, vdot, wdot);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [u, sigma] = GetBoundedInputs(u_desired, kF, kM, l, sigmamax)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   MODIFY
    %
    %   - Compute an input u (4x1) that is close to the input u_desired but
    %     that can be realized by spin rates that are between 0 and sigmamax
    %   - Compute the spin rates sigma (4x1) - these are not used by the rest
    %     of the simulation code, but will be necessary in your onboard C code,
    %     so it is useful to make sure you know how to compute them here
    %
    %     (See HW2.1.2)
    W = [l*kF -l*kF 0 0;
         0 0 l*kF -l*kF;
         kM kM -kM -kM;
         kF kF kF kF];
    % calculate desired spin rate
    s_desired = inv(W)*u_desired;
    for i = 1:length(s_desired)
        if s_desired(i) < 0
            s_desired(i) = 0;
        elseif s_desired(i) > sigmamax^2
            s_desired = sigmamax^2;
        end
    end
    u = W*s_desired;
    sigma = s_desired.^.5;
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
function x =  pack_state(o,theta,v,w)
    x = [o; theta; v; w];
end
function [o,theta,v,w] = unpack_state(x)
    o = x(1:3,:);
    theta = x(4:6,:);
    v = x(7:9,:);
    w = x(10:12,:);
end
