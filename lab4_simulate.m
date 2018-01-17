function [data, params] = lab4_simulate
    
    % Parameters
    % gravity
    params.g = 9.81;
    % - mass
    params.m = (19.365+5.845)*.0283;
    % - moment of inertia
    params.J = diag([0.0023,0.0023,0.008]);
    % - spar length (m)
    params.l = .384/2;
    % - aerodynamic force and moment coefficients
    params.kF = 6.616040610072618e-06;
    params.kM = 1.1e-7; 
    % - maximum spin rate
    params.sigmamax = (1e3);
    % - radius of bounding volume (sphere) around quadrotor
    params.r = 3*params.l;
    % - sample time
    params.dt = (1/50);
    
    % Simulation
    % - initial time
    t0 = 0;
    % - initial state
    o0 = [-1.5; 1; -2];
    theta0 = [0; 0; 0];
    v0 = [0; 0; 0];
    w0 = [0; 0; 0];
    x0 = [o0; theta0; v0; w0];
    % - final time
    tf = 10;
    
    % Problem
    % - desired position
    o_desired = [-1.5; 1.5; -2];
    % - goal position
    o_goal = [1.5; -1.0; -1.5];
    % - obstacles
    %   * create an empty cell array to hold obstacles
    obst = {};
    %   * uncomment this line to add a spherical obstacle (center, radius)
    % obst = AddObstacle_Sphere(obst, [0; 1; -2], 0.5);
    %   * add n random spherical obstacles
    n = 1;
    obst = AddObstacle_RandomSpheres(obst, n, 0.1, 0.5, 2.5, 2.5, 2.5, ...
                                     o_desired, o_goal, params);
	
    % DESIGN
    %
    %     %%%%%%%%%
    %     % FIXME %
    %     %%%%%%%%%
    %
    % - controller
    
    params.xe = [0;0;0;
                 0;0;0;
                 0;0;0;
                 0;0;0];
    params.ue = [0;0;0; params.m*params.g];
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
        f = [0; 0; params.m*params.g]+R*[0; 0; -u(4)];
        % linear acceleration
        vdot = 1/params.m*f;
        % angular acceleration
        tau = [u(1); u(2); u(3)];
        w_wedge = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
        wdot = inv(params.J)*(tau-w_wedge*params.J*w);
        xdot = pack_state(odot, thetadot, vdot, wdot);
        hsym = xdot;
        % compute Ac and Bc
        Ac = double((subs(jacobian(hsym,x),[x;u],[params.xe;params.ue])));
        Bc = double((subs(jacobian(hsym,u),[x;u],[params.xe;params.ue])));
        % dicretize (compute Ad and Bd)
        [Ad, Bd] = ssdata(c2d(ss(Ac,Bc,[],[]),params.dt));
        Qd = diag([20,20,20,10,10,10,1,1,1,1,1,1]);
        Rd = eye(4);
    params.K = dlqr(Ad,Bd,Qd,Rd);
    
  %  - planner
    params.k_att = 20; % speed of sphere
    params.b_att = 0.1;
    params.k_rep = 20;   % pushing away from obstacles
    params.b_rep = 1;
    params.k_des = 0.02;  % step
    params.b_des = 0.02;
    
    % Data to store
    data.t = [t0];
    data.x = [x0];
    data.u = [];
    data.o_desired = [o_desired];
    data.o_goal = [o_goal];
    data.obst = {obst};
    
    % Iterate over t1/dt sample intervals.
    for i = 1:round(tf/params.dt)

        % Current time and state
        t = data.t(:, end);
        x = data.x(:, end);
        %%%%% Move obstacle as a function of time %%%%%
        obst{1}.p = [cos(t); sin(t); -0.5];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%% Follow the 50cm above the object %%%%%
        o_goal = obst{1}.p + [0;0;-0.5];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Planner (update o_desired)
        o_desired = planner(t, x, o_desired, o_goal, obst, params);
        
        % Controller (update u)
        u = controller(t, x, o_desired, params);
        
        % Simulator (update t and x)
        [tsol, xsol] = ode45(@(t, x) h(t, x, u, params.g, params.m, params.J), [t t+params.dt], x);
        
        % Store data
        data.o_desired(:, end+1) = o_desired;
        data.u(:, end+1) = u;
        data.t(:, end+1) = tsol(end, :)';
        data.x(:, end+1) = xsol(end, :)';
        data.obst{:, end+1} = obst;

    end

    % Get position and orientation
    data.o = data.x(1:3, :);
    data.theta = data.x(4:6, :);

end

function obst = AddObstacle_Sphere(obst, p, s)
    obst{end+1} = struct('type', 1, 'p', p, 's', s);
end

function obst = ...
    AddObstacle_RandomSpheres(obst, ...
                              n, smin, smax, dx, dy, dz, ...
                              o_start, o_goal, params)

	% n = number of spheres to add
    % smin = minimum radius
    % smax = maximum radius
    % o_start = start position
    % o_goal = goal position
    % (-dx, dx), (-dy, dy), (-dz, 0) = room dimensions

    % Loop to add n spheres one at a time
    for i=1:n

        % Keep trying until we add a sphere that does not interfere
        % with the start or goal position of the quadrotor.
        while(1)
            % Random center position in the room
            p = [dx;dy;dz].*([-1;-1;-1]+[2;2;1].*rand(3,1));
            % Random radius
            s = 0.2;
            % Check for non-interference
            ds = norm(p-o_start)-(params.r+s);
            dg = norm(p-o_goal)-(params.r+s);
            if ((ds>0)&&(dg>0.5))
                obst = AddObstacle_Sphere(obst, p, s);
                break;
            end
        end

    end
end

function xdot = h(t, x, u, g, m, J)
    
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
       
end

function u = controller(t, x, o_desired, params)
    xe = [o_desired;zeros(9,1)];
    u_desired = params.ue-params.K*(x-xe);
    [u,sigma] = GetBoundedInputs(u_desired,params);
end

function o_desired = planner(t, x, o_desired, o_goal, obst, params)
    
    q = o_desired;
    q_goal = o_goal;
    r = params.r;
    % compute attractive gradient
    if norm(q-q_goal)<= params.b_att
        gradf = params.k_att*(q-q_goal);
    else
        gradf = params.k_att*params.b_att*(q-q_goal)/norm(q-q_goal);
    end
 
    for i=1:length(obst)
        % get center (p) and radius (s) of i'th spherical obstacle
        p = obst{i}.p;
        s = obst{i}.s;
        d = norm(q-p)-(r+s);
        dgrad = (q-p)/norm(q-p);
        % compute repulsive gradient
        if d <= params.b_rep
           gradf = gradf - params.k_rep*(1/d-1/params.b_rep)*((1/d)^2)*dgrad;
        end
    end

    if norm(params.k_des*gradf) <= params.b_des
        q = q - params.k_des*gradf;
    else
        q = q - params.b_des*gradf/norm(gradf);
    end
    o_desired = q;

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
function [u, sigma] = GetBoundedInputs(u_desired, params)

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
    l = params.l;
    kF = params.kF;
    kM = params.kM;
    W = [l*kF -l*kF 0 0;
         0 0 l*kF -l*kF;
         kM kM -kM -kM;
         kF kF kF kF];
    % calculate desired spin rate
    s_desired = inv(W)*u_desired;
    for i = 1:length(s_desired)
        if s_desired(i) < 0
            s_desired(i) = 0;
        elseif s_desired(i) > params.sigmamax^2
            s_desired = params.sigmamax^2;
        end
    end
    u = W*s_desired;
    sigma = s_desired.^.5;
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



