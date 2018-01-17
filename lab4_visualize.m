function lab4_visualize(data, params, moviefile)
%
%   data, params    the output of lab4_simulate
%
%   moviefile       a filename (e.g., 'movie.mp4') where a movie of the
%                   simulation will be saved - if this filename is empty, i.e.,
%                   if it is [], then no movie will be saved

% Set flag if we are making a movie
makemovie = ~isempty(moviefile);

% Load a description of the quadrotor (points and faces that make up
% a triangular mesh)
[pQuad_InBody, fQuad] = GetQuadModel('quadmodel.mat');

% Scale the quadrotor to match the spar length (model is scaled to l=0.5)
pQuad_InBody = pQuad_InBody*(params.l / 0.5);

% Load a description of the mocap system (one point for each camera)
pMocap_InRoom = GetMocapModel('mocapmodel.mat');

% Create a description of the room frame and the body frame (four points -
% at the origin and then at the end of the x, y, z unit vectors)
pRoomFrame_InRoom = 0.25*[zeros(3,1) eye(3)];
pBodyFrame_InBody = 0.25*[zeros(3,1) eye(3)];

% Create an empty figure
fig = [];

% Start making movie, if necessary.
if (makemovie)
    myV = VideoWriter(moviefile,'MPEG-4');
    myV.Quality = 100;
    myV.FrameRate = 30  ;
    open(myV);
    movie_nframes = 0;
    movie_tstep = 1/myV.FrameRate;
end

% Loop over all time
tic;
firsttime = 1;
while (1)
    
    % Get the current time
    if (makemovie)
        t_cur = movie_nframes * movie_tstep;
    else
        t_cur = toc;
    end
    
    % Find the last time step that has time less than the current time
    i = find(t_cur >= data.t, 1, 'last');
    
    % Update geometry
    [pQuad_InRoom, pBodyFrame_InRoom] = ...
        UpdateGeometry(data.t(i), data.o(:, i), data.theta(:, i), ...
                       pQuad_InBody, pBodyFrame_InBody);
    
    % Update figure
    fig = UpdateFigure(fig, data.t(i), data.o(:, i), ...
                            data.o_desired(:, i), data.o_goal, params.r, ...
                            data.obst{:, i}, ...
                            pMocap_InRoom, pQuad_InRoom, fQuad, ...
                            pRoomFrame_InRoom, pBodyFrame_InRoom);
	
    % Pause if not making a movie and if this is the first time step
	if (~makemovie && firsttime)
        firsttime = 0;
        pause;
        tic
	end
    
    % If making a movie, store the current figure as a frame
    if (makemovie)
        frame = getframe(gcf);
        writeVideo(myV,frame);
        movie_nframes = movie_nframes + 1;
    end
    
    % Check if done
    if (i == length(data.t))
        break;
    end
end

% Close and save the movie, if necessary.
if (makemovie)
    for i=1:myV.FrameRate
        frame = getframe(gcf);
        writeVideo(myV,frame);
    end
    close(myV);
end

end

function [pQuad_InRoom, pBodyFrame_InRoom] = ...
    UpdateGeometry(t, o, theta, pQuad_InBody, pBodyFrame_InBody)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE (COORDINATE TRANSFORMATION)
%
%  Inputs:
%
%   t                   time
%   o                   position of body frame in coordinates of room frame
%   theta               ZYX Euler Angles describing orientation of body
%                       frame in coordinates of room frame
%   pQuad_InBody        points in coordinates of body frame (3xM)
%   pBodyFrame_InBody   points in coordinates of body frame (3xN)
%
%  Outputs:
%
%   pQuad_InRoom        points in coordinates of room frame (3xM)
%   pBodyFrame_InRoom   points in coordinates of room frame (3xN)

%-----------Define Rotation Matrix------------
% yaw = theta(1)
% pitch = theta(2)
% roll = theta(3)
rz = [cos(theta(1)) -sin(theta(1)) 0; sin(theta(1)) cos(theta(1)) 0; 0 0 1];
ry = [cos(theta(2)) 0 sin(theta(2)); 0 1 0; -sin(theta(2)) 0 cos(theta(2))];
rx = [1 0 0; 0 cos(theta(3)) -sin(theta(3)); 0 sin(theta(3)) cos(theta(3))];
RQuadinRoom = rz*ry*rx;



% p_in1 = o0in1 + R0in1*p_in1
for c = 1:size(pQuad_InBody,2)
    pQuad_InRoom(:,c) = o + RQuadinRoom*pQuad_InBody(:,c);
end
for d = 1:size(pBodyFrame_InBody,2)
    pBodyFrame_InRoom(:,d) = o + RQuadinRoom*pBodyFrame_InBody(:,d);
end

end

function [fig] = UpdateFigure(fig, t, o, o_desired, o_goal, r, obst, ...
                              pMocap_InRoom, pQuad_InRoom, fQuad, ...
                              pRoomFrame_InRoom, pBodyFrame_InRoom)

    if isempty(fig)

        % Create figure
        clf;
        set(gcf,'renderer','opengl');
        set(gcf,'color','w');

        % Create an axis for text (it's important this is in the back,
        % so you can rotate the view and other stuff!)
        axes('position',[0 0 1 1]);
        axis([0 1 0 1]);
        hold on;
        axis off;
        fig.text=text(0.05,0.1,sprintf('t = %6.2f\n',t),'fontsize',10,...
                      'verticalalignment','top','fontname','monaco');

        % Create an axis for the view from the room frame
        axes();
        axis equal;
        axis([-4 4 -4 4 -3.5 0.1]);
        axis manual;
        hold on;

        % Reverse the y and z axes to get the "z down" view
        set(gca, 'ydir', 'reverse');
        set(gca, 'zdir', 'reverse');

        % Draw lights and labels
        lighting flat
        light('Position',[0 -2 -1])
        light('Position',[0 -2 1])
        xlabel('x');
        ylabel('y');
        zlabel('z');

        % Draw the room frame and the body frame
        fig.roomframe = DrawFrame([], pRoomFrame_InRoom);
        fig.bodyframe = DrawFrame([], pBodyFrame_InRoom);

        % Draw the floor
        fig.floor = patch([-2.5 2.5 2.5 -2.5 -2.5], ...
                          [-2.5 -2.5 2.5 2.5 -2.5], ...
                          0.9*[1 1 1]);

        % Draw the mocap system (put a point at the location of each camera)
        fig.mocap = scatter3(pMocap_InRoom(1,:), ...
                             pMocap_InRoom(2,:), ...
                             pMocap_InRoom(3,:), ...
                             15, 'k', 'filled');

        % Draw the quadrotor
        fig.quad = DrawMesh([], pQuad_InRoom, fQuad, 'y', 0.6);

        % Draw trace of position and of desired position
        fig.odes = DrawTrace([], nan(3,1), 'k', 200);
        fig.o = DrawTrace([], nan(3,1), [1, 0.6, 0], 200);

        % Create unit sphere
        [fig.xSphere,fig.ySphere,fig.zSphere]=sphere(16);
        [m,n]=size(fig.xSphere);

        % Create template array for sphere color
        c = ones(m,n,3);

        % Draw obstacles (all spheres)
        for i=1:length(obst)
            if (obst{i}.type == 1)
                c(:,:,1) = 0.75;
                c(:,:,2) = 0.25;
                c(:,:,3) = 0.25;
                fig.obst(i) = ...
                    DrawBubble([], obst{i}. p,obst{i}.s, ...
                                   fig.xSphere, fig.ySphere, fig.zSphere, ...
                                   c, 0.8);
            end
        end

        % Draw bubble at desired position
        c(:,:,1) = 0.25;
        c(:,:,2) = 0.25;
        c(:,:,3) = 0.75;
        fig.robotbubble = ...
            DrawGlowingBubble([], o_desired, r, ...
                                  fig.xSphere, fig.ySphere, fig.zSphere, ...
                                  c, 0.4);

        % Draw bubble at goal position
        c(:,:,1) = 0.25;
        c(:,:,2) = 0.75;
        c(:,:,3) = 0.25;
        fig.goalbubble = ...
            DrawBubble([], o_goal, r, ...
                           fig.xSphere, fig.ySphere, fig.zSphere, ...
                           c, 0.2);

        % Make everything look pretty
        h = findobj('Type','surface');
        set(h,'FaceLighting','gouraud',...
              'FaceColor','interp',...
              'EdgeColor',[.4 .4 .4],...
              'LineStyle','none',...
              'BackFaceLighting','lit',...
              'AmbientStrength',0.4,...
              'DiffuseStrength',0.6,...
              'SpecularStrength',0.5);
        material default

    else

        % Update figure
        % - update time
        set(fig.text,'string',sprintf('t = %6.2f', t));
        % - update body frame
        fig.bodyframe = DrawFrame(fig.bodyframe, pBodyFrame_InRoom);
        % - update quadrotor
        fig.quad = DrawMesh(fig.quad, pQuad_InRoom);
        % - update traces
        fig.odes = DrawTrace(fig.odes, o_desired);
        fig.o = DrawTrace(fig.o, o);
        % - update obstacles
        for i=1:length(obst)
            if (obst{i}.type == 1)
                fig.obst(i) = ...
                    DrawBubble(fig.obst(i), obst{i}.p, obst{i}.s, ...
                               fig.xSphere, fig.ySphere, fig.zSphere);
            end
        end
        % - update bubble at desired position
        fig.robotbubble = ...
            DrawGlowingBubble(fig.robotbubble, o_desired, r, ...
                              fig.xSphere, fig.ySphere, fig.zSphere);
    end

    % Tell MATLAB to update the figure window so we see what we just drew
    % on the screen immediately
    drawnow;

end

% Creates or updates a triangular mesh
function mesh = DrawMesh(mesh, p, f, color, alpha)
    if isempty(mesh)
        mesh = patch('Vertices',p','Faces',f,...
                     'FaceColor',color,'FaceAlpha',alpha,'EdgeAlpha',alpha);
    else
        set(mesh,'vertices',p');
    end
end

% Creates or updates three lines that describe a frame
function frame = DrawFrame(frame, p)
    if isempty(frame)
        frame.x = plot3(p(1,[1 2]),p(2,[1 2]),p(3,[1 2]),'r-','linewidth',3);
        frame.y = plot3(p(1,[1 3]),p(2,[1 3]),p(3,[1 3]),'g-','linewidth',3);
        frame.z = plot3(p(1,[1 4]),p(2,[1 4]),p(3,[1 4]),'b-','linewidth',3);
    else
        set(frame.x,'xdata',p(1,[1 2]),'ydata',p(2,[1 2]),'zdata',p(3,[1 2]));
        set(frame.y,'xdata',p(1,[1 3]),'ydata',p(2,[1 3]),'zdata',p(3,[1 3]));
        set(frame.z,'xdata',p(1,[1 4]),'ydata',p(2,[1 4]),'zdata',p(3,[1 4]));
    end
end

% Creates or updates position trace
function trace = DrawTrace(trace, p, c, n)
    if (isempty(trace))
        trace = line(p(1), p(2), p(3), ...
                     'color', c, 'marker', '.', 'markersize', 12, ...
                     'linestyle', 'none', 'UserData', n);
    else
        x = get(trace,'xdata');
        y = get(trace,'ydata');
        z = get(trace,'zdata');
        n = get(trace,'UserData');
        if (length(x)>=n)
            x = x(2:end);
            y = y(2:end);
            z = z(2:end);
        end
        x(:,end+1) = p(1);
        y(:,end+1) = p(2);
        z(:,end+1) = p(3);
        set(trace,'xdata',x,'ydata',y,'zdata',z);
    end
end

% Loads quadrotor geometry from a file
function [p, f] = GetQuadModel(filename)
    load(filename);
end

% Loads mocap geometry from a file
function p = GetMocapModel(filename)
    load(filename);
end


function bubble = DrawGlowingBubble(bubble,o,r,x,y,z,c,a)
    if (isempty(bubble))
        bubble.surface = surf(o(1)+r*x,o(2)+r*y,o(3)+r*z,c,'FaceAlpha',a);
        bubble.light = light('position',o,'style','local');
    else
        set(bubble.surface,'xdata',o(1)+r*x,'ydata',o(2)+r*y,'zdata',o(3)+r*z);
        set(bubble.light,'position',o);
    end
end

function bubble = DrawBubble(bubble,o,r,x,y,z,c,a)
    if (isempty(bubble))
        bubble = surf(o(1)+r*x,o(2)+r*y,o(3)+r*z,c,'FaceAlpha',a);
    else
        set(bubble,'xdata',o(1)+r*x,'ydata',o(2)+r*y,'zdata',o(3)+r*z);
    end
end

