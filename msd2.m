%% Deformable object with interconnected mass-spring-damper
%
%  Author : Auralius Manurung (manurung.auralius@gmail.com)
%  Note   : Allow user interaction using Euler Integration
% 

%% Data structures for the nodes
%  nodes.r
%  nodes.c
%  nodes.node.intialPos
%  nodes.node.pos
%  nodes.node.force
%  nodes.node.vel
%  nodes.node.acc
%  nodes.node.force_ext
%
% 8 interconnections:
%  O   O   O
%    \ | /
%  O---O---O
%    / | \
%  O   O   O


%%
function msd2()
% This is the main function

    % Warning! This clears everything!
	close all;
    clc;
    clear all;
    
    % For plotting purpose, we will aslo display the current simulation
    % time
    S.f = figure;
    S.a = axes;
    S.h = plot(0, 0);
    S.mText = uicontrol('style','text');
    
    % Data structures that are embedded to S, to allow data exchange among
    % callbacks
    S.pt1 = [];      % Mouse drag start position
    S.pt2 = [];      % Mouse drag current position 
    S.k = 0;         % 1 when key is down, 0 when key is released
    S.trigger = 0;   % 1 when key is down    
    guidata(S.f, S);
    
    xlabel('meter');
    ylabel('meter');
        
    % The calbacks
    set(S.h,'ButtonDownFcn',@(varargin)startDragFcn(varargin, S))
    set(S.f, 'WindowButtonUpFcn', @(varargin)stopDragFcn(varargin, S));
        

    % Parameters
    row = 5;
    col = 5;
    stiffness = 10;           % N/m
    damping =   1;            % Ns/m
    mass = 0.01;              % Kg
    ts = 0.005;               % Seconds
    virtualSpringConst = 10;
    
    % Build the nodes and the canvas
    nodes = buildNodes(row, col);
    canvas = createCanvas(nodes);
    canvas = drawNodes(S, canvas, nodes, 0);
    
    % Selected node by mouse
    rs = 1;
    cs = 1;
           
    % This is the main iterations
    for i = 0 : 10000
        % Do update
        nodes = updateNode(nodes, mass, stiffness, damping, ts);
        
        if mod(i, 10) == 0
            canvas = drawNodes(S, canvas, nodes, ts*i);
            
            S = guidata(S.f);  
                        
            if S.trigger == 1
                S.trigger = 0;
                [rs, cs] = findClosestNodeFromPt(nodes, S.pt1);
            end
            
            if nodes.node(rs, cs).isFixed == 0                
                nodes.node(rs, cs).force_ext =  (S.pt2 -nodes.node(rs, cs).pos) .* (S.k * virtualSpringConst);
            end
                        
            guidata(S.f, S);
        end
    end
        
end

%%
function nodes = buildNodes(row, col)
% Build the nodes 
    nodes.row = row;
    nodes.col = col;
    
    for c = 1: col
        for r = 1 : row
            node(r,c).initalPos = [(c - 1) / 100 (r - 1) / 100 ]; % 1 cm step
            node(r,c).pos = node(r,c).initalPos;
            node(r,c).acc = [0 0];
            node(r,c).vel = [0 0];
            node(r,c).force_ext = [0 0];

            % The last row is fixed
            if (r == 1)
                node(r,c).isFixed = 1;
            else
                node(r,c).isFixed = 0;
            end
        end
    end
    
    nodes.node = node;
end

%% 
function nodes = updateNode(nodes, mass, stiffness, damping, ts)
% Update all nodes per time sampling
    row = nodes.row;
    col = nodes.col;
    node = nodes.node;
    
    % Force update
    % Calculate force on each node
    for r = 1 : row
        nextRow = r + 1;
        prevRow = r - 1;
        
        for c = 1 : col
            nextCol = c + 1;
            prevCol = c - 1;
            
            f1 = [0 0];
            f2 = [0 0];
            f3 = [0 0];
            f4 = [0 0];
            f5 = [0 0];
            f6 = [0 0];
            f7 = [0 0];
            f8 = [0 0];

            % Link 1
            if (r < row && c > 1)
                l0 = node(r, c).initalPos - node(nextRow, prevCol).initalPos;
                lt = node(r, c).pos - node(nextRow, prevCol).pos;
                n = norm(lt, 2);                
                f1 = stiffness * (norm(l0, 2) - n) * lt / n;
            end

            % Link 2
            if (r < row)
                l0 = node(r, c).initalPos - node(nextRow, c).initalPos;
                lt = node(r, c).pos - node(nextRow, c).pos;
                n = norm(lt, 2);
                f2 = stiffness * (norm(l0, 2) - n) * lt / n;
            end

            % Link 3
            if (c < col)
                l0 = node(r, c).initalPos - node(r, nextCol).initalPos;
                lt = node(r, c).pos - node(r, nextCol).pos;
                n = norm(lt, 2);
                f3 = stiffness * (norm(l0, 2) - n) * lt / n;
            end

            % Link 4
            if (r > 1 && c < col)
                l0 = node(r, c).initalPos - node(prevRow, nextCol).initalPos;
                lt = node(r, c).pos - node(prevRow, nextCol).pos;
                n = norm(lt, 2);
                f4 = stiffness * (norm(l0, 2) - n) * lt / n;
            end

            % Link 5
            if (r > 1)
                l0 = node(r, c).initalPos - node(prevRow, c).initalPos;
                lt = node(r, c).pos - node(prevRow, c).pos;
                n = norm(lt, 2);
                f5 = stiffness * (norm(l0, 2) - n) * lt / n;    
            end

            % Link 6
            if (c > 1)
                l0 = node(r, c).initalPos - node(r, prevCol).initalPos;                        
                lt = node(r, c).pos - node(r, prevCol).pos; 
                n = norm(lt, 2);
                f6 = stiffness * (norm(l0, 2) - n) * lt / n;                     
            end
            
            % Link 7
            if (r < row && c < col)
                l0 = node(r, c).initalPos - node(nextRow, nextCol).initalPos;                        
                lt = node(r, c).pos - node(nextRow, nextCol).pos; 
                n = norm(lt, 2);
                f7 = stiffness * (norm(l0, 2) - n) * lt / n;                     
            end
            
            % Link 8
            if (r > 1 && c > 1)
                l0 = node(r, c).initalPos - node(prevRow, prevCol).initalPos;                        
                lt = node(r, c).pos - node(prevRow, prevCol).pos; 
                n = norm(lt, 2);
                f8 = stiffness * (norm(l0, 2) - n) * lt / n;                     
            end

            node(r,c).force =  f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 - ... 
                               damping * node(r,c).vel + mass * [0 -9.81] ...
                               + node(r,c).force_ext;

        end
    end

    % Position, velocity, and acceelleration update
    for r = 1 : row        
        for c = 1: col
            if  node(r,c).isFixed ~= 1            
                node(r,c).acc = node(r,c).force ./ mass;
                node(r,c).vel = node(r,c).vel + node(r,c).acc .* ts;
                node(r,c).pos = node(r,c).pos + node(r,c).vel .* ts;           
            end               
        end
    end
    
    nodes.node = node;
end

%%
function canvas = createCanvas(nodes)
    % Graphic thingy     
    index = 1;
    for c = 1 : nodes.col
        for r = 1 : nodes.row
            canvas(index,:) = nodes.node(r, c).pos;
            index = index + 1;
        end
    end

    canvas_min = min(canvas);
    canvas_max = max(canvas);
    range = canvas_max - canvas_min;

    xlim([canvas_min(1)-range(1) canvas_max(1)+range(1)])
    ylim([canvas_min(2)-range(2)*3 canvas_max(2)+range(2)])
end

%% 
function canvas = drawNodes(S, canvas, nodes, timestamp)    
% Draw the nodes
    index = 1;    
    for c = 1 : nodes.col
        % Vertical line, going down
        for r = nodes.row : -1 : 1
            canvas(index, :) = nodes.node(r, c).pos;
            index = index + 1;
        end

        % Zig-zag line, going up
        for r = 1 : nodes.row
            canvas(index,:) = nodes.node(r,c).pos;
            index = index + 1;
            if (c < nodes.col)
                canvas(index ,:) = nodes.node(r, c + 1).pos;
                index = index + 1;
            end                          
        end

    end

    set(S.h, 'XData', canvas(:,1));
    set(S.h, 'YData', canvas(:,2));    
    set(S.mText,'String', timestamp);

    drawnow;
end

%% 
function [r_, c_] = findClosestNodeFromPt(nodes, pt)
    d = inf;
    for r = 1 : nodes.row
        for c = 1 : nodes.col
            d_ = norm(nodes.node(r,c).pos - pt);
            if (d_ < d)
                d = d_;
                r_ = r;
                c_ = c;
            end
        end
    end
end

%%
function startDragFcn(varargin, S)   
     S = guidata(S.f);  
     set( S.f, 'WindowButtonMotionFcn', @(varargin)draggingFcn(varargin, S) );
     pt = get(S.a, 'CurrentPoint');
     S.pt1 = pt(1,1:2);
     S.pt2 = S.pt1;
     S.k = 1;
     S.trigger = 1;
     guidata(S.f, S);
end

%%
function draggingFcn(varargin, S)
    S = guidata(S.f);  
    pt = get(S.a, 'CurrentPoint');
    S.pt2 = pt(1,1:2);    
    guidata(S.f, S);
end

%%
function stopDragFcn(varargin, S)
    S = guidata(S.f);  
    S.k = 0;
    guidata(S.f, S);
end
