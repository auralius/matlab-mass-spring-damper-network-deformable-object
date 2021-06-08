%% Deformable object with interconnected mass-spring-damper
%
%  Author : Auralius Manurung (manurung.auralius@gmail.com)
%  Note   : Collapsing under gravity.
%

%% Data structures for the container
%  particles.initial_pos
%  particles.pos
%  particles.vel
%  particles.acc


%%
% This is the main function
function msd_up_tree()
% Warning! This clears everything!
close all
clc;
clear all;

% Parameters
stiffness = 10;    % N/m
damping = 1;      % Ns/m
mass = 0.1;          % Kg
k_wall = 50000;
ts = 0.001;        % Seconds
t_end = 15;

% Build the model and the canvas
[particles, sorted_idx] = build_model();

N = length(particles);
K = ceil(t_end/ts);
X = zeros(N,2,K);

% This is the main iterations
for i = 1 : K
   particles = update_states(particles, sorted_idx, mass, stiffness, damping, k_wall, ts);
   
    % store the results for animation purpose
    for k = 1:length(particles)
        X(k,:,i) = particles(k).pos;
    end
end

save('msd_up_tree_sim_result.mat', 'X');
GenerateGIF('msd_up_tree.gif', X, 0.1, ts)

end

%%
% Build the container
function [particles, sorted_idx] = build_model()
% load the particles
x = dlmread('up_tree_logo.mat');

% k = 1;
% h = 0.01;
% for px = 0.4 : 2*h : 0.6
%     for py = 0.4 : 2*h : 0.6
%         x(k,:)=[px py];
%         k = k + 1;     
%     end
% end

N = length(x);

for k = 1 : N
    particles(k).initial_pos = x(k,:);
    particles(k).pos = x(k,:);
    particles(k).vel = [0 0];
    particles(k).acc = [0 0];
end

sorted_idx = GetSortedIndex(x);

end

%%
% Update all container per time sampling
function particles = update_states(particles, sorted_idx, mass, stiffness, damping, k_wall, ts)
N = length(particles);
f = zeros(N, 2);
    
for i = 1:N
    % find neighboring particles
    neighbour_idx = FindNeighbours(i, sorted_idx, 3);

    for j = 1 : length(neighbour_idx)
        pair_idx = neighbour_idx(j);
        
        l0 = particles(i).initial_pos - particles(pair_idx).initial_pos;
        lt = particles(i).pos - particles(pair_idx).pos;
        %n = norm(lt, 2);
        %f_ = stiffness * (norm(l0, 2) - n) * lt / n;
        f_ = stiffness * (l0-lt);
        f(i,:) = f(i,:) + f_;
        %f(,:) = f(j,:) - f_;
    end
    
    f_contact = CalculateContactForceOneParticle(particles(i).pos, k_wall);
    
    f(i,:) = f(i,:)  - damping * particles(i).vel + ...
        mass * [0 -9.81] + f_contact;
end


% Position, velocity, and acceelleration update
for i = 1:N
    particles(i).acc = f(i,:) ./ mass;
    particles(i).vel = particles(i).vel + particles(i).acc .* ts;
    particles(i).pos = particles(i).pos + particles(i).vel .* ts;
end

end

%%
function f = CalculateContactForceOneParticle(x, k_wall)
% The walls are located at
% x<0, x>1, and y<0

f = [0 0];


if x(1,1) < 0
    f(1,1) = -k_wall*x(1,1);    % x-potivie force
elseif x(1,1) > 1
    f(1,1) = k_wall*(1-x(1,1)); % x-negative force
end

if x(1,2) < 0
    f(1,2) = -k_wall*x(1,2);    % y-positive force
end
end

%%
function GenerateGIF(fn, X, every_n_secs, data_sampling_time)
% Initial drawing
h_fig = figure;
hold on

h_plot = plot(X(:,1,1), X(:,2,1), 'b.');
axis equal
xlim([0 1]);
ylim([0 1]);

h_text = uicontrol('style','text');

K = size(X , 3);

for k = 1 : K
    if mod(k-1, every_n_secs/data_sampling_time) == 0
        set(h_plot,'XData', X(:,1,k), 'YData', X(:,2,k));
        set(h_text,'String', (k-1)*data_sampling_time);
        drawnow;
        write2gif(h_fig, k, fn);
    end
end
end

%%
function neighbour_idx = FindNeighbours(i, sorted_idx, n_neighbours)
N = length(sorted_idx);

% find where is the location of particles i in the sorted_idx
i_ = find(sorted_idx, i);

% the neighbours are the adjacent particles within a certain range
neighbour_idx = sorted_idx(max(i_-n_neighbours,1) : min(i_+n_neighbours,N));
neighbour_idx = neighbour_idx(neighbour_idx~=i);
end

%%
function sorted_idx = GetSortedIndex(x)
[~,sorted_idx] = sort(vecnorm(x,2,2));
end