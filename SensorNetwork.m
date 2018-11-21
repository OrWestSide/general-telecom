close all;
clear all;
clc;

N = 1000;                     % Define constant N
sensors = [1:N]';             % Create sensor vector   
region = zeros(N,N);          % Create region with zeros
iter = 50;                    % Number of iterations
a = 1;                        % Min range of sensor position
b = 1000;                     % Max range of sensor position
R = 50;                       % Sensors range
old_cost = zeros(iter,1);
new_cost = zeros(iter,1);

% Place randomly 4 sink nodes
for ii = 1:4
    x_sink(ii,1) = round((b-a).*rand + a);
    y_sink(ii,1) = round((b-a).*rand + a);
end
% Store positions
xs = x_sink;
ys = y_sink;

for i = 1:iter
    disp(i)
    x_sink = xs;
    y_sink = ys;
    % Find a random coordinate in 1000x1000 plex and assign a sensor
    for ii = 1:N
        x(ii,1) = round((b-a).*rand + a);
        y(ii,1) = round((b-a).*rand + a);
        region(x,y) = sensors(ii);
    end
    % Plot sensors in area
%     figure()
%     plot(x,y,'ks')
%     text(x,y,num2str(sensors))
%     title('Randomly placed nodes in the 1000x1000 grid')
    % Find nodes that can be connected
    neighbours = zeros(N,N);      % Initialize neighbour array
    for ii = 1:N
        for jj = 1:N
            if (ii ~= jj)
                dist = sqrt( (x(ii)-x(jj))^2 + (y(ii)-y(jj))^2);
                if (dist < R)
                    neighbours(ii,jj) = 1;
                end
            end
        end
    end
%     % Draw dots between connected sensors
%     figure()
%     plot(x,y,'ks')
%     text(x,y,num2str(sensors))
%     hold on;
%     for ii = 1:N
%         for jj = 1:N
%             if (neighbours(ii,jj) == 1)
%                 plot([x(ii) x(jj)],[y(ii) y(jj)],'b:');
%                 hold on;
%             end
%         end
%     end
%     title('After random placement of sink nodes')

    % Connect extra nodes to rest of network
    for ii = 1:4
        for jj = 1:N
            dist = sqrt( (x_sink(ii)-x(jj))^2 + (y_sink(ii)-y(jj))^2);
            if (dist < R)
                neighbours(ii+N,jj) = 1;
                neighbours(jj,ii+N) = 1;
            else
                neighbours(ii+N,jj) = 0;
                neighbours(jj,ii+N) = 0;
            end
        end
    end

    % Conver to graph.object, for easier manipulation
    G = graph(neighbours);

    costs = 0;      % Total cost
    penalites = 0;  % Total penalties
    P = 10;         % Penalty cost
    
    for ii = 1:N
        C = 0;      % Cost for specific node
        % Compute paths for every sink node
        path1 = shortestpath(G,ii,N+1);
        path2 = shortestpath(G,ii,N+2);
        path3 = shortestpath(G,ii,N+3);
        path4 = shortestpath(G,ii,N+4);

        % See it there is at least one path
        if (length(path1) > 0 || length(path2) > 0 || length(path3) > 0 || length(path4) > 0)
            % Find shortest path
            min_path = 0;
            min_length = 20000;
            if (length(path1) < min_length && length(path1) > 0)
                min_path = path1;
                min_length = length(path1);
            end
            if (length(path2) < min_length && length(path2) > 0)
                min_path = path2;
                min_length = length(path2);
            end
            if (length(path3) < min_length && length(path3) > 0)
                min_path = path3;
                min_length = length(path3);
            end
            if (length(path4) < min_length && length(path4) > 0)
                min_path = path4;
                min_length = length(path4);
            end
            C = min_length-1;
        else
            penalties = penalites + P;
        end
        costs = costs + C;

        % Plot path for some nodes
%         if (ii == 44)
%             % Start figure and plot existing nodes
%             figure()
%             plot(x,y,'ks');
%             hold on;
%             text(x,y,num2str(sensors));
%             plot(x_sink,y_sink,'rs');
%             for iii = 1:N
%                 for jjj = 1:N
%                     if (neighbours(iii,jjj) == 1)
%                         plot([x(iii) x(jjj)],[y(iii) y(jjj)],'b:');
%                         hold on;
%                     end
%                 end
%             end
%             % Plot path
%             for lll = 1:min_length-1
%                 if (lll == min_length-1)
%                     plot([x(min_path(lll)) x_sink(mod(min_path(lll+1),100))],[y(min_path(lll)) y_sink(mod(min_path(lll+1),100))],'k');
%                     hold on;
%                 else
%                     plot([x(min_path(lll)) x(min_path(lll+1))],[y(min_path(lll)) y(min_path(lll+1))],'k');
%                     hold on;
%                 end
%             end
%             hold off;
%             title('Shortest path from node number 44 to closed sink node');
%         end
    end
    old_cost(i) = costs + penalites;
    costs = 0;
    penalites = 0;
    
    % The optimal position of the sink nodes, independently of the position of
    % the position of the rest of the nodes is so that they cover equal area of
    % the plex. So we have:
    x_sink(1) = 250;    % Sink node 1
    y_sink(1) = 250;
    x_sink(2) = 250;    % Sink node 2
    y_sink(2) = 750;
    x_sink(3) = 750;    % Sink node 3
    y_sink(3) = 250;
    x_sink(4) = 750;    % Sink node 4
    y_sink(4) = 750;

    % Plot again, but this time, sink nodes will be optimal
%     figure()
%     plot(x,y,'ks')
%     text(x,y,num2str(sensors))
%     hold on;
%     plot(x_sink,y_sink,'rs')
%     hold on;
%     for ii = 1:N
%         for jj = 1:N
%             if (neighbours(ii,jj) == 1)
%                 plot([x(ii) x(jj)],[y(ii) y(jj)],'b:');
%                 hold on;
%             end
%         end
%     end
%     hold off;
%     title('New sink node placement')
    % Find new neighbours of the new sink nodes replacing the old ones
    for ii = 1:4
        for jj = 1:N
            dist = sqrt( (x_sink(ii)-x(jj))^2 + (y_sink(ii)-y(jj))^2);
            if (dist < R)
                neighbours(ii+N,jj) = 1;
                neighbours(jj,ii+N) = 1;
            else
                neighbours(ii+N,jj) = 0;
                neighbours(jj,ii+N) = 0;
            end
        end
    end

    % Create new graph.object with new sink node placement
    G_new = graph(neighbours);

    costs_new = 0;      % Total cost
    penalites_new = 0;  % Total penalties
    P = 10;         % Penalty cost
    for ii = 1:N
        C = 0;      % Cost for specific node
        % Compute paths for every sink node
        path1 = shortestpath(G_new,ii,N+1);
        path2 = shortestpath(G_new,ii,N+2);
        path3 = shortestpath(G_new,ii,N+3);
        path4 = shortestpath(G_new,ii,N+4);

        % See it there is at least one path
        if (length(path1) > 0 || length(path2) > 0 || length(path3) > 0 || length(path4) > 0)
            % Find shortest path
            min_path = 0;
            min_length = 20000;
            if (length(path1) < min_length && length(path1) > 0)
                min_path = path1;
                min_length = length(path1);
            end
            if (length(path2) < min_length && length(path2) > 0)
                min_path = path2;
                min_length = length(path2);
            end
            if (length(path3) < min_length && length(path3) > 0)
                min_path = path3;
                min_length = length(path3);
            end
            if (length(path4) < min_length && length(path4) > 0)
                min_path = path4;
                min_length = length(path4);
            end
            C = min_length-1;
        else
            penalties_new = penalites_new + P;
        end
        costs_new = costs_new + C;
    end
    new_cost(i) = costs_new + penalites_new;
    costs_new = 0;
    penalites_new = 0;
end


% Display costs
disp('Average old cost')
disp(sum(old_cost)/iter)
disp('Average new cost')
disp(sum(new_cost)/iter)
disp('Improvement')
disp((sum(old_cost)/iter) - (sum(new_cost)/iter))