close all;
clear all;
clc;

t = [1:10000];      % Define number of time slots
N = 6;              % Define number of stations
pf = 0.1;           % Define probability that a station generates a frame
frSize = 10;        % Define frame duration
maxCol = 16;        % Define maximum times a station will try to retransmitt

nodes = zeros(N,1);    % At t = 0, 0 stations want to transmitt

% Declaring some metric vectors and counters to ensure functionallity
channelState = zeros(length(t),1);
transmittingStation = zeros(length(t),1);
collisionNumber = zeros(N,1);
tryingToTransmitt = zeros(N,1);
tryAgainIn = zeros(N,1);
alreadyTrying = zeros(N,1);
gotNext = zeros(N,1);
ii = 1;
                    
while ii <= length(t)    % For all time slots   
    nodes = zeros(N,1);
    gotNext = zeros(N,1);
    % We generate a probability for every station
    probs = rand(N,1);
    % If the probability is less than pf then the station tries to
    % transmitt. So the corresponding position at the node vector becomes
    % 1. Given that the station is not already trying to transmit
    for jj = 1:N
        if (probs(jj,1) < pf && alreadyTrying(jj) ~= 1)
            nodes(jj) = 1;
        end
    end
    
    % We find how many stations that were waiting will try to transmit in
    % this slot
    for kk = 1:N
        if (alreadyTrying(kk) == 1)     % If this station was trying before, update its waiting time
            if (tryAgainIn(kk) == 0)
                gotNext(kk) = 1;            % Denote intention to transmit
            else
                tryAgainIn(kk) = tryAgainIn(kk) - 1;
            end            
        end
    end
    
    % Find how many nodes want to transmit
    nodesToTransmit = sum(nodes == ones(N,1)) + sum(gotNext == ones(N,1));
    
    % If only one station tries to transmitt, then he gets the channel. If
    % zero stations try to transmit, then the channel is idle. If more than
    % one stations try to transmitt, then we have collision.
    if (nodesToTransmit == 1)           % Case 1: One station want to transmitt
        jj = 1;
        while jj <= frSize
            channelState(ii,1) = 1;     % For 10 time slots the station that got the channel transmitts
            if (sum(nodes) == 1 && sum(gotNext) == 0)
                transmittingStation(ii,1) = find(nodes == 1);       % Find the station transmitting
            elseif (sum(nodes) == 0 && sum(gotNext) == 1)
                transmittingStation(ii,1) = find(gotNext == 1);       % Find the station transmitting
            end
            
            
            for kk = 1:N
                if (alreadyTrying(kk,1) == 1)                   
                    if (kk == transmittingStation(ii,1))        % If station transmitting now was one that had collided before,
                        alreadyTrying(kk,1) = 0;                % update his data.
                        collisionNumber(kk) = 0;
                        tryAgainIn(kk) = 0;
                    end
                end
            end
            
            % Check if someone from before was trying to transmitt. If his
            % time has come, since the slot is full, let him try again later.
            for kk=1:N
                if (alreadyTrying(kk) == 1)
                    tryAgainIn(kk) = tryAgainIn(kk) - 1;        % Update waiting time
                    if (tryAgainIn(kk) <= 0)                    % If waiting time get to zero, station cannot transmitt so waits again
                        collisionNumber(kk) = collisionNumber(kk) + 1;          % Update numbers the station have collided
                        if (collisionNumber(kk) > maxCol)   % If we have more than 16 collisions, station drops the packet
                            tryingToTransmitt(kk,1) = 0;
                            alreadyTrying(kk) = 0;
                            collisionNumber(kk) = 0;
                            continue;
                        end
                        n = min(10,collisionNumber(kk,1));                      % Calculate n
                        tryAgainIn(kk,1) = round(((2^n - 1) - 0)*rand(1,1));    % Generate random waiting period
                    end
                end
            end
            
            % Check if someone new wants to transmitt
            probs = rand(N,1);
            for kk = 1:N            % If he wants again AND he is not already waiting AND he is not transmitting now
                if (probs(kk,1) < pf && alreadyTrying(kk) == 0 && kk ~= transmittingStation(ii,1))
                    tryingToTransmitt(kk,1) = 1;
                end
            end
            
            if (sum(tryingToTransmitt) == 0)        % If noone wants to transmit, continue
                jj = jj + 1;
                ii = ii + 1;
                continue;
            else                                    % One or more stations want to transmitt!!
                for kk = 1:N
                    if (tryingToTransmitt(kk,1) == 1)       % Find who wants to transmitt
                        collisionNumber(kk) = collisionNumber(kk) + 1;          % Update numbers they have collided
                        if (collisionNumber(kk) > maxCol)   % If we have more than 16 collisions, station drops the packet
                            tryingToTransmitt(kk,1) = 0;
                            collisionNumber(kk) = 0;
                            alreadyTrying(kk) = 0;
                            continue;
                        end
                        n = min(10,collisionNumber(kk,1));                      % Calculate n
                        tryAgainIn(kk,1) = round(((2^n - 1) - 0)*rand(1,1));    % Generate random waiting period
                        alreadyTrying(kk) = 1;                                  % Denote that this station is trying to transmitt
                    end
                end
            end
            
            tryingToTransmitt = zeros(N,1);        % Clear values of vector
            
            jj = jj + 1;        % Update counter for this frame transmitted
            ii = ii + 1;        % Update counter holding time slots
        end
    elseif (nodesToTransmit == 0)            % Case 2: Noone wants to transmitt
        channelState(ii,1) = 0;              % For this slot, the channel is idle
        transmittingStation(ii,1) = 0;       % No station is transmitting
        
        ii = ii + 1;                         % Update counter holding time slots
        
    else                                     % Case 3: More than one station wants to transmitt -> collision
        channelState(ii,1) = 2;              % For this slot, the channel is not used, because of collision
        transmittingStation(ii,1) = 0;       % No station is transmitting
        
        % Find the stations that collided and make them wait
        for kk = 1:N        % Case 1: A station that was waiting, tried to transmit and collided
            if (alreadyTrying(kk) == 1 && tryAgainIn(kk) == 0)     % If this station was waiting to transmit earlier and his waiting time ended
                collisionNumber(kk) = collisionNumber(kk) + 1;     % Update number of collisions
                if (collisionNumber(kk) > maxCol)                  % If we have more than 16 collisions, station drops the packet
                    collisionNumber(kk) = 0;
                    alreadyTrying(kk) = 0;
                    continue;
                end
                n = min(10,collisionNumber(kk,1));                      % Calculate n
                tryAgainIn(kk,1) = round(((2^n - 1) - 0)*rand(1,1));    % Generate random waiting period and make station wait again
            end
        end
        
        for kk = 1:N        % Case 2: One or more of the stations tried to transmitt without waiting first and collided
            if (nodes(kk) == 1)     % If this station tried to transmitt, make it wait
                collisionNumber(kk) = collisionNumber(kk) + 1;          % In this case, we don't need to check the number of collisions
                n = min(10,collisionNumber(kk));                        % Calculate n
                tryAgainIn(kk,1) = round(((2^n - 1) - 0)*rand(1,1));    % Generate random waiting period
                alreadyTrying(kk) = 1;                                  % Denote that this station is trying to transmitt
            end
        end
        
       ii = ii + 1;         % Update counter holding time slots
    end
end