%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%									  %
%	Tel416 Homework 2 - Exercise 6	  %
%									  %
%	Orestis Zekai - 2011030021		  %
%	Petros Toupas - 2011030125		  %
%									  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

% Transition matrix
% Assuming that every LINE is a pmf
P = [1/2 1/8 1/8 1/8 1/8 ; 1/4 1/8 1/16 1/16 1/2 ; 1/4 1/8 1/8 1/4 1/4 ; 1/8 0 1/2 1/4 1/8 ; 0 1/2 1/4 1/4 0];
% Number of samples for each experiment
samples = 100;
% Monte Carlo experiments
MC = 1000;

%% Eigenvalue decomposition
[Q,L] = eig(P');
eigenvalues = diag(L);

% Find eigenvalue=1 and normalize its corresponding eigenvector
for i=1:length(eigenvalues)
	if (norm(eigenvalues(i)-1) < 10^(-6))
		p_inf = Q(:,i)/sum(Q(:,i));
		break;
	end
end

% Entropy of convergence distribution
p_infEntropy = -sum(p_inf.*log2(p_inf));

% Compute the condition entropies and then entropy rate
Xn_Given_1 = -sum(P(1,:).*log2(P(1,:)));
Xn_Given_2 = -sum(P(2,:).*log2(P(2,:)));
Xn_Given_3 = -sum(P(3,:).*log2(P(3,:)));
% For states 4,5 there are two probabilities that are zero so we have to check
Xn_Given_4 = 0;
for j=1:length(P(4,:))
	if (P(4,j) ~= 0)
		Xn_Given_4 = Xn_Given_4 - P(4,j)*log2(P(4,j));
	end
end
Xn_Given_5 = 0;
for j=1:length(P(5,:))
	if (P(5,j) ~= 0)
		Xn_Given_5 = Xn_Given_5 - P(5,j)*log2(P(5,j));
	end
end

entropyRate = p_inf(1)*Xn_Given_1 + p_inf(2)*Xn_Given_2 + p_inf(3)*Xn_Given_3 + p_inf(4)*Xn_Given_4 + p_inf(5)*Xn_Given_5;


%% Markov chain simulation

% a)Huffman encoding based on p_inf
% We are gonna use the following huffman code
%
%    1 -> 00
%	 2 -> 010
%	 3 -> 10
%    4 -> 011
%	 5 -> 11
%
% We are gonna assume that each state has the same probability to start from
startingprobs = [0.2 0.2 0.2 0.2 0.2];

% For 1000 Monte Carlo Experiments
for j=1:MC
	
	% Compute the starting state
	helping = rand(1,1);
	if (helping < 0.2)
		starting = 1;
	elseif (helping < 0.4)
		starting = 2;
	elseif (helping < 0.6)
		starting = 3;
	elseif (helping < 0.8)
		starting = 4;
	else
		starting = 5;
	end	
	
	currentState = starting;
	bitHelp = 1;
	
	% Simulate one chain
	for i=1:samples
		% Compute the next state
		helping = rand(1,1);
		if (helping >= 0 && helping <= P(currentState,1))
			nextState = 1;
		elseif (helping >= P(currentState,1) && helping <= ( P(currentState,1)+P(currentState,2) ))
			nextState = 2;
		elseif (helping >= (P(currentState,1)+P(currentState,2)) && helping <= ( P(currentState,1)+P(currentState,2)+P(currentState,3) ))
			nextState = 3;
		elseif (helping >= (P(currentState,1)+P(currentState,2)+P(currentState,3)) && helping <= ( P(currentState,1)+P(currentState,2)+P(currentState,3)+P(currentState,4) ))
			nextState = 4;
		else
			nextState = 5;
		end

		% Encode current state
		if (currentState == 1)
			coded(bitHelp,j) = 0;
			coded(bitHelp+1,j) = 0;
			bitHelp = bitHelp + 2;
		elseif (currentState == 2)
			coded(bitHelp,j) = 0;
			coded(bitHelp+1,j) = 1;
			coded(bitHelp+2,j) = 0;
			bitHelp = bitHelp + 3;
		elseif (currentState == 3)
			coded(bitHelp,j) = 1;
			coded(bitHelp+1,j) = 0;
			bitHelp = bitHelp + 2;
		elseif (currentState == 4)
			coded(bitHelp,j) = 0;
			coded(bitHelp+1,j) = 1;
			coded(bitHelp+2,j) = 1;
			bitHelp = bitHelp + 3;
		else
			coded(bitHelp,j) = 1;
			coded(bitHelp+1,j) = 1;
			bitHelp = bitHelp + 2;
		end	
		
		% Save and update variables
		x(i,j)=currentState;
		currentState = nextState;	
	end
end

% b)Huffman encoding based on distibution x(n)/x(n-1)
% The matrix of codes we are gonna use is the following
%
%     |  0      00     00    1000   0000  |
%     |  100    010    010   1001   1     |
% P = |  101    0110   011   01     011   |
%     |  110    0111   10    11     01    |
%     |  111    1      11    101    0001  |
%
% We are going to use the same starting probabilities
% Plus we are gonna assume that the state before the starting one is state no 1

% This is gonna help us estimate the average code word length
l = 0;

% For 1000 Monte Carlo experiments
for j=1:MC

	% Starting state
	helping = rand(1,1);
	if (helping < 0.2)
		starting = 1;
	elseif (helping < 0.4)
		starting = 2;
	elseif (helping < 0.6)
		starting = 3;
	elseif (helping < 0.8)
		starting = 4;
	else
		starting = 5;
	end	
	
	currentState = starting;
	previousState = 1;
	bitHelp = 1;

	for i=1:samples
		helping = rand(1,1);
		if (helping >= 0 && helping <= P(currentState,1))
			nextState = 1;
		elseif (helping >= P(currentState,1) && helping <= ( P(currentState,1)+P(currentState,2) ))
			nextState = 2;
		elseif (helping >= (P(currentState,1)+P(currentState,2)) && helping <= ( P(currentState,1)+P(currentState,2)+P(currentState,3) ))
			nextState = 3;
		elseif (helping >= (P(currentState,1)+P(currentState,2)+P(currentState,3)) && helping <= ( P(currentState,1)+P(currentState,2)+P(currentState,3)+P(currentState,4) ))
			nextState = 4;
		else
			nextState = 5;
		end

		% Encode current state based on previous state
		if (currentState == 1)
			if (previousState == 1)
				codedb(bitHelp,j) = 0;
				bitHelp = bitHelp + 1;
				l = l + 1;
			elseif (previousState == 2)
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 0;
				bitHelp = bitHelp + 2;
				l = l + 2;
			elseif (previousState == 3)
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 0;
				bitHelp = bitHelp + 2;
				l = l + 2;
			elseif (previousState == 4)
				codedb(bitHelp,j) = 1;
				codedb(bitHelp+1,j) = 0;
				codedb(bitHelp+2,j) = 0;
				codedb(bitHelp+3,j) = 0;
				bitHelp = bitHelp + 4;
				l = l + 4;
			else
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 0;
				codedb(bitHelp+2,j) = 0;
				codedb(bitHelp+3,j) = 0;
				bitHelp = bitHelp + 4;
				l = l + 4;
			end				
		elseif (currentState == 2)
			if (previousState == 1)
				codedb(bitHelp,j) = 1;
				codedb(bitHelp+1,j) = 0;
				codedb(bitHelp+2,j) = 0;
				bitHelp = bitHelp + 3;
				l = l + 3;
			elseif (previousState == 2)
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 1;
				codedb(bitHelp+2,j) = 0;
				bitHelp = bitHelp + 3;
				l = l + 3;
			elseif (previousState == 3)
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 1;
				codedb(bitHelp+2,j) = 0;
				bitHelp = bitHelp + 3;
				l = l + 3;
			elseif (previousState == 4)
				codedb(bitHelp,j) = 1;
				codedb(bitHelp+1,j) = 0;
				codedb(bitHelp+2,j) = 0;
				codedb(bitHelp+3,j) = 1;
				bitHelp = bitHelp + 4;
				l = l + 4;
			else
				codedb(bitHelp,j) = 1;
				bitHelp = bitHelp + 1;
				l = l + 1;
			end				
		elseif (currentState == 3)
			if (previousState == 1)
				codedb(bitHelp,j) = 1;
				codedb(bitHelp+1,j) = 0;
				codedb(bitHelp+2,j) = 1;
				bitHelp = bitHelp + 3;
				l = l + 3;
			elseif (previousState == 2)
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 1;
				codedb(bitHelp+2,j) = 1;
				codedb(bitHelp+3,j) = 0;
				bitHelp = bitHelp + 4;
				l = l + 4;
			elseif (previousState == 3)
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 1;
				codedb(bitHelp+2,j) = 1;
				bitHelp = bitHelp + 3;
				l = l + 3;
			elseif (previousState == 4)
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 1;
				bitHelp = bitHelp + 2;
				l = l + 2;
			else
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 0;
				codedb(bitHelp+2,j) = 1;
				bitHelp = bitHelp + 3;
				l = l + 3;
			end				
		elseif (currentState == 4)
			if (previousState == 1)
				codedb(bitHelp,j) = 1;
				codedb(bitHelp+1,j) = 1;
				codedb(bitHelp+2,j) = 0;
				bitHelp = bitHelp + 3;
				l = l + 3;
			elseif (previousState == 2)
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 1;
				codedb(bitHelp+2,j) = 1;
				codedb(bitHelp+3,j) = 1;
				bitHelp = bitHelp + 4;
				l = l + 4;
			elseif (previousState == 3)
				codedb(bitHelp,j) = 1;
				codedb(bitHelp+1,j) = 0;
				bitHelp = bitHelp + 2;
				l = l + 2;
			elseif (previousState == 4)
				codedb(bitHelp,j) = 1;
				codedb(bitHelp+1,j) = 1;
				bitHelp = bitHelp + 2;
				l = l + 2;
			else
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 1;
				bitHelp = bitHelp + 2;
				l = l + 2;
			end				
		else
			if (previousState == 1)
				codedb(bitHelp,j) = 1;
				codedb(bitHelp+1,j) = 1;
				codedb(bitHelp+2,j) = 1;
				bitHelp = bitHelp + 3;
				l = l + 3;
			elseif (previousState == 2)
				codedb(bitHelp,j) = 1;
				bitHelp = bitHelp + 1;
				l = l + 1;
			elseif (previousState == 3)
				codedb(bitHelp,j) = 1;
				codedb(bitHelp+1,j) = 1;
				bitHelp = bitHelp + 2;
				l = l + 2;
			elseif (previousState == 4)
				codedb(bitHelp,j) = 1;
				codedb(bitHelp+1,j) = 0;
				codedb(bitHelp+2,j) = 1;
				bitHelp = bitHelp + 3;
				l = l + 3;
			else
				codedb(bitHelp,j) = 0;
				codedb(bitHelp+1,j) = 0;
				codedb(bitHelp+2,j) = 0;
				codedb(bitHelp+3,j) = 1;
				bitHelp = bitHelp + 4;
				l = l + 4;
			end				
		end
		xb(i,j)=currentState;
		previousState = currentState;
		currentState = nextState;	
	end
end

%% Markov chain decoding

% a)Dencoding based on p_inf
for j=1:MC
	bitHelp = 1;
	for i=1:samples
		if (coded(bitHelp,j) == 0)
			bitHelp = bitHelp + 1;
			if (coded(bitHelp,j) == 0)
				decoded(i,j) = 1;
				bitHelp = bitHelp + 1;
			else
				if (coded(bitHelp+1,j) == 1)
					decoded(i,j) = 4;
					bitHelp = bitHelp + 2;
				else
					decoded(i,j) = 2;
					bitHelp = bitHelp + 2;
				end
			end
		else
			bitHelp = bitHelp + 1;
			if (coded(bitHelp,j) == 1)
				decoded(i,j) = 5;
				bitHelp = bitHelp + 1;
			else
				decoded(i,j) = 3;
				bitHelp = bitHelp + 1;
			end
		end		
	end
end

if (decoded == x)
	fprintf('Dencoding based on p_inf successful \n');
end

% Computing mean length of coded / decoded symbol based on p_inf distribution
for j=1:MC 
	la = 0;
	for i=1:samples
		if (x(i,j) == 1)
			la = la + 2;
		elseif (x(i,j) == 2)
			la = la + 3;
		elseif (x(i,j) == 3)
			la = la + 2;
		elseif (x(i,j) == 4)
			la = la + 3;
		else
			la = la + 2;
		end	
	end
	mean_length(j) = la/samples;
end

total = 0;
for i=1:length(mean_length)
	total = total + mean_length(i);
end
meanLength = total/length(mean_length);
text = 'Entropy of p_inf is: %4.4f \n';
fprintf(text,p_infEntropy);
text = 'Mean Lenght based on p_inf coding is: %4.4f \n\n';
fprintf(text,meanLength)


% b)Decoding based on distibution x(n)/x(n-1)
% We assusme that the reciever know the state that the sender starts

for j=1:MC
	state = 1;
	bitHelp = 1;
	for i=1:samples
		if (state == 1)
			if (codedb(bitHelp,j) == 0)
				decodedb(i,j) = 1;
				state = 1;
				bitHelp = bitHelp + 1;
			else
				bitHelp = bitHelp + 1;
				if (codedb(bitHelp,j) == 1)
					bitHelp = bitHelp + 1;
					if (codedb(bitHelp,j) == 0)
						decodedb(i,j) = 4;
						state = 4;
						bitHelp = bitHelp + 1;
					else
						decodedb(i,j) = 5;
						state = 5;
						bitHelp = bitHelp + 1;
					end
				else
					bitHelp = bitHelp + 1;
					if (codedb(bitHelp,j) == 0)
						decodedb(i,j) = 2;
						state = 2;
						bitHelp = bitHelp + 1;
					else
						decodedb(i,j) = 3;
						state = 3;
						bitHelp = bitHelp + 1;
					end
				end
			end
		elseif (state == 2)
			if (codedb(bitHelp,j) == 1)
				decodedb(i,j) = 5;
				state = 5;
				bitHelp = bitHelp + 1;
			else
				bitHelp = bitHelp + 1;
				if (codedb(bitHelp,j) == 0)
					decodedb(i,j) = 1;
					state = 1;
					bitHelp = bitHelp + 1;
				else
					if (codedb(bitHelp+1,j) == 0)
						decodedb(i,j) = 2;
						state = 2;
						bitHelp = bitHelp + 2;
					else
						if (codedb(bitHelp+2,j) == 0)
							decodedb(i,j) = 3;
							state = 3;
							bitHelp = bitHelp + 3;
						else
							decodedb(i,j) = 4;
							state = 4;
							bitHelp = bitHelp + 3;
						end
					end
				end
			end
		elseif (state == 3)
			if (codedb(bitHelp,j) == 1)
				bitHelp = bitHelp + 1;
				if (codedb(bitHelp,j) == 0)
					decodedb(i,j) = 4;
					state = 4;
					bitHelp = bitHelp + 1;
				else
					decodedb(i,j) = 5;
					state = 5;
					bitHelp = bitHelp + 1;
				end
			else
				bitHelp = bitHelp + 1;
				if (codedb(bitHelp,j) == 0)
					decodedb(i,j) = 1;
					state = 1;
					bitHelp = bitHelp + 1;
				else
					bitHelp = bitHelp + 1;
					if (codedb(bitHelp,j) == 0)
						decodedb(i,j) = 2;
						state = 2;
						bitHelp = bitHelp + 1;
					else
						decodedb(i,j) = 3;
						state = 3;
						bitHelp = bitHelp + 1;
					end
				end
			end
		elseif (state == 4)
			if (codedb(bitHelp,j) == 0)
				decodedb(i,j) = 3;
				state = 3;
				bitHelp = bitHelp + 2;
			else
				bitHelp = bitHelp + 1;
				if (codedb(bitHelp,j) == 1)
					decodedb(i,j) = 4;
					state = 4;
					bitHelp = bitHelp + 1;
				else
					bitHelp = bitHelp + 1;
					if (codedb(bitHelp,j) == 1)
						decodedb(i,j) = 5;
						state = 5;
						bitHelp = bitHelp + 1;
					else
						bitHelp = bitHelp + 1;
						if (codedb(bitHelp,j) == 0)
							decodedb(i,j) = 1;
							state = 1;
							bitHelp = bitHelp + 1;
						else
							decodedb(i,j) = 2;
							state = 2;
							bitHelp = bitHelp + 1;
						end
					end
				end
			end
		else
			if (codedb(bitHelp,j) == 1)
				decodedb(i,j) = 2;
				state = 2;
				bitHelp = bitHelp + 1;
			else
				bitHelp = bitHelp + 1;
				if (codedb(bitHelp,j) == 1)
					decodedb(i,j) = 4;
					state = 4;
					bitHelp = bitHelp + 1;
				else
					bitHelp = bitHelp + 1;
					if (codedb(bitHelp,j) == 1)
						decodedb(i,j) = 3;
						state = 3;
						bitHelp = bitHelp + 1;
					else
						bitHelp = bitHelp + 1;
						if (codedb(bitHelp,j) == 1)
							decodedb(i,j) = 5;
							state = 5;
							bitHelp = bitHelp + 1;
						else
							decodedb(i,j) = 1;
							state = 1;
							bitHelp = bitHelp + 1;
						end
					end
				end
			end
		end
	end
end

if (decodedb == xb)
	disp('Decoding based on distibution x(n)/x(n-1) successful');
end

text = 'Entropy rate of probability distribution x(n)/x(n-1) is: %4.4f \n';
fprintf(text,entropyRate);
text = 'Mean Lenght based on probability distribution x(n)/x(n-1) coding is: %4.4f \n\n';
fprintf(text,l/(MC*samples));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&%
%--------------------------------------CONCLUSIONS------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                           %
% Coding with the convergence distibution has the same results as exercise 5. We get really %
% close to entropy, but the probabilities are not a "power of two" so we never actually     %
% reach it.                                                                                 %
%                                                                                           %
% Entropy rate is lesser than entropy in this case.So, taking advantage that the source has %
% memory, when we compute the average code word length we get a number that is of course    %
% greater than entropy rate but lesser than entropy.                                        %
%                                                                                           %
% So making use of the memory of a source we can decrease the average code word length      %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%