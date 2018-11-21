%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%									  %
%	Tel415 Homework 3 				  %
%									  %
%	Orestis Zekai - 2011030021		  %
%	Petros Toupas - 2011030125		  %
%									  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% Programm parameters

% Bits transmited
N = 10^4;
% Number of users
n = 10;
% Signature length
m = 31;
% Power of users
% Zero represents the power of user 1.It changes after each iteration.
P = [0 5 5 5 10 10 10 15 15 15];
% Power interval for user 1
p_user1 = [0:20];
% Samples used to determine SMI and SA-SMI filters
K = 200;

% Generate random bits for each user
for i=1:n
	for j=1:N
		helping = rand();
		if (helping <= 0.5)
			B(i,j) = -1;
		else
			B(i,j) = 1;
		end
	end
end

% Create signature for each user
value = 1/sqrt(m);
for i=1:m
	for j=1:n
		helping = rand();
		if (helping <= 0.5)
			s(i,j) = -value;
		else
			s(i,j) = value;
		end
	end
end

% Vector of zeros for BERs
BER_MF = zeros(1,length(p_user1));
BER_MVDR = zeros(1,length(p_user1));
BER_ZF = zeros(1,length(p_user1));
BER_SMI = zeros(1,length(p_user1));
BER_SA_SMI = zeros(1,length(p_user1));
BER_SA_SMI_loop = zeros(1,length(p_user1));

% Random noise
noise_matrix = randn(m,N);

%% Compute BER for each value of p_user1

for i=1:length(p_user1)

	% Change value of the 1st place of power vector
	P(1) = p_user1(i);
	Watts = db2pow(P);

	% Power matrix creation
	power_matrix = diag(sqrt(Watts));

	% Transmit and add noise
	x = s*power_matrix*B + noise_matrix;

	% MF equals with the signature of the 1st user
	w_matched = s(:,1);
	% Compute MVDR filter
	R = s*diag(Watts)*s.' + eye(m);
	w_mvdr = inv(R)*s(:,1);
	% Compute ZF filter
	zf_matrix = (s*(inv(s'*s)));
	w_zf = zf_matrix(:,1);
	% Compute SMI filter
	smi_matrix = (x(:,1:K)*x(:,1:K)')./K;
	w_smi = inv(smi_matrix)*s(:,1);
	
	% Compute SA-SMI filter
	x_new = x - s(:,1)*sqrt(Watts(1))*B(1,:);
 	sa_smi_matrix = (x_new(:,1:K) * (x_new(:,1:K))')./K;
    w_sa_smi = inv(sa_smi_matrix)*s(:,1);
	
	% Compute SA-SMI loop filter
	R_y = (x(:,1:K)*x(:,1:K)')./K;
	w_loop = inv(R_y)*s(:,1);
	
	x_loop = w_loop'*x;
	est_loop = sign(x_loop);
	while(1)
		
		z = x - s(:,1)*sqrt(Watts(1))*est_loop;

		R_z = (z(:,1:K)*z(:,1:K)')./K;
		w_loop = inv(R_z)*s(:,1);

		x_loop2 = w_loop'*x;
		est_loop2 = sign(x_loop2);

		if (est_loop == est_loop2)
			break;
		else
			est_loop = est_loop2;
		end
	end

	% Filter the signal with matched filter
	x_matched = w_matched.'*x;
	% Filter the signal with MVDR filter
	x_mvdr = w_mvdr.'*x;
	% Filter the signal with ZF filter
	x_zf = w_zf.'*x;
	% % Filter the signal with SMI filter
	x_smi = w_smi.'*x;
	% % Filter the signal with SA-SMI filter
	x_sa_smi = w_sa_smi.'*x;

	% Make a decision based on x_matched
	est_matched = sign(x_matched);
	% Make a decision based on x_mvdr
	est_mvdr = sign(x_mvdr);
	% Make a decision based on x_zf
	est_zf = sign(x_zf);
	% Make a decision based on x_smi
	est_smi = sign(x_smi);
	% Make a decision based on x_sa_smi
	est_sa_smi = sign(x_sa_smi);



	for k=1:N
		if (norm(est_matched(k) - B(1,k)) >= 10^(-6))
			BER_MF(i) = BER_MF(i) + 1;
		end
		if (norm(est_mvdr(k) - B(1,k)) >= 10^(-6))
			BER_MVDR(i) = BER_MVDR(i) + 1;
		end
		if (norm(est_zf(k) - B(1,k)) >= 10^(-6))
			BER_ZF(i) = BER_ZF(i) + 1;
		end
		if (norm(est_smi(k) - B(1,k)) >= 10^(-6))
			BER_SMI(i) = BER_SMI(i) + 1;
		end
		if (norm(est_sa_smi(k) - B(1,k)) >= 10^(-6))
			BER_SA_SMI(i) = BER_SA_SMI(i) + 1;
		end
		if (norm(est_loop2(k) - B(1,k)) >= 10^(-6))
			BER_SA_SMI_loop(i) = BER_SA_SMI_loop(i) + 1;
		end
	end
	BER_MF(i) = BER_MF(i)/N;
	BER_MVDR(i) = BER_MVDR(i)/N;
	BER_ZF(i) = BER_ZF(i)/N;
	BER_SMI(i) = BER_SMI(i)/N;
	BER_SA_SMI(i) = BER_SA_SMI(i)/N;
	BER_SA_SMI_loop(i) = BER_SA_SMI_loop(i)/N;
	
end

figure()
semilogy(p_user1,BER_MF,'b');
hold on;
semilogy(p_user1,BER_MVDR,'r');
hold on;
semilogy(p_user1,BER_ZF,'g');
hold on;
semilogy(p_user1,BER_SMI,'y');
hold on;
semilogy(p_user1,BER_SA_SMI,'c');
hold on;
semilogy(p_user1,BER_SA_SMI_loop,'k');
hold on;
xlabel('SNR');
ylabel('BER');
legend('MF','MVDR','ZF','SMI','SA-SMI','SA_SMI_loop');


%% Optimize signature for user 1
% Rz = R - Watts(1)*s(:,1)*s(:,1).';

% [Q,L] = eig(Rz);

% s(:,1) = Q(:,1)./sum(Q(:,1));

%% Compute BERs for new signature

% Vector of zeros for BERs
BER_MF = zeros(1,length(p_user1));
BER_MVDR = zeros(1,length(p_user1));
BER_ZF = zeros(1,length(p_user1));
BER_SMI = zeros(1,length(p_user1));
BER_SA_SMI = zeros(1,length(p_user1));
BER_SA_SMI_loop = zeros(1,length(p_user1));

% Random noise
noise_matrix = randn(m,N);

% Noise variance
variance = var(var(noise_matrix));

%% Compute BER for each value of p_user1

for i=1:length(p_user1)

	% Change value of the 1st place of power vector
	P(1) = p_user1(i);
	Watts = db2pow(P);

	R = s*diag(Watts)*s.' + eye(m);
	R_sign = R - sqrt(Watts(1))*s(:,1)*s(:,1).';
	[Q,L] = eig(R_sign);
	s(:,1) = Q(:,1)./sum(Q(:,1));

	% Power matrix creation
	power_matrix = diag(sqrt(Watts));

	% Transmit and add noise
	x = s*power_matrix*B + noise_matrix;

	% MF equals with the signature of the 1st user
	w_matched = s(:,1);
	% Compute MVDR filter
	R = s*diag(Watts)*s.' + eye(m);
	w_mvdr = inv(R)*s(:,1);
	% Compute ZF filter
	zf_matrix = (s*(inv(s'*s)));
	w_zf = zf_matrix(:,1);
	% Compute SMI filter
	smi_matrix = (x(:,1:K)*x(:,1:K)')./K;
	w_smi = inv(smi_matrix)*s(:,1);
	
	% Compute SA-SMI filter
	x_new = x - s(:,1)*sqrt(Watts(1))*B(1,:);
 	sa_smi_matrix = (x_new(:,1:K) * (x_new(:,1:K))')./K;
    w_sa_smi = inv(sa_smi_matrix)*s(:,1);
	
	% Compute SA-SMI loop filter
	R_y = (x(:,1:K)*x(:,1:K)')./K;
	w_loop = inv(R_y)*s(:,1);
	
	x_loop = w_loop'*x;
	est_loop = sign(x_loop);
	while(1)
		
		z = x - s(:,1)*sqrt(Watts(1))*est_loop;

		R_z = (z(:,1:K)*z(:,1:K)')./K;
		w_loop = inv(R_z)*s(:,1);

		x_loop2 = w_loop'*x;
		est_loop2 = sign(x_loop2);

		if (est_loop == est_loop2)
			break;
		else
			est_loop = est_loop2;
		end
	end

	% Filter the signal with matched filter
	x_matched = w_matched.'*x;
	% Filter the signal with MVDR filter
	x_mvdr = w_mvdr.'*x;
	% Filter the signal with ZF filter
	x_zf = w_zf.'*x;
	% % Filter the signal with SMI filter
	x_smi = w_smi.'*x;
	% % Filter the signal with SA-SMI filter
	x_sa_smi = w_sa_smi.'*x;

	% Make a decision based on x_matched
	est_matched = sign(x_matched);
	% Make a decision based on x_mvdr
	est_mvdr = sign(x_mvdr);
	% Make a decision based on x_zf
	est_zf = sign(x_zf);
	% Make a decision based on x_smi
	est_smi = sign(x_smi);
	% Make a decision based on x_sa_smi
	est_sa_smi = sign(x_sa_smi);



	for k=1:N
		if (norm(est_matched(k) - B(1,k)) >= 10^(-6))
			BER_MF(i) = BER_MF(i) + 1;
		end
		if (norm(est_mvdr(k) - B(1,k)) >= 10^(-6))
			BER_MVDR(i) = BER_MVDR(i) + 1;
		end
		if (norm(est_zf(k) - B(1,k)) >= 10^(-6))
			BER_ZF(i) = BER_ZF(i) + 1;
		end
		if (norm(est_smi(k) - B(1,k)) >= 10^(-6))
			BER_SMI(i) = BER_SMI(i) + 1;
		end
		if (norm(est_sa_smi(k) - B(1,k)) >= 10^(-6))
			BER_SA_SMI(i) = BER_SA_SMI(i) + 1;
		end
		if (norm(est_loop2(k) - B(1,k)) >= 10^(-6))
			BER_SA_SMI_loop(i) = BER_SA_SMI_loop(i) + 1;
		end
	end
	BER_MF(i) = BER_MF(i)/N;
	BER_MVDR(i) = BER_MVDR(i)/N;
	BER_ZF(i) = BER_ZF(i)/N;
	BER_SMI(i) = BER_SMI(i)/N;
	BER_SA_SMI(i) = BER_SA_SMI(i)/N;
	BER_SA_SMI_loop(i) = BER_SA_SMI_loop(i)/N;
	
end

% figure();
semilogy(p_user1,BER_MF,'--b');
hold on;
semilogy(p_user1,BER_MVDR,'--mo');
hold on;
% semilogy(p_user1,BER_ZF,'--g');
% hold on;
% semilogy(p_user1,BER_SMI,'--y');
% hold on;
% semilogy(p_user1,BER_SA_SMI,'--c');
% hold on;
% semilogy(p_user1,BER_SA_SMI_loop,'--k');
hold off;
xlabel('SNR');
ylabel('BER');
legend('MF','MVDR','ZF','SMI','SA-SMI','SA_SMI_loop');