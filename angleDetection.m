%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%									  %
%	Tel415 Homework 4 				  %
%									  %
%	Orestis Zekai - 2011030021		  %
%	Petros Toupas - 2011030125		  %
%									  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% 1a

P = [7 5 7 5 7];		% Powers in db
P_watt = db2pow(P);		% Convert db in watts
m = 20;					% Elements of antenna array

% Create portions a(theta)

for i=0:(m-1)
	a_theta1(i+1,1) = exp(-1j*pi*i*sin(deg2rad(-75)));
	a_theta2(i+1,1) = exp(-1j*pi*i*sin(deg2rad(-60)));
	a_theta3(i+1,1) = exp(-1j*pi*i*sin(deg2rad(0)));
	a_theta4(i+1,1) = exp(-1j*pi*i*sin(deg2rad(60)));
	a_theta5(i+1,1) = exp(-1j*pi*i*sin(deg2rad(75)));
end
a_theta3 = complex(a_theta3,0);


% Compute autocorrelation matrix
R = P_watt(1)*(a_theta1*a_theta1') + P_watt(2)*(a_theta2*a_theta2') + P_watt(3)*(a_theta3*a_theta3') + P_watt(4)*(a_theta4*a_theta4') + P_watt(5)*(a_theta5*a_theta5') + eye(m);

% Singular value decomposition
[U,S,V] = svd(R);
Qn = U(:,6:20);

angles = [-90:1:90];

for i=1:length(angles)
	
	% Compute a(theta)
	for ii=0:(m-1)
		a_theta(ii+1,1) = exp(-1j*ii*pi*sin(deg2rad(angles(i))));
	end

	P_mf(i) = real(a_theta'*R*a_theta);
	P_mvdr(i) = real(1/(a_theta'*inv(R)*a_theta));
	P_music(i) = 1/(norm(Qn'*a_theta).^2);

end

P_mf = P_mf./max(P_mf);
P_mvdr = P_mvdr./max(P_mvdr);
P_music = P_music./max(P_music);

% ESPIRIT method
var = diag(S);
var = mean(var(6:m));
RR = (R-var*eye(m));
R_esp = [RR(1:m-1,:);RR(2:m,:)];

[Ue,Se,Ve] = svd(R_esp);
U1 = Ue(1:38/2,1:5);
U2 = Ue(38/2+1:38,1:5);

R1 = U1'*U1;
R2 = U1'*U2;

[P,D] = eig(inv(R1)*R2);

w = diag(D);
omega = -1j*log(((w)));
esp_th = real(asind(-omega/pi));


figure();
plot(angles,P_mf,'b');
hold on;
plot(angles,P_mvdr,'r');
hold on;
plot(angles,P_music,'m');
hold on;
line([esp_th(1) esp_th(1)],[0 1],'Color','k');
line([esp_th(2) esp_th(2)],[0 1],'Color','k');
line([esp_th(3) esp_th(3)],[0 1],'Color','k');
line([esp_th(4) esp_th(4)],[0 1],'Color','k');
line([esp_th(5) esp_th(5)],[0 1],'Color','k');
plot([75 75],[0 1],'g--');
hold on;
plot([-75 -75],[0 1],'g--');
hold on;
plot([60 60],[0 1],'g--');
hold on;
plot([-60 -60],[0 1],'g--');
hold on;
plot([0 0],[0 1],'g--');
hold on;
legend('MF','MVDR','MUSIC','ESPRIT','true');
title('P=[7 5 7 5 7] R known');
xlabel('angle');
ylabel('P(angle)');
axis([-90 90 0 1]);

%% 1b
N = 50;

A = [a_theta1 a_theta2 a_theta3 a_theta4 a_theta5];
Po = diag(sqrt(P_watt));
B = sign(randn(5,N));
W = (randn(m,N)+1j*randn(m,N))/sqrt(2);

Y = A*Po*B+W;
R = (1/N).*(Y*Y');

% Singular value decomposition
[U,S,V] = svd(R);
Qn = U(:,6:20);

angles = [-90:1:90];

for i=1:length(angles)
	
	% Compute a(theta)
	for ii=0:(m-1)
		a_theta(ii+1,1) = exp(-1j*ii*pi*sin(deg2rad(angles(i))));
	end

	P_mf(i) = real(a_theta'*R*a_theta);
	P_mvdr(i) = real(1/(a_theta'*inv(R)*a_theta));
	P_music(i) = 1/(norm(Qn'*a_theta).^2);

end

P_mf = P_mf./max(P_mf);
P_mvdr = P_mvdr./max(P_mvdr);
P_music = P_music./max(P_music);

% ESPIRIT method
var = diag(S);
var = mean(var(6:m));
RR = (R-var*eye(m));
R_esp = [RR(1:m-1,:);RR(2:m,:)];

[Ue,Se,Ve] = svd(R_esp);
U1 = Ue(1:38/2,1:5);
U2 = Ue(38/2+1:38,1:5);

R1 = U1'*U1;
R2 = U1'*U2;

[P,D] = eig(inv(R1)*R2);

w = diag(D);
omega = -1j*log(((w)));
esp_th = real(asind(-omega/pi));


figure();
plot(angles,P_mf,'b');
hold on;
plot(angles,P_mvdr,'r');
hold on;
plot(angles,P_music,'m');
hold on;
line([esp_th(1) esp_th(1)],[0 1],'Color','k');
line([esp_th(2) esp_th(2)],[0 1],'Color','k');
line([esp_th(3) esp_th(3)],[0 1],'Color','k');
line([esp_th(4) esp_th(4)],[0 1],'Color','k');
line([esp_th(5) esp_th(5)],[0 1],'Color','k');
plot([75 75],[0 1],'g--');
hold on;
plot([-75 -75],[0 1],'g--');
hold on;
plot([60 60],[0 1],'g--');
hold on;
plot([-60 -60],[0 1],'g--');
hold on;
plot([0 0],[0 1],'g--');
hold on;
legend('MF','MVDR','MUSIC','ESPRIT','true');
title('P=[7 5 7 5 7] N=50');
xlabel('angle');
ylabel('P(angle)');
axis([-90 90 0 1]);

%% 1c
N = 200;

A = [a_theta1 a_theta2 a_theta3 a_theta4 a_theta5];
Po = diag(sqrt(P_watt));
B = sign(randn(5,N));
W = (randn(m,N)+1j*randn(m,N))/sqrt(2);

Y = A*Po*B+W;
R = (1/N).*(Y*Y');

% Singular value decomposition
[U,S,V] = svd(R);
Qn = U(:,6:20);

angles = [-90:1:90];

for i=1:length(angles)
	
	% Compute a(theta)
	for ii=0:(m-1)
		a_theta(ii+1,1) = exp(-1j*ii*pi*sin(deg2rad(angles(i))));
	end

	P_mf(i) = real(a_theta'*R*a_theta);
	P_mvdr(i) = real(1/(a_theta'*inv(R)*a_theta));
	P_music(i) = 1/(norm(Qn'*a_theta).^2);

end

P_mf = P_mf./max(P_mf);
P_mvdr = P_mvdr./max(P_mvdr);
P_music = P_music./max(P_music);

% ESPIRIT method
var = diag(S);
var = mean(var(6:m));
RR = (R-var*eye(m));
R_esp = [RR(1:m-1,:);RR(2:m,:)];

[Ue,Se,Ve] = svd(R_esp);
U1 = Ue(1:38/2,1:5);
U2 = Ue(38/2+1:38,1:5);

R1 = U1'*U1;
R2 = U1'*U2;

[P,D] = eig(inv(R1)*R2);

w = diag(D);
omega = -1j*log(((w)));
esp_th = real(asind(-omega/pi));


figure();
plot(angles,P_mf,'b');
hold on;
plot(angles,P_mvdr,'r');
hold on;
plot(angles,P_music,'m');
hold on;
line([esp_th(1) esp_th(1)],[0 1],'Color','k');
line([esp_th(2) esp_th(2)],[0 1],'Color','k');
line([esp_th(3) esp_th(3)],[0 1],'Color','k');
line([esp_th(4) esp_th(4)],[0 1],'Color','k');
line([esp_th(5) esp_th(5)],[0 1],'Color','k');
plot([75 75],[0 1],'g--');
hold on;
plot([-75 -75],[0 1],'g--');
hold on;
plot([60 60],[0 1],'g--');
hold on;
plot([-60 -60],[0 1],'g--');
hold on;
plot([0 0],[0 1],'g--');
hold on;
legend('MF','MVDR','MUSIC','ESPRIT','true');
title('P=[7 5 7 5 7] N=200');
xlabel('angle');
ylabel('P(angle)');
axis([-90 90 0 1]);

%% 1d
N = 1000;

A = [a_theta1 a_theta2 a_theta3 a_theta4 a_theta5];
Po = diag(sqrt(P_watt));
B = sign(randn(5,N));
W = (randn(m,N)+1j*randn(m,N))/sqrt(2);

Y = A*Po*B+W;
R = (1/N).*(Y*Y');

% Singular value decomposition
[U,S,V] = svd(R);
Qn = U(:,6:20);

angles = [-90:1:90];

for i=1:length(angles)
	
	% Compute a(theta)
	for ii=0:(m-1)
		a_theta(ii+1,1) = exp(-1j*ii*pi*sin(deg2rad(angles(i))));
	end

	P_mf(i) = real(a_theta'*R*a_theta);
	P_mvdr(i) = real(1/(a_theta'*inv(R)*a_theta));
	P_music(i) = 1/(norm(Qn'*a_theta).^2);

end

P_mf = P_mf./max(P_mf);
P_mvdr = P_mvdr./max(P_mvdr);
P_music = P_music./max(P_music);

% ESPIRIT method
var = diag(S);
var = mean(var(6:m));
RR = (R-var*eye(m));
R_esp = [RR(1:m-1,:);RR(2:m,:)];

[Ue,Se,Ve] = svd(R_esp);
U1 = Ue(1:38/2,1:5);
U2 = Ue(38/2+1:38,1:5);

R1 = U1'*U1;
R2 = U1'*U2;

[P,D] = eig(inv(R1)*R2);

w = diag(D);
omega = -1j*log(((w)));
esp_th = real(asind(-omega/pi));


figure();
plot(angles,P_mf,'b');
hold on;
plot(angles,P_mvdr,'r');
hold on;
plot(angles,P_music,'m');
hold on;
line([esp_th(1) esp_th(1)],[0 1],'Color','k');
line([esp_th(2) esp_th(2)],[0 1],'Color','k');
line([esp_th(3) esp_th(3)],[0 1],'Color','k');
line([esp_th(4) esp_th(4)],[0 1],'Color','k');
line([esp_th(5) esp_th(5)],[0 1],'Color','k');
plot([75 75],[0 1],'g--');
hold on;
plot([-75 -75],[0 1],'g--');
hold on;
plot([60 60],[0 1],'g--');
hold on;
plot([-60 -60],[0 1],'g--');
hold on;
plot([0 0],[0 1],'g--');
hold on;
legend('MF','MVDR','MUSIC','ESPRIT','true');
title('P=[7 5 7 5 7] N=1000');
xlabel('angle');
ylabel('P(angle)');
axis([-90 90 0 1]);










%% 2

P = [-7 -5 -7 -5 -7];		% Powers in db
P_watt = db2pow(P);		% Convert db in watts
m = 20;					% Elements of antenna array

% Create portions a(theta)

for i=0:(m-1)
	a_theta1(i+1,1) = exp(-1j*pi*i*sin(deg2rad(-75)));
	a_theta2(i+1,1) = exp(-1j*pi*i*sin(deg2rad(-60)));
	a_theta3(i+1,1) = exp(-1j*pi*i*sin(deg2rad(0)));
	a_theta4(i+1,1) = exp(-1j*pi*i*sin(deg2rad(60)));
	a_theta5(i+1,1) = exp(-1j*pi*i*sin(deg2rad(75)));
end
a_theta3 = complex(a_theta3,0);


% Compute autocorrelation matrix
R = P_watt(1)*(a_theta1*a_theta1') + P_watt(2)*(a_theta2*a_theta2') + P_watt(3)*(a_theta3*a_theta3') + P_watt(4)*(a_theta4*a_theta4') + P_watt(5)*(a_theta5*a_theta5') + eye(m);

% Singular value decomposition
[U,S,V] = svd(R);
Qn = U(:,6:20);

angles = [-90:1:90];

for i=1:length(angles)
	
	% Compute a(theta)
	for ii=0:(m-1)
		a_theta(ii+1,1) = exp(-1j*ii*pi*sin(deg2rad(angles(i))));
	end

	P_mf(i) = real(a_theta'*R*a_theta);
	P_mvdr(i) = real(1/(a_theta'*inv(R)*a_theta));
	P_music(i) = 1/(norm(Qn'*a_theta).^2);

end

P_mf = P_mf./max(P_mf);
P_mvdr = P_mvdr./max(P_mvdr);
P_music = P_music./max(P_music);

% ESPIRIT method
var = diag(S);
var = mean(var(6:m));
RR = (R-var*eye(m));
R_esp = [RR(1:m-1,:);RR(2:m,:)];

[Ue,Se,Ve] = svd(R_esp);
U1 = Ue(1:38/2,1:5);
U2 = Ue(38/2+1:38,1:5);

R1 = U1'*U1;
R2 = U1'*U2;

[P,D] = eig(inv(R1)*R2);

w = diag(D);
omega = -1j*log(((w)));
esp_th = real(asind(-omega/pi));


figure();
plot(angles,P_mf,'b');
hold on;
plot(angles,P_mvdr,'r');
hold on;
plot(angles,P_music,'m');
hold on;
line([esp_th(1) esp_th(1)],[0 1],'Color','k');
line([esp_th(2) esp_th(2)],[0 1],'Color','k');
line([esp_th(3) esp_th(3)],[0 1],'Color','k');
line([esp_th(4) esp_th(4)],[0 1],'Color','k');
line([esp_th(5) esp_th(5)],[0 1],'Color','k');
plot([75 75],[0 1],'g--');
hold on;
plot([-75 -75],[0 1],'g--');
hold on;
plot([60 60],[0 1],'g--');
hold on;
plot([-60 -60],[0 1],'g--');
hold on;
plot([0 0],[0 1],'g--');
hold on;
legend('MF','MVDR','MUSIC','ESPRIT','true');
title('P=[-7 -5 -7 -5 -7] R known');
xlabel('angle');
ylabel('P(angle)');
axis([-90 90 0 1]);

%% 1b
N = 50;

A = [a_theta1 a_theta2 a_theta3 a_theta4 a_theta5];
Po = diag(sqrt(P_watt));
B = sign(randn(5,N));
W = (randn(m,N)+1j*randn(m,N))/sqrt(2);

Y = A*Po*B+W;
R = (1/N).*(Y*Y');

% Singular value decomposition
[U,S,V] = svd(R);
Qn = U(:,6:20);

angles = [-90:1:90];

for i=1:length(angles)
	
	% Compute a(theta)
	for ii=0:(m-1)
		a_theta(ii+1,1) = exp(-1j*ii*pi*sin(deg2rad(angles(i))));
	end

	P_mf(i) = real(a_theta'*R*a_theta);
	P_mvdr(i) = real(1/(a_theta'*inv(R)*a_theta));
	P_music(i) = 1/(norm(Qn'*a_theta).^2);

end

P_mf = P_mf./max(P_mf);
P_mvdr = P_mvdr./max(P_mvdr);
P_music = P_music./max(P_music);

% ESPIRIT method
var = diag(S);
var = mean(var(6:m));
RR = (R-var*eye(m));
R_esp = [RR(1:m-1,:);RR(2:m,:)];

[Ue,Se,Ve] = svd(R_esp);
U1 = Ue(1:38/2,1:5);
U2 = Ue(38/2+1:38,1:5);

R1 = U1'*U1;
R2 = U1'*U2;

[P,D] = eig(inv(R1)*R2);

w = diag(D);
omega = -1j*log(((w)));
esp_th = real(asind(-omega/pi));


figure();
plot(angles,P_mf,'b');
hold on;
plot(angles,P_mvdr,'r');
hold on;
plot(angles,P_music,'m');
hold on;
line([esp_th(1) esp_th(1)],[0 1],'Color','k');
line([esp_th(2) esp_th(2)],[0 1],'Color','k');
line([esp_th(3) esp_th(3)],[0 1],'Color','k');
line([esp_th(4) esp_th(4)],[0 1],'Color','k');
line([esp_th(5) esp_th(5)],[0 1],'Color','k');
plot([75 75],[0 1],'g--');
hold on;
plot([-75 -75],[0 1],'g--');
hold on;
plot([60 60],[0 1],'g--');
hold on;
plot([-60 -60],[0 1],'g--');
hold on;
plot([0 0],[0 1],'g--');
hold on;
legend('MF','MVDR','MUSIC','ESPRIT','true');
title('P=[-7 -5 -7 -5 -7] N=50');
xlabel('angle');
ylabel('P(angle)');
axis([-90 90 0 1]);

%% 1c
N = 200;

A = [a_theta1 a_theta2 a_theta3 a_theta4 a_theta5];
Po = diag(sqrt(P_watt));
B = sign(randn(5,N));
W = (randn(m,N)+1j*randn(m,N))/sqrt(2);

Y = A*Po*B+W;
R = (1/N).*(Y*Y');

% Singular value decomposition
[U,S,V] = svd(R);
Qn = U(:,6:20);

angles = [-90:1:90];

for i=1:length(angles)
	
	% Compute a(theta)
	for ii=0:(m-1)
		a_theta(ii+1,1) = exp(-1j*ii*pi*sin(deg2rad(angles(i))));
	end

	P_mf(i) = real(a_theta'*R*a_theta);
	P_mvdr(i) = real(1/(a_theta'*inv(R)*a_theta));
	P_music(i) = 1/(norm(Qn'*a_theta).^2);

end

P_mf = P_mf./max(P_mf);
P_mvdr = P_mvdr./max(P_mvdr);
P_music = P_music./max(P_music);

% ESPIRIT method
var = diag(S);
var = mean(var(6:m));
RR = (R-var*eye(m));
R_esp = [RR(1:m-1,:);RR(2:m,:)];

[Ue,Se,Ve] = svd(R_esp);
U1 = Ue(1:38/2,1:5);
U2 = Ue(38/2+1:38,1:5);

R1 = U1'*U1;
R2 = U1'*U2;

[P,D] = eig(inv(R1)*R2);

w = diag(D);
omega = -1j*log(((w)));
esp_th = real(asind(-omega/pi));


figure();
plot(angles,P_mf,'b');
hold on;
plot(angles,P_mvdr,'r');
hold on;
plot(angles,P_music,'m');
hold on;
line([esp_th(1) esp_th(1)],[0 1],'Color','k');
line([esp_th(2) esp_th(2)],[0 1],'Color','k');
line([esp_th(3) esp_th(3)],[0 1],'Color','k');
line([esp_th(4) esp_th(4)],[0 1],'Color','k');
line([esp_th(5) esp_th(5)],[0 1],'Color','k');
plot([75 75],[0 1],'g--');
hold on;
plot([-75 -75],[0 1],'g--');
hold on;
plot([60 60],[0 1],'g--');
hold on;
plot([-60 -60],[0 1],'g--');
hold on;
plot([0 0],[0 1],'g--');
hold on;
legend('MF','MVDR','MUSIC','ESPRIT','true');
title('P=[-7 -5 -7 -5 -7] N=200');
xlabel('angle');
ylabel('P(angle)');
axis([-90 90 0 1]);

%% 1d
N = 1000;

A = [a_theta1 a_theta2 a_theta3 a_theta4 a_theta5];
Po = diag(sqrt(P_watt));
B = sign(randn(5,N));
W = (randn(m,N)+1j*randn(m,N))/sqrt(2);

Y = A*Po*B+W;
R = (1/N).*(Y*Y');

% Singular value decomposition
[U,S,V] = svd(R);
Qn = U(:,6:20);

angles = [-90:1:90];

for i=1:length(angles)
	
	% Compute a(theta)
	for ii=0:(m-1)
		a_theta(ii+1,1) = exp(-1j*ii*pi*sin(deg2rad(angles(i))));
	end

	P_mf(i) = real(a_theta'*R*a_theta);
	P_mvdr(i) = real(1/(a_theta'*inv(R)*a_theta));
	P_music(i) = 1/(norm(Qn'*a_theta).^2);

end

P_mf = P_mf./max(P_mf);
P_mvdr = P_mvdr./max(P_mvdr);
P_music = P_music./max(P_music);

% ESPIRIT method
var = diag(S);
var = mean(var(6:m));
RR = (R-var*eye(m));
R_esp = [RR(1:m-1,:);RR(2:m,:)];

[Ue,Se,Ve] = svd(R_esp);
U1 = Ue(1:38/2,1:5);
U2 = Ue(38/2+1:38,1:5);

R1 = U1'*U1;
R2 = U1'*U2;

[P,D] = eig(inv(R1)*R2);

w = diag(D);
omega = -1j*log(((w)));
esp_th = real(asind(-omega/pi));


figure();
plot(angles,P_mf,'b');
hold on;
plot(angles,P_mvdr,'r');
hold on;
plot(angles,P_music,'m');
hold on;
line([esp_th(1) esp_th(1)],[0 1],'Color','k');
line([esp_th(2) esp_th(2)],[0 1],'Color','k');
line([esp_th(3) esp_th(3)],[0 1],'Color','k');
line([esp_th(4) esp_th(4)],[0 1],'Color','k');
line([esp_th(5) esp_th(5)],[0 1],'Color','k');
plot([75 75],[0 1],'g--');
hold on;
plot([-75 -75],[0 1],'g--');
hold on;
plot([60 60],[0 1],'g--');
hold on;
plot([-60 -60],[0 1],'g--');
hold on;
plot([0 0],[0 1],'g--');
hold on;
legend('MF','MVDR','MUSIC','ESPRIT','true');
title('P=[-7 -5 -7 -5 -7] N=1000');
xlabel('angle');
ylabel('P(angle)');
axis([-90 90 0 1]);