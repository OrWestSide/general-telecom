function [s_errors,b_errors]=Exer3(SNR)

% clear all;
% close all;
% clc;

%% 1

%Dimiourgia 3*N tuxaiwn bits
N=50;
bits=(sign(randn(3*N,1))+1)/2;

%% 2

%Metatropi twn bits se 8-PSK
A=1;
X=bits_to_PSK_8(bits,A);
%scatterplot(X)

%% 3

%Dimiourgia palmou srrc
T=1;
over=10;
Ts=T/over;
a=0.8;

[pulse,t]=srrc_pulse(T,Ts,6,a);

%Filtrarisma tis diastasis Xi
counter=0;
counter2=0;
for i=0.1:Ts:N
    counter2=counter2+1;
    counter=counter+1;
    if (counter==10)
        X_i(counter2)=X(counter2/10,1);  %#ok<*SAGROW>
        counter=0;
    else
        X_i(counter2)=0;
    end 
end

tX_i=[0:Ts:N-Ts];

XI=conv(pulse,X_i);
tXI = [min(t)+min(tX_i):Ts:max(t)+max(tX_i)];

%Sxediasmos kumatomorfis meta tin suneliksi
figure();
plot(tXI,XI);
title('(3) - XI(t) - Kymatomorfh Eksodou') 
xlabel('t (sec)');
ylabel('XI(t)');
grid on;

%Upologismos periodogrammatos
Fs=1/Ts;
Ttotal = length(tXI)*Ts; 
FXI = fftshift(fft(XI,2048)*Ts); 
PXI = (abs(FXI).^2)/Ttotal; 
FXI_axis = (-Fs/2) : (Fs/length(FXI)) : (Fs/2)-(Fs/length(FXI));

%Sxediasmos periodogrammatos
figure();
plot(FXI_axis,PXI);
title('(3) - XI(t) - Periodogramma');
xlabel('F (Hz)');
ylabel('PXI(F)');
grid on;




%Filtrarisma tis diastasis Xq
counter=0;
counter2=0;
for i=0.1:Ts:N
    counter2=counter2+1;
    counter=counter+1;
    if (counter==10)
        X_q(counter2)=X(counter2/10,2);  %#ok<*SAGROW>
        counter=0;
    else
        X_q(counter2)=0; 
    end
end

tX_q=[0:Ts:N-Ts];

XQ=conv(pulse,X_q);
tXQ = [min(t)+min(tX_q):Ts:max(t)+max(tX_q)];

%Sxediasmos kumatomorfis meta tin suneliksi
figure();
plot(tXQ,XQ);
title('(3) - XQ(t) - Kymatomorfh Eksodou') 
xlabel('t (sec)');
ylabel('XQ(t)');
grid on;

%Upologismos periodogrammatos
Ttotal = length(tXQ)*Ts; 
FXQ = fftshift(fft(XQ,2048)*Ts); 
PXQ = (abs(FXQ).^2)/Ttotal; 
FXQ_axis = (-Fs/2) : (Fs/length(FXQ)) : (Fs/2)-(Fs/length(FXQ));

%Sxediasmos periodogrammatos
figure();
plot(FXQ_axis,PXQ);
title('(3) - XQ(t) - Periodogramma');
xlabel('F (Hz)');
ylabel('PXQ(F)');
grid on;

%% 4

F0=2;

%XI
%Pollaplasiasmos me ton katallilo forea
XI=XI.*(2*cos(2*pi*F0*tXI));

%Sxediasmos kumatomorfis meta apo ton pollaplasiasmo
figure();
plot(tXI,XI);
title('(4) - XI(t) - Pollaplasiasmeno me ton forea')
xlabel('t (sec)');
ylabel('XI(t)');
grid on;

%Upologismos periodogrammatos
Ttotal = length(tXI)*Ts;
XFI = fftshift(fft(XI,2048)*Ts);
PXI = (abs(XFI).^2)/Ttotal;
FXI_axis = (-Fs/2) : (Fs/length(XFI)) : (Fs/2)-(Fs/length(XFI));

%Sxediasmos periodogrammatos
figure(); 
plot(FXI_axis,PXI)
title('(4) - XI(t) - Periodogramma Pollaplasiasmeno me ton forea') 
xlabel('F (Hz)'); 
ylabel('PXI(F)');
grid on;


%XQ
%Pollaplasiasmos me ton katallilo forea
XQ=XQ.*((-1)*2*sin(2*pi*F0*tXQ));

%Sxediasmos kumatomorfis meta apo ton pollaplasiasmo
figure();
plot(tXQ,XQ);
title('(4) - XQ(t) - Pollaplasiasmeno me ton forea')
xlabel('t (sec)');
ylabel('XQ(t)');
grid on;

%Upologismos periodogrammatos
Ttotal = length(tXQ)*Ts;
XFQ = fftshift(fft(XQ,2048)*Ts);
PXQ = (abs(XFQ).^2)/Ttotal;
FXQ_axis = (-Fs/2) : (Fs/length(XFQ)) : (Fs/2)-(Fs/length(XFQ));

%Sxediasmos periodogrammatos
figure(); 
plot(FXQ_axis,PXQ)
title('(4) - XQ(t) - Periodogramma Pollaplasiasmeno me ton forea') 
xlabel('F (Hz)'); 
ylabel('PXQ(F)');
grid on;

%% 5

Xt=XI+XQ;
tXt=tXI;

figure();
plot(tXt,Xt);
title('(5) - X(t) - Kymatomorfh Eksodou Athroismatos');
xlabel('t (sec)');
ylabel('X(t)');
grid on;

Ttotal=length(tXt)*Ts;
XF = fftshift(fft(Xt,2048)*Ts);
PX = (abs(XF).^2)/Ttotal;
FX_axis = (-Fs/2) : (Fs/length(XF)) : (Fs/2)-(Fs/length(XF));
% 
% figure();
plot(FX_axis,PX)
title('(5) - X(t) - Periodogramma Athroismatos')
xlabel('F (Hz)'); 
ylabel('PX(F)');
grid on;

%% 6

%Idaniko kanali

%% 7

%SNR=10;
%Prosthesi leukou gaussianou thorubou stin eksodo tou kanaliou
Yt=awgn(Xt,SNR);
tYt=tXt;
%Yt=Xt;

%% 8

%XI
%Pollaplasiasmos me ton katallilo forea
YI=Yt.*(cos(2*pi*F0*tYt));
tYI=tYt;

%Sxediasmos kumatomorfis meta apo ton pollaplasiasmo
figure();
plot(tYI,YI);
title('(8) - YI(t) - Pollaplasiasmeno me ton forea kai me thorubo')
xlabel('t (sec)');
ylabel('YI(t)');
grid on;

%Upologismos periodogrammatos
Ttotal = length(tYI)*Ts;
YFI = fftshift(fft(YI,2048)*Ts);
PYI = (abs(YFI).^2)/Ttotal;
FYI_axis = (-Fs/2) : (Fs/length(YFI)) : (Fs/2)-(Fs/length(YFI));

%Sxediasmos periodogrammatos
figure(); 
plot(FYI_axis,PYI)
title('(8) - YI(t) - Periodogramma Pollaplasiasmeno me ton forea kai me thorubo') 
xlabel('F (Hz)'); 
ylabel('PYI(F)');
grid on;


%XQ
%Pollaplasiasmos me ton katallilo forea
YQ=Yt.*((-1)*sin(2*pi*F0*tYt));
tYQ=tYt;

%Sxediasmos kumatomorfis meta apo ton pollaplasiasmo
figure();
plot(tYQ,YQ);
title('(8) - YQ(t) - Pollaplasiasmeno me ton forea kai me thorubo')
xlabel('t (sec)');
ylabel('YQ(t)');
grid on;

%Upologismos periodogrammatos
Ttotal = length(tYQ)*Ts;
YFQ = fftshift(fft(YQ,2048)*Ts);
PYQ = (abs(YFQ).^2)/Ttotal;
FYQ_axis = (-Fs/2) : (Fs/length(YFQ)) : (Fs/2)-(Fs/length(YFQ));

%Sxediasmos periodogrammatos
figure(); 
plot(FYQ_axis,PYQ)
title('(8) - YQ(t) - Periodogramma Pollaplasiasmeno me ton forea kai me thorubo') 
xlabel('F (Hz)'); 
ylabel('PYQ(F)');
grid on;

%% 9

%Filtraroume tin eksodo
YIt=conv(YI,pulse)*Ts;
tYIt = [min(t)+min(tYI):Ts:max(t)+max(tYI)];

%Sxediazoume tin kumatomorfi tis eksodou
figure();
plot(tYIt,YIt);
title('(9) - YI(t) - Kymatomorfh Eksodou me thorubo');
xlabel('t (sec)');
ylabel('YI(t)');
grid on;

%Upologizoume to periodogramma
Ttotal = length(tYIt)*Ts;
YFI = fftshift(fft(YIt,2048)*Ts);
PYI = (abs(YFI).^2)/Ttotal;
FYI_axis = (-Fs/2) : (Fs/length(YFI)) : (Fs/2)-(Fs/length(YFI));

%Sxediazoume to periodogramma 
figure(); 
plot(FYI_axis,PYI)
title('(9) - YI(t) - Periodogramma')
xlabel('F (Hz)');
ylabel('PYI(F)'); 
grid on;


YQt=conv(YQ,pulse)*Ts;
tYQt = [min(t)+min(tYQ):Ts:max(t)+max(tYQ)];

%Sxediazoume tin kumatomorfi tis eksodou
figure();
plot(tYQt,YQt);
title('(9) - YQ(t) - Kymatomorfh Eksodou me thorubo');
xlabel('t (sec)');
ylabel('YQ(t)');
grid on;

%Upologizoume to periodogramma
Ttotal = length(tYQt)*Ts;
YFQ = fftshift(fft(YQt,2048)*Ts);
PYQ = (abs(YFQ).^2)/Ttotal;
FYQ_axis = (-Fs/2) : (Fs/length(YFQ)) : (Fs/2)-(Fs/length(YFQ));

%Sxediazoume to periodogramma 
figure(); 
plot(FYQ_axis,PYQ)
title('(9) - YQ(t) - Periodogramma')
xlabel('F (Hz)');
ylabel('PYQ(F)'); 
grid on;

%% 10

%Deigmatoliptoume
YI_s=YIt(130:T/Ts:620);
YQ_s=YQt(130:T/Ts:620);

Y=[YI_s',YQ_s'];

%Sxediazoume
scatterplot(Y)

%% 11

[est_X,est_bit_seq]=detect_PSK_8(Y,1);

%% 12

s_errors=symbol_errors(est_X,X);

%% 13

b_errors=bit_errors(est_bit_seq,bits);