clear all;
close all;
clc;

%%A.1
[pulse1,t1]=srrc_pulse(1,1/10,4,0);
[pulse2,t2]=srrc_pulse(1,1/10,4,0.5);
[pulse3,t3]=srrc_pulse(1,1/10,4,1);

figure(1);
plot(t1,pulse1,t2,pulse2,t3,pulse3)
legend('a = 0','a = 0.5','a = 1');
title('(A.1) Palmoi SRRC');
xlabel('t (sec)');
ylabel('pulse(t)');
grid on;
axis([-4 4 -0.4 1.4])

%%A.2
fft1 = fftshift(fft(pulse1,1024))*1/10;
fft2 = fftshift(fft(pulse2,1024))*1/10;
fft3 = fftshift(fft(pulse3,1024))*1/10;

fft1Coe=abs(fft1.^2);
fft2Coe=abs(fft2.^2);
fft3Coe=abs(fft3.^2);

Fs=10;
F_axis = (-Fs/2) : (Fs/1024) : (Fs/2)-(Fs/1024);

figure(2)
plot(F_axis,fft1Coe,F_axis,fft2Coe,F_axis,fft3Coe)
legend('a = 0','a = 0.5','a = 1');
title('(A.2) Fasmatikh Pyknothta Energeias SRRC gia N = 1024 (ME PLOT)')
xlabel('F (Hz)');
ylabel('|pulse(F)|^2')
grid on

figure(3)
semilogy(F_axis,fft1Coe,F_axis,fft2Coe,F_axis,fft3Coe);
legend('a = 0','a = 0.5','a = 1');
title('(A.2) Fasmatikh Pyknothta Energeias SRRC gia N = 1024 (ME SEMILOGY)')
xlabel('F (Hz)');
ylabel('|pulse(F)|^2')
grid on

fft1b = fftshift(fft(pulse1,2048))*1/10;
fft2b = fftshift(fft(pulse2,2048))*1/10;
fft3b = fftshift(fft(pulse3,2048))*1/10;

Fs=10;

fft1bCoe=abs(fft1b.^2);
fft2bCoe=abs(fft2b.^2);
fft3bCoe=abs(fft3b.^2);
F_axis = (-Fs/2) : (Fs/2048) : (Fs/2)-(Fs/2048);

figure(4)
plot(F_axis,fft1bCoe,F_axis,fft2bCoe,F_axis,fft3bCoe)
legend('a = 0','a = 0.5','a = 1');
title('(A.2) Fasmatikh Pyknothta Energeias SRRC gia N = 2048 (ME PLOT)')
xlabel('F (Hz)');
ylabel('|pulse(F)|^2')
grid on

figure(5)
semilogy(F_axis,fft1bCoe,F_axis,fft2bCoe,F_axis,fft3bCoe);
legend('a = 0','a = 0.5','a = 1');
title('(A.2) Fasmatikh Pyknothta Energeias SRRC gia N = 2048 (ME SEMILOGY)')
xlabel('F (Hz)');
ylabel('|pulse(F)|^2')
grid on

%%A.3
W1=(1+0)/(2*1);
W2=(1+0.5)/(2*1);
W3=(1+1)/(2*1);

figure(6)
semilogy(F_axis,fft1bCoe,F_axis,fft2bCoe,F_axis,fft3bCoe,F_axis,1/1000,F_axis,1/10^6);
legend('a = 0','a = 0.5','a = 1','c = T/1000','c = T/(10^6)');
title('(A.2) Fasmatikh Pyknothta Energeias SRRC gia N = 2048 (ME SEMILOGY)')
xlabel('F (Hz)');
ylabel('|pulse(F)|^2')
grid on

%%B.1
T = 1;
over = 10;
Ts = T/over;
 for k=0:2*4
        f1=dfilt.delay(k*T*(1/Ts));
        pulse1_kT=filter(f1,pulse1);
        prod=pulse1_kT.*pulse1;
        olok1(k+1) = sum(prod)*Ts;
         if (k==1)        
                 figure(7)
                 subplot(2,1,1)
                 plot(t1, prod);
                 legend(['k=', num2str(k) ,', a=0'])
                 title(['(B.1) Emfanisi se koino plot ton pulse(t) kai pulse(t-', num2str(k) ,'T) gia a=0'])
                 ylabel(['pulse(t)*pulse(t-', num2str(k) ,'T)'])
                 xlabel('t')
                 grid on
         end
        if (k==6)
                 subplot(2,1,2)
                 plot(t1, prod);
                 legend(['k=', num2str(k) ,', a=0'])
                 title(['(B.1) Emfanisi se koino plot ton pulse(t) kai pulse(t-', num2str(k) ,'T) gia a=0'])
                 ylabel(['pulse(t)*pulse(t-', num2str(k) ,'T)'])
                 xlabel('t')
                 grid on
         end
 end

  for k=0:2*4
        f2=dfilt.delay(k*T*(1/Ts));
        pulse2_kT=filter(f2,pulse2);
        prod=pulse2_kT.*pulse2;
        olok2(k+1) = sum(prod)*Ts;
         if (k==1)        
                 figure(8)
                 subplot(2,1,1)
                 plot(t2, prod);
                 legend(['k=', num2str(k) ,', a=0'])
                 title(['(B.1) Emfanisi se koino plot ton pulse(t) kai pulse(t-', num2str(k) ,'T) gia a=0.5'])
                 ylabel(['pulse(t)*pulse(t-', num2str(k) ,'T)'])
                 xlabel('t')
                 grid on
         end
        if (k==6)
                 subplot(2,1,2)
                 plot(t2, prod);
                 legend(['k=', num2str(k) ,', a=0'])
                 title(['(B.1) Emfanisi se koino plot ton pulse(t) kai pulse(t-', num2str(k) ,'T) gia a=0.5'])
                 ylabel(['pulse(t)*pulse(t-', num2str(k) ,'T)'])
                 xlabel('t')
                 grid on
         end
  end
 
   for k=0:2*4
        f3=dfilt.delay(k*T*(1/Ts));
        pulse3_kT=filter(f3,pulse3);
        prod=pulse3_kT.*pulse3;
        olok3(k+1) = sum(prod)*Ts;
         if (k==1)        
                 figure(9)
                 subplot(2,1,1)
                 plot(t3, prod);
                 legend(['k=', num2str(k) ,', a=0'])
                 title(['(B.1) Emfanisi se koino plot ton pulse(t) kai pulse(t-', num2str(k) ,'T) gia a=1'])
                 ylabel(['pulse(t)*pulse(t-', num2str(k) ,'T)'])
                 xlabel('t')
                 grid on
         end
        if (k==6)
                 subplot(2,1,2)
                 plot(t3, prod);
                 legend(['k=', num2str(k) ,', a=0'])
                 title(['(B.1) Emfanisi se koino plot ton pulse(t) kai pulse(t-', num2str(k) ,'T) gia a=1'])
                 ylabel(['pulse(t)*pulse(t-', num2str(k) ,'T)'])
                 xlabel('t')
                 grid on
         end
   end

%%C.1

N1=50;
N2=100;
b1=(sign(randn(N1,1))+1)/2;
b2=(sign(randn(N2,1))+1)/2;

%%C.2

X1=bits_to_2PAM(b1);
X2=bits_to_2PAM(b2);

counter=0;
counter2=0;
for i=0.1:Ts:50
    counter2=counter2+1;
    counter=counter+1;
    if (counter==10)
        Xd1(counter2)=X1(counter2/10);
        counter=0;
    else
        Xd1(counter2)=0;
        
    end
   
end

counter=0;
counter2=0;
for i=0.1:Ts:100
    counter2=counter2+1;
    counter=counter+1;
    if (counter==10)
        Xd2(counter2)=X2(counter2/10);
        counter=0;
    else
        Xd2(counter2)=0;
        
    end
   
end

tXd1 = [0:Ts:N1*T-Ts];
X=conv(Xd1,pulse3);
tX=[t3(1)+tXd1(1):Ts:t3(end)+tXd1(end)];

tXd2 = [0:Ts:N2*T-Ts];
Xa=conv(Xd2,pulse3);
tXa=[t3(1)+tXd2(1):Ts:t3(end)+tXd2(end)];

figure()
plot(tX, X);
title('(C.2.b) X(t)*pulse3(t) me 2-PAM gia N=50');
xlabel('t (sec)');
ylabel('X_d(t)');
grid on

figure()
plot(tXa, Xa);
title('(C.2.b) X(t)*pulse3(t) me 2-PAM gia N=100');
xlabel('t (sec)');
ylabel('X_d(t)');
grid on

figure()
plot(tXd1, Xd1);
title('(C.2.b) X_d(t) me 2-PAM gia N=50');
xlabel('t (sec)');
ylabel('X_d(t)');
grid on

figure()
plot(tXd2, Xd2);
title('(C.2.b) X_d(t) me 2-PAM gia N=100');
xlabel('t (sec)');
ylabel('X_d(t)');
grid on

Z=conv(X,pulse3)*Ts;
tz=[tX(1)+t3(1):Ts:tX(end)+t3(end)];
figure()
plot(tz,Z)
hold on
stem([1:N1]*T,X1)

Z=conv(Xa,pulse3)*Ts;
tz=[tXa(1)+t3(1):Ts:tXa(end)+t3(end)];
figure()
plot(tz,Z)
hold on
stem([1:N2]*T,X2)

%%4PAM
X1=bits_to_4PAM(b1);
X2=bits_to_4PAM(b2);

counter=0;
counter2=0;
for i=0.1:Ts:25
    counter2=counter2+1;
    counter=counter+1;
    if (counter==10)
        Xd1a(counter2)=X1(counter2/10);
        counter=0;
    else
        Xd1a(counter2)=0;
        
    end
   
end

counter=0;
counter2=0;
for i=0.1:Ts:50
    counter2=counter2+1;
    counter=counter+1;
    if (counter==10)
        Xd2a(counter2)=X2(counter2/10);
        counter=0;
    else
        Xd2a(counter2)=0;
        
    end
   
end

tXd1 = [0:Ts:(N1/2)*T-Ts];
X=conv(Xd1a,pulse3);
tX=[t3(1)+tXd1(1):Ts:t3(end)+tXd1(end)];

tXd2 = [0:Ts:(N2/2)*T-Ts];
Xa=conv(Xd2a,pulse3);
tXa=[t3(1)+tXd2(1):Ts:t3(end)+tXd2(end)];

figure()
plot(tX, X);
title('(C.2) X(t)*pulse3(t) me 4-PAM gia N=50');
xlabel('t (sec)');
ylabel('X_d(t)');
grid on

figure()
plot(tXa, Xa);
title('(C.2) X(t)*pulse3(t) me 4-PAM gia N=100');
xlabel('t (sec)');
ylabel('X_d(t)');
grid on

figure()
plot(tXd1, Xd1a);
title('(C.2.b) X_d(t) me 2-PAM gia N=50');
xlabel('t (sec)');
ylabel('X_d(t)');
grid on

figure()
plot(tXd2, Xd2a);
title('(C.2.b) X_d(t) me 2-PAM gia N=100');
xlabel('t (sec)');
ylabel('X_d(t)');
grid on

Z=conv(X,pulse3)*Ts;
tz=[tX(1)+t3(1):Ts:tX(end)+t3(end)];
figure()
plot(tz,Z)
hold on
stem([1:(N1/2)]*T,X1)

Z=conv(Xa,pulse3)*Ts;
tz=[tXa(1)+t3(1):Ts:tXa(end)+t3(end)];
figure()
plot(tz,Z)
hold on
stem([1:(N2/2)]*T,X2)


%-------------------------------------
CX1=bits_to_4PAM(b1);
CX2=bits_to_4PAM(b2);

counter=0;
counter2=0;
for i=0.1:Ts:25*2
    counter2=counter2+1;
    counter=counter+1;
    if (counter==10*2)
        CXd1a(counter2)=CX1(counter2/(10*2));
        counter=0;
    else
        CXd1a(counter2)=0;
        
    end
   
end

counter=0;
counter2=0;
for i=0.1:Ts:50*2
    counter2=counter2+1;
    counter=counter+1;
    if (counter==10*2)
        CXd2a(counter2)=CX2(counter2/(10*2));
        counter=0;
    else
       CXd2a(counter2)=0;
        
    end
   
end

CtXd1 = [0:Ts:(N1/2)*T*2-Ts];
CX=conv(CXd1a,pulse3);
CtX=[t3(1)+CtXd1(1):Ts:t3(end)+CtXd1(end)];

CtXd2 = [0:Ts:(N2/2)*T*2-Ts];
CXa=conv(CXd2a,pulse3);
CtXa=[t3(1)+CtXd2(1):Ts:t3(end)+CtXd2(end)];


figure()
plot(CtXd1, CXd1a);
title('(C.2.b) X_d(t) me 2-PAM gia N=50');
xlabel('t (sec)');
ylabel('X_d(t)');
grid on

figure()
plot(CtXd2, CXd2a);
title('(C.2.b) X_d(t) me 2-PAM gia N=100');
xlabel('t (sec)');
ylabel('X_d(t)');
grid on

CZ=conv(CX,pulse3)*Ts;
Ctz=[CtX(1)+t3(1):Ts:CtX(end)+t3(end)];
figure()
plot(Ctz,CZ)
hold on
stem([1:(N1/2)]*T*2,X1)

CZ=conv(CXa,pulse3)*Ts;
Ctz=[CtXa(1)+t3(1):Ts:CtXa(end)+t3(end)];
figure()
plot(Ctz,CZ)
hold on
stem([1:(N2/2)]*T*2,X2)