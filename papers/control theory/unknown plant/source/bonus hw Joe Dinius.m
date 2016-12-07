close all; clear all; clc
%bonus hw code
%Joe Dinius
%ame 455
load data;

freq = w;
modG = a;
phaG = b*pi/180;

x = modG.*cos(phaG);
y = modG.*sin(phaG);

figure(1);
plot(x,y,'b',x,-y,'b')
axis([-10 10 -10 10]);
grid on;
xlabel('Re')
ylabel('Im')
title('Nyquist Plot of Sampled Data')

figure(2);
modGdB = 20*log10(a);
phaGdeg= b;
subplot(2,1,1)
grid on;
semilogx(freq,modGdB)
ylabel('M_{dB}');
title('Bode Plot of Sampled Data');
subplot(2,1,2);
grid on;
semilogx(freq,phaGdeg);
ylabel('\phi (deg)');
xlabel('\omega (rad/sec)');

%generate guess
o1 = 5.51;
o2 = 2.07172;
x1 = -1.5266;
y2 = 9.21991;

A1  = [1 -o2^2 0;
       1 -o1^2 1/x1];
rfA1= rref(A1);
a0  = rfA1(1,3);
a2  = rfA1(2,3);

A2  = [1 -o1^2 0;
       1 -o2^2 1/(o2*y2)];
rfA2= rref(A2);
a1  = rfA2(1,3);
a3  = rfA2(2,3);

num = 1; 
den = [a3 a2 a1 a0];

TF = tf(num,den);

[mag,phase,wG1] = bode(TF);
magG(1:42) = mag(1,1,:);
phaseG(1:42) = phase(1,1,:);

[reG,imG,wG] = nyquist(TF);
reG1(1:135) = reG(1,1,:);
imG1(1:135) = imG(1,1,:);

%compare guess to model
figure(3)
plot(x,y,'b',x,-y,'b',reG1,imG1,'g*',reG1,-imG1,'g*')
title('Nyquist Comparison')
xlabel('Re')
ylabel('Im')
legend('Act','Act','Model','Model');

figure(4)
modGdB = 20*log10(a);
phaGdeg= b;
subplot(2,1,1)
grid on;
semilogx(freq,modGdB,wG1,20*log10(magG),'g*')
ylabel('M_{dB}');
title('Bode Comparison');
subplot(2,1,2);
grid on;
semilogx(freq,phaGdeg,wG1,phaseG,'g*');
ylabel('\phi (deg)');
xlabel('\omega (rad/sec)');

%generate open and closed loop responses
figure(5)
step(TF)
figure(6)
step(feedback(TF,1));

%generate root locus
figure(7)
rlocus(TF)

%generate final plot
figure(8)
step(feedback(TF,.05))