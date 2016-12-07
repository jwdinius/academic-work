close all; clear all; clc

omn = 5.97; %rad/sec
damp= .95;

den = [1/omn^3 3*damp^2/omn^2 3*damp/omn 1];

num = 25;

TF = tf(num,den);

figure(1);
bode(num,den);