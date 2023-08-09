clc;
clear;
close all;

%========square wave input=============
dt = 1; 
t = 0: dt: 100;
Vrest = 0;
R = 10; %Mohm
C = 1; %nF
Vthr = 5; %mV
Vspk = 70; %mV

V(1) = Vrest;


I=zeros(1,length(t));
for t0=1:length(t)
    if t(t0)>10 && t(t0)<(60)
        I(t0)=1;
    end
end

firing = 0;

figure(1)
plot(t, I);
title('Current Stimulus');
xlabel('Time in ms');
ylabel('Current (nA)');

for i = 1: length(t)-1
    dVdt = (I(i) - V(i)/R)/C;

    if V(i) < Vthr
        V(i+1) = V(i) + dVdt*dt;
    else
        V(i) = Vspk;
        V(i+1) = 0;
        firing = firing + 1;
    end
end

figure(2)
plot(t, V);
title('Neuronal firing for square wave stimulation');
xlabel('Time in ms');
ylabel('Voltage (mV)');

%Uncomment clc, clear, close to run above simulations
% clc;
% clear;
% close all;

%========sinusoidal input with frequencies=============
dt = 1; 
t = 0: dt: 1000;
Vrest = 0;
R = 10; %Mohm
C = 1; %nF
Vthr = 5; %mV
Vspk = 70; %mV

V(1,1) = Vrest; 
V(2,1) = Vrest; 
V(3,1) = Vrest; 
V(4,1) = Vrest;
V(5,1) = Vrest;
V(6,1) = Vrest;
V(7,1) = Vrest;


f = [1,2,5,10,20,50,100];



firing = 0;
firing_freq = zeros(1,length(f));

for i = 1: length(f)
    I = sin(2*pi*f(i)*(t/1000));
    firing = 0;

    for j = 1: length(t)-1
        Isin = I(i);
        dVdt = (Isin - V(i, j)/R)/C;

        if V(i,j) < 0
            V(i,j) = 0;
        end
    
        if V(i, j) < Vthr
            V(i, j+1) = V(i, j) + dVdt*dt;
        else
            V(i, j) = Vspk;
            V(i, j+1) = 0;
            firing = firing + 1;
        end
    end
    firing_freq(i) = firing;
end

%Uncomment to see sinusoidal stimulus graph and much more
figure(1)
plot (t,30 + sin(2*pi*f(7)*(t/1000)))
hold all
plot (t,25 +sin(2*pi*f(6)*(t/1000)))
plot (t,20 +sin(2*pi*f(5)*(t/1000)))
plot (t,15 +sin(2*pi*f(4)*(t/1000)))
plot (t,10 +sin(2*pi*f(3)*(t/1000)))
plot (t,5 +sin(2*pi*f(2)*(t/1000)))
plot (t,sin(2*pi*f(1)*(t/1000)))
title('Sinusoidal Stimulus with differentiating frequency')
ylabel('Current (nA)')
xlabel('time (ms)')

figure(2)
plot (t,600 + V(7, :))
hold all
plot (t,500 + V(6, :))
plot (t,400 + V(5, :))
plot (t,300 + V(4, :))
plot (t,200 + V(3, :))
plot (t,100 + V(2, :))
plot (t, V(1, :))
title('IAF responses due to Sinusoidal Stimulation with differentiating frequency')
ylabel('Voltage (mV)')
xlabel('time (ms)')
% 
figure(3)
x = logspace(0, 2);
Y = interp1(1:length(firing_freq), firing_freq, linspace(1, length(firing_freq), length(x)), 'linear');
semilogx(x, Y)
set(gca,'YGrid','on', 'XGrid', 'on')
title('Spike Count VS Stimulus Frequency')
ylabel('Total Spikes')
xlabel('Stimulus frequency (Hz)')


%Uncomment clc, clear, close to run above simulations
clc;
clear;
close all;



% %=================Problem 2b===================
dt = 1; 
t = 0: dt: 1000;
f = [1,2,5,10,20,50,100];
a = 0.02;
c = -65;
d = 8;
b = 0.2;
v = -65;
u = b*v;
firing = 0;
firing_freq = zeros(1,length(f));
Vthr = 30; %mV
v_rec = zeros(length(f), length(t));
u_rec = zeros(length(f), length(t));
v_rec (1,1) = c; 
v_rec (2,1) = c; 
v_rec (3,1) = c; 
v_rec (4,1) = c;
v_rec (5,1) = c;
v_rec (6,1) = c;
v_rec (7,1) = c;
u_rec (1,1) = u; 
u_rec (2,1) = u; 
u_rec (3,1) = u; 
u_rec (4,1) = u; 
u_rec (5,1) = u; 
u_rec (6,1) = u; 
u_rec (7,1) = u; 


for i = 1: length(f)
    I =  10 * sin(2*pi*f(i)*(t/1000));
    firing = 0;
    j = 2;

    while j <= length(t)-1
        Isin = I(j-1);
        dv = 0.04 * v_rec(i,j-1)^2 + 5 * v_rec(i,j-1) + 140 - u_rec(i,j-1) + Isin;
        du = a * (b*v_rec(i,j-1) - u_rec(i,j-1));
        v_rec(i,j) = v_rec(i,j-1) + dt * dv;
        u_rec(i,j) = u_rec(i,j-1) + dt * du;

        if v_rec(i,j) >= Vthr
            v_rec(i,j) = Vthr;
            v_rec(i,j+1) = c;
            u_rec(i,j+1) = u_rec(i,j) + d;
            firing = firing + 1;
            j = j + 2;
        else
            j = j + 1;
        end
    end
    firing_freq(i) = firing;
end

% figure(1)
% plot (t,600 + v_rec(7, :))
% hold all
% plot (t,500 + v_rec(6, :))
% plot (t,400 + v_rec(5, :))
% plot (t,300 + v_rec(4, :))
% plot (t,200 + v_rec(3, :))
% plot (t,100 + v_rec(2, :))
% plot (t, v_rec(1, :))
% title('IAF responses due to Sinusoidal Stimulation with differentiating frequency')
% ylabel('Voltage (mV)')
% xlabel('time (ms)')
% 
% figure(2)
% x = logspace(0, 2);
% Y = interp1(1:length(firing_freq), firing_freq, linspace(1, length(firing_freq), length(x)), 'linear');
% semilogx(x, Y)
% set(gca,'YGrid','on', 'XGrid', 'on')
% title('Spike Count VS Stimulus Frequency')
% ylabel('Total Spikes')
% xlabel('Stimulus frequency (Hz)')
% % 
% figure(3)
% plot (t,600 + u_rec(7, :))
% hold all
% plot (t,500 + u_rec(6, :))
% plot (t,400 + u_rec(5, :))
% plot (t,300 + u_rec(4, :))
% plot (t,200 + u_rec(3, :))
% plot (t,100 + u_rec(2, :))
% plot (t, u_rec(1, :))
% title('neuronal responses (u) due to Sinusoidal Stimulation with differentiating frequency')
% ylabel('u - scale')
% xlabel('time (ms)')

%Uncomment clc, clear, close to run above simulations
clc;
clear;
close all;



% %=================Problem 2c===================
C = 1; %nF
R = 10; %MOhm
Vrest = 0; %mV
Vspk = 70;
Tthresh = 50; %threshold time constant ms
Einh = -15; % synpatic reversal potential
Esyn = Einh; % inhibitory reversal potential
e = 2.7183; 
Tsyn = 15; %synaptic time constant
gpeak = 0.1; %peak synaptic conductance
dt = 1; 
t = 0: dt: 1500; %total simulation time3

theta1 = zeros(1,length(t));
theta2 = zeros(1,length(t));
z1 = zeros(1,length(t));
g1 = zeros(1,length(t));
z2 = zeros(1,length(t));
g2 = zeros(1,length(t));
v1 = zeros(1,length(t));
v2 = zeros(1,length(t));


%"Both neurons should receive constant (or constantly?) current injection.
%Specifically, inject 1.1 nA into neuron 1 and 0.9 nA into neuron 2."

%Constant Current input for neuron 1
I1=zeros(1,length(t));
for t0=1:length(t)
    if t(t0)>0 && t(t0)<1501
        I1(t0)=1.1;
    end
end


%Constant Current input for neuron 2
I2=zeros(1,length(t));
for t0=1:length(t)
    if t(t0)>0 && t(t0)<1501
        I2(t0)=0.9;
    end
end

u1 = 0;
u2 = 0;

for i = 2 : length(t)-1
    % Voltage update for neuron 1 and 2
    dv1 = (-v1(i-1)/R - g1(i-1) * (v1(i-1) - Esyn) + I1(i-1) )/C;
    dv2 = (-v2(i-1)/R - g2(i-1) * (v2(i-1) - Esyn) + I2(i-1) )/C;

    if v1(i) ~= Esyn
        v1(i) = v1(i-1) + dt * dv1;
    end

    if v2(i) ~= Esyn
        v2(i) = v2(i-1) + dt * dv2;
    end

    % Threshold update for n1 and 2
    dtheta1 = (-theta1(i-1) + v1(i-1)) / Tthresh;
    dtheta2 = (-theta2(i-1) + v2(i-1)) / Tthresh;
    theta1(i) = theta1(i-1) + dtheta1*dt;
    theta2(i) = theta2(i-1) + dtheta2*dt;

    %Conductance update 1
    dz1 = (-z1(i-1) / Tsyn) + (gpeak / (Tsyn/e)) * u1;
    dz2 = (-z2(i-1) / Tsyn) + (gpeak / (Tsyn/e)) * u2;
    z1(i) = z1(i-1) + dz1*dt;
    z2(i) = z2(i-1) + dz2*dt;

    %Conductance update 2
    dg1 = (-g1(i-1) / Tsyn) + z1(i);
    dg2 = (-g2(i-1) / Tsyn) + z2(i);    
    g1(i) = g1(i-1) + dg1*dt;
    g2(i) = g2(i-1) + dg2*dt;
    
   
    if u1 == 1
        u1 = 0;
    end

    if u2 == 1
        u2 = 0;
    end

%When the neuron fires an action potential, reset the membrane voltage to Einh 
%on the next time step. 
    if v1(i) >= theta1(i)
        v1(i) = Vspk;
        v1(i+1) = Einh;
        u2 = 1; % 
    end

    if v2(i) >= theta2(i)
        v2(i) = Vspk;
        v2(i+1) = Einh;
        u1 = 1;
    end
end

figure(1)
plot(t, I1)
hold all
plot(t, 100 + I2)
plot(t, v1)
plot(t, 100 + v2)
plot(t, theta1)
plot(t, 100 + theta2)
ylim([-10 150])
title('Reciprocal Inhibition')
legend('1.1nA Constant Current','0.9 nA Constant Current', ...
    'neuron 1 AP','neuron 2 AP', 'theta1','theta2', 'Location','bestoutside')
xlabel('time (ms)')
ylabel('Scaled AP (mV)')






