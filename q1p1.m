clc;
clear;
close all;
% MAIN PROGRAM P1, P2

%Initial parameters
y0(1,1) = -65; % Initial Membrane Voltage
y0(2,1) = 0.05; %m - sodium gate
y0(3,1) = 0.32; %n - pottasium
y0(4,1) = 0.59; %h - leak

gNa = 120; %mS/cm^2
gK = 36; %mS/cm^2
gL = 0.3; %ms/cm^2
mC = 1; %uF/cm^2

%Na, K Nernst Potential
nL = -61; %mV
nNa = 55;
nK = -77;

dt=[0,10];  % time of integration in ms. Change according to each simulations

options=odeset('RelTol',1e-4,'AbsTol',[1e-8, 1e-8, 1e-8, 1e-8],'MaxStep',0.01);
[t,y]=ode45('hh_diff_eq',dt,y0,options);

iK = ((y(:,3).^4) .* (y(:,1)-nK)) .* gK ;
iNa = ((y(:,2).^3) .* y(:,4) .* (y(:,1)-nNa)) .* gNa ;
iL = (y(:,1)-nL) .* gL ;

%====Preparation Done. Comment/Uncomment to view Results from here=======


%Uncomment to simulate without threshold & rebound behavior observance
%setup, must be in same condition as hh_diff_eq
is = 5; %square wave(current) input start time
id = 0.35;
im = 50; %current input magnitude
I=zeros(1,length(t));
for t0=1:length(t)
    if t(t0)>is && t(t0)<(is+id)
        I(t0)=im;
    end
end

%Uncomment to test Threshold & Rebound behavior
% is = 2;
% id = 0.5;
% I=zeros(1,length(t));
% for t0=1:length(t)
%     if t(t0)>is && t(t0)<(is+id)
%         I(t0)=1;
%     elseif t(t0)>7 && t(t0)<(7+id)
%         I(t0)=10;
%     elseif t(t0)>25 && t(t0)<(25+id)
%         I(t0)=25;
%     elseif t(t0)>40 && t(t0)<(40+id)
%         I(t0)=25;
%     elseif t(t0)>60 && t(t0)<(60+id)
%         I(t0)=25;
%     elseif t(t0)>70 && t(t0)<(70+id)
%         I(t0)=-150;
%     end
% end

V=y(:,1);
m=y(:,2);
n=y(:,3);
h=y(:,4);

iNa = (m.^3 .* h .* (V-nNa)) .* gNa ;
iK = (n.^4 .* (V-nK)) .* gK;
iL = (V-nL) .* gL;
iC = V;

%Uncomment below for all plots, HAS TO BE SIMULATED SEPEARATELY 
%For Threshold & Rebound behavior
%Action Potential Evolution
% figure(1);
% plot(t, I-98, t,y(:,1));
% title('Square pulse triggering action potential')
% xlabel('time (ms)')
% ylabel('Voltage (mV)')
% legend('Current Input', 'Membrane Pot. (mV)')

% %Conductance Evolution 
% figure(2);
% p1 = plot(t, gK .* n.^4 );
% hold on
% p2 = plot(t,gNa .* (m.^3) .* h);
% legend([p1, p2], 'Potassium Conductance', 'Sodium Conductance')
% ylabel('Conductance, mS/cm^2')
% xlabel('time (ms)')
% title('Conductance Evolution for Gna and GK')
% 
% %m,n,h evolution
% figure(3);
% plot(t,y(:,2)) %m
% hold all
% plot(t,y(:,3)) %n
% plot(t,y(:,4)) %h
% title('Change in m,n,h gates')
% legend('m (Na+ Channel Opening)','n (K+ Channel Opening)','h')
% xlabel('time (ms)')
% ylabel('HH Variables')
% 
%Current type evolution
% figure(4);
% plot(t,iNa,t,iK,t, iL, t, 65 + iC)
% ylabel('Current Density (J)')
% xlabel('time (ms)')
% title('Current Evolution')
% legend('iNa','iK', 'iL', 'iC')

% Threhold and Rebound behavior
% figure(5);
% plot(t, I-98, t,y(:,1));
% title('Threhold and Rebound behavior')
% xlabel('time (ms)')
% ylabel('Voltage (mV)')
% legend('Current Input', 'Membrane Pot. (mV)', 'Location','southwest')

 