function [dydt] = hhdiff_reduced(t,y)

%must be in same condition as q1p1
%Uncomment to simulate without threshold & rebound behavior observance
% is = 5;
% id = 0.35;
% im = 25;
% 
% I=zeros(1,length(t));
% for t0=1:length(t)
%     if t(t0)>is && t(t0)<(is+id)
%         I(t0)=im;
%     end
% end

% Problem1 - (iv) Threshold & Rebound behavior
is = 2;
id = 0.5;
I=zeros(1,length(t));
for t0=1:length(t)
    if t(t0)>is && t(t0)<(is+id)
        I(t0)=1;
    elseif t(t0)>7 && t(t0)<(7+id)
        I(t0)=10;
    elseif t(t0)>25 && t(t0)<(25+id)
        I(t0)=25;
    elseif t(t0)>40 && t(t0)<(40+id)
        I(t0)=25;
    elseif t(t0)>60 && t(t0)<(60+id)
        I(t0)=25;
    elseif t(t0)>70 && t(t0)<(70+id)
        I(t0)=-150;
    end
end


V=y(1,1);
m=y(2,1);
n=y(3,1);
h=y(4,1);
vRest=-65;

vMem=V-vRest; %1st vMem is -65 - (-65) = 0
gNa = 120; %mS/cm^2
gK = 36; %mS/cm^2
gL = 0.3; %ms/cm^2

mC = 1; %uF/cm^2
lNp = -61; %mV
nNa = 55;
nK = -77;


a_n=.01*(10-vMem)./(exp((10-vMem)/10)-1);
b_n=.125*exp(-vMem/80);

a_m=.1*(25-vMem)./(exp((25-vMem)/10)-1);
b_m=4*exp(-vMem/18);

a_h=0.07*exp(-vMem/20);
b_h=1/(exp((30-vMem)/10)+1);

m_inf = a_m / (a_m+b_m);

%ode V, n, m, h
dydt(1,1)= I - (gNa * m_inf^3 * vMem * (0.89-(1.1.*n)) * (V - nNa) ) - (gK*n^4*(V-nK)) - (gL*(V-lNp)) ; 
% dydt(1,1)= I- (gNa*m^3*h*(V-nNa)) - (gK*n^4*(V-nK)) - (gL*(V-lNp)) ; 
dydt(2,1)= (a_m*(1-m)-b_m*m);
dydt(3,1)= (a_n*(1-n)-b_n*n); 
dydt(4,1)= (a_h*(1-h)-b_h*h);

end
