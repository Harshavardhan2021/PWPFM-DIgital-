function [kp,kd]=controller(theta,i)
global I
xi=0.7;
th=atan(sqrt(1-xi^2)/xi);
k=1;
Mx=0.0135; % Max moment possible by the thrusters about each axis
kp=Mx/(theta*pi/180);
wn=sqrt(kp*k/(I(i,i)));
kd=2*xi*wn*I(i,i);
tr=(pi-th)/(wn*sqrt(1-xi^2));
wb=wn*sqrt(1-2*xi^2+sqrt(2-4*xi^2+4*xi^4));
end  


