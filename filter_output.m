function f=filter_output(time,u,f0,C)
global Km  Tm
 f=f0+(Km*(C-u)-f0)*(1-exp(-time/Tm));
end 

