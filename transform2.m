function T_1= transform2(T_0, x)

TR_Thr= [0.8324 -0.3398 0.4378;
    0.4511 -0.0434 -0.8914;
    0.3219 0.9395 0.1172;];

TR_Skew= [1 0 .364;
    0 1 0;
    0 0 1.064;];

% R3=[1/sin(thr3) 0 0;
%     0 1 0;
%     -1/tan(thr3) 0 1;];
% R3inv=[1/sin(thr3) 0 0;
%     0 1 0;
%     -1/tan(thr3) 0 1;];
if x==1
    
T_1= TR_Skew*TR_Thr*T_0;
elseif x==0
        T_1= TR_Thr^(-1)*TR_Skew^(-1)*T_0;
% elseif x==2
%     T_1= R2*R1*T_0*(R2*R1)';
elseif x==2 
    T_1=TR_Thr*T_0; 
    
    
else
    T_1=NaN;
end

end