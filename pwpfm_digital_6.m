
clear
clc
t_sample=1e-3; % time step for the pwpfm
Isp=60; % Specific Impulse
g=9.8; 

%% Initial and Desired Conditions
r=[1 0 0];
r=r/norm(r); 
theta=0;
q=[cosd(theta/2),r(1)*sind(theta/2),r(2)*sind(theta/2),r(3)*sind(theta/2)];
omg0= [0;0;0]; % Starting from rest. 
q0_nm= colmnmtrx(q); % Initial Quaternion
r=[1 0 0 ];
r=r/norm(r); 
 theta=5; % Rotation Angle
q=[cosd(theta/2),r(1)*sind(theta/2),r(2)*sind(theta/2),r(3)*sind(theta/2)];
qd_nm= colmnmtrx(q); % Final Quaternion
%S/C MOI matrix in Global BC
global I
I= [.5518, 0, 0;
    0, .5796, 0;
    0, 0, .6156];
%% Control Variables and parameters

qi=q0_nm;
omgi=omg0; %PI frame
q_doti= 0.5*Xi(qi)*omgi; % Quaternion Kinematics
%qed= [0.01, 0.01, 0.01, 0.01]'; % For termination by error
d_theta=abs((Xi(qd_nm))'*qi*2);
[kp(1),kd(1)]=controller(5,1);
[kp(2),kd(2)]=controller(5,2);
[kp(3),kd(3)]=controller(5,3);
kp=[kp(1) 0 0;
    0 kp(2) 0;
    0 0 kp(3)];
kd=[kd(1) 0 0;
    0 kd(2) 0;
    0 0 kd(3)]; % Controller values calculated according to saturation limits. 

delt= 1; % time step for the controller algorithm
%% 
%PWPFM
global Km Tm
Km=4.5;
min_on=1e-3*3; % Valve limits 

Uon=0.1;
Cdz= Uon/Km;
Uoff=0.05;
Csat= 1+ Uoff/Km;
h=Uon-Uoff;
%Tm=-(min_on)/log(1-h/Km);
Tm=0.3;
U= .02;
No=100; % Number of Iterations
quat = q0_nm;
u=[0;0;0]; % starting from a trigger-off state
f0=[0;0;0];
f=f0;
T=[t_sample;t_sample;t_sample];
trigger=[0;0;0];
trig=0;
 tr=[0;0;0];
 m_prop=0;
 peak=[0;0;0];
 DA=eul(qd_nm);
%% Control Loop
i=1;
for it=1:No
    % Lyapunov stabilised torque function
    dltq_vm=(Xi(qd_nm))'*qi; 
    d_theta=eul(qi);
    dqq(it)=dltq_vm(1);
    dq(:,it)=dltq_vm;
    dlq4= qi'*qd_nm;
    M_p = -kp*sign(dlq4)*dltq_vm*2- kd*omgi; % Control torque, in Inertial Frame
    M_tr=transform2(M_p,1); % Transforming into orthogonal Thruster frame
    thr(it,1) = M_tr(1)/(2*0.260); % Transforming into skewed  Thruster Frame
    thr(it,2) = M_tr(2)/(2*0.260);
    thr(it,3) = M_tr(3)/(0.260*2);
    Tc = (M_tr./abs(M_tr))*U; %direction determiner
    C(it,:)=abs(thr(it,:))./U; % The input to the PWPFM
    thrust=u.*Tc;
     for j=1:3
        if C(it,j)>=Csat
             C(it,j)=Csat;
          end 
       if C(it,j)<=Cdz
              C(it,j)=Cdz;
         end 
     end 
    %% PWPFM SIMULATION
    t0=0;
    
  th=[0;0;0];
    for time=t_sample:t_sample:delt
        for j=1:3
            f(j)=filter_output(t_sample,u(j),f(j),C(it,j));
            if f0(j)==0 || f0(j)==Uoff
                if f(j)>=Uon
                    f(j)=Uon;
                    f0(j)=Uon;
                    u(j)=1;
                    trigger(j)=1;
                     tr(j)=tr(j)+1;
                     tim(j,tr(j))=T(j);
                    T(j)=t_sample;
                    trig=1; 
                else 
                    T(j)=T(j)+t_sample;
                end
            else
                if f(j)<=Uoff
                    f(j)=Uoff;
                    f0(j)=Uoff;
                    u(j)=0;
                    trigger(j)=1;
                     tr(j)=tr(j)+1;
                     tim(j,tr(j))=T(j);
                    T(j)=t_sample; 
                    trig=1; 
                else
                    T(j)=T(j)+t_sample;
                end
            end
        end
       
                
             M_tr=[thrust(1);thrust(2);thrust(3)]*0.260*2;
           M_p=transform2(M_tr,0); 
           w_dotr=inv(I)*(M_p-cross(omgi,I*omgi)); % Euler equations
           q_doti= 0.5*Xi(qi)*omgi;
           omgi=omgi+w_dotr*t_sample; 
           % Updating the quaternions
           qi(1)= qi(1)+ q_doti(1)*t_sample;
           qi(2)= qi(2)+ q_doti(2)*t_sample;
           qi(3)= qi(3)+ q_doti(3)*t_sample;
           qi(4)= qi(4)+ q_doti(4)*t_sample;
           qi=qi/norm(qi);
           
          %% 
          % To update quaternion through ode solver. Currently, it is
          % linearly propagated. Comment the above updation section and
          % uncomment the below section for ode solver. 
%             t_span=[t0,time];
%              w0=[omgi(1),omgi(2),omgi(3)];
%             [t,w]=ode45(@(t,w) odefcn(t,w,M_p,I), t_span, w0);
%             omgi=[w(size(w,1),1);w(size(w,1),2);w(size(w,1),3)];
%             
%             for a=1:size(t)-1
%                 wi=w(a,:);
%                 q_doti= 0.5*Xi(qi)*wi';
%                 qi=qi+q_doti*(t(a+1)-t(a));
%                 qi=qi/norm(qi);
%             end
%             t0=time;
%             if t0==delt
%                 t0=0;
%             end 
            thrust=u.*Tc; 
        t1(i)=thrust(1);
        t2(i)=thrust(2);
        t3(i)=thrust(3);
        T1(i)=T(1);
        U1(i)=u(1);
        C1(i)=C(it,1); 
        f1(i)=f(1);
        f2(i)=f(2);
        f3(i)=f(3);
        i=i+1;
        trig=0; 
        m_prop=m_prop+(abs(thrust(1))+abs(thrust(2))+abs(thrust(3)))*t_sample/(g*Isp); % Propellant Mass calculation
        th=th+omgi*t_sample;
    end
    np(it,:)=tr; 
    if it>1
    ang_disp(it,:)=th'+ang_disp(it-1,:);
    else 
        ang_disp(it,:)=th; 
    end 
%% 
quat=[quat,qi];
    q1(it)=qi(1);
    q2(it)=qi(2);
    q3(it)=qi(3);
    q4(it)=qi(4); 
    w_rel(it,:)=omgi;
    it 
    A(it,:)=eul(qi); % Euler Angle calculation
end 

 it=1:No;
% plot(it*delt/No,q1,it*delt/No,q2,it*delt/No,q3,it*delt/No,q4);
% xlabel('time(s)');
% ylabel('quaternions');
% legend('q1','q2','q3','q4');
% figure;
% plot(it*delt,w_rel);
% xlabel('time(s)');
% ylabel('angular velocities (rad/s)')
% legend('w1','w2','w3')
% figure;
%% 
% Various Plots to use 

%  plot(it*delt,A(:,1)*180/pi)
%  xlabel('time(s)');
%  ylabel('euler angle corresponding to the X-axis (degrees)')
%  title(['Sampling rate =',num2str(t_sample)])
%  figure; 
%  clear all;
%  sampling_rate=[1e-3;5e-3;10e-3;100e-3;1];
% legend('1','2','3')
% yline(DA(1)*180/pi);
% yline(DA(2)*180/pi);
% yline(DA(3)*180/pi);
% figure;
% % plot(it*delt,C(:,1),it*delt,C(:,2),it*delt,C(:,3));
% % figure; 
% % i=1:150000;
% % subplot(1,3,1);
% % plot(i,t1);
% % subplot(1,3,2);
% % plot(i,t2);
% % subplot(1,3,3);
% % plot(i,t3);
% % figure;
% % plot(i*t_sample,f1);
% % figure;
% % plot(i*t_sample,f2);
% % figure;
% % plot(i*t_sample,f3);
% % figure;
%  plot(it*delt,ang_disp*180/pi);
%  xlabel('time(s)');
%  ylabel('angular-displacement (degrees)');
%   legend('1','2','3');
% plot(study,TP(:,1));
% xlabel('step input (Degrees)')
% ylabel('Rise Time (s)')
%            
%             