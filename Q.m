function V=Q(V0,qi,x)
q0=qi(4);
q1=qi(1);q2=qi(2);q3=qi(3);
R=[q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3), 2*(q0*q2+q1*q3);
    2*(q0*q3+q1*q2), q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1);
    2*(q1*q3-q0*q2), 2*(q0*q1+q2*q3), q0^2-q1^2-q2^2+q3^2];
if x==0
V=R*V0;
end
if x==1
    V=R^(-1)*V0;
end 
if x==2
    V=R^(-1)*V0;
    V=V(1);
end
if x==3
    V=R^(-1)*V0;
    V=V(2);
end
if x==4
    V=R^(-1)*V0;
    V=V(3);
end
if x==5
    V=inv(R)*V0*(inv(R))'; 
end 

end 