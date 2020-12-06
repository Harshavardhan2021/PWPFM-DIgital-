function dwdt=odefcn(t,w,M,I)
dwdt=zeros(3,1);
I1=I(1,1); I2=I(2,2); I3=I(3,3); 
dwdt(1)=(M(1)-w(2)*w(3)*(I3-I2))/I1;
dwdt(2)=(M(2)-w(1)*w(3)*(I1-I3))/I2;
dwdt(3)=(M(3)-w(2)*w(1)*(I2-I1))/I3;

end 

