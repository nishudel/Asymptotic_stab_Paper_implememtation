function Fc= closed_loop(t,q)

m=5;
l=0.5;
r=1.0;
Mt=10;
Mh=15;
g=9.81;
qd=[q(4);q(5);q(6)];

  
D(1,1)=(r^2*(4*Mh + 4*Mt + 5*m))/4;
D(1,2)=-(m*r^2*cos(q(1) - q(2)))/2;
D(1,3)=Mt*l*r*cos(q(1) - q(3));
D(2,1)=D(1,2);
D(2,2)= (m*r^2)/4;
D(2,3)=0;
D(3,1)=D(1,3);
D(3,2)=D(2,3);
D(3,3)= Mt*l^2;

G =[ -(g*r*sin(q(1))*(2*Mh + 2*Mt + 3*m))/2; (g*m*r*sin(q(2)))/2; -Mt*g*l*sin(q(3))];

C(1,1)=0;
C(1,2)= -(m*q(5)*r^2*sin(q(1) - q(2)))/2;
C(1,3)=Mt*l*q(6)*r*sin(q(1) - q(3));
C(2,1)=(m*q(4)*r^2*sin(q(1) - q(2)))/2;
C(2,2)=0;
C(2,3)=0;
C(3,1)=-Mt*l*q(4)*r*sin(q(1) - q(3));
C(3,2)=0;
C(3,3)=0;

B =[-1,0; 0,-1; 1,1];

Fx= [qd; D\(-C*qd - G)];


Gx=[zeros(3,2);D\B];

Jac_LFH=[0 0 0 0 0 1;0 0 0 1 1 0];

%controller

LF2h= Jac_LFH*Fx;
LGLFH=Jac_LFH*Gx;

Y=[q(3)-pi/6;q(1)+q(2)];
Yd=[q(6);q(4)+q(5)];

%gains
Kp=64*eye(2);
Kd=16*eye(2);
V=-Kp*Y-Kd*Yd;

U = (LGLFH)\(-LF2h+V);

 Fc=Fx+Gx*U; 
 
