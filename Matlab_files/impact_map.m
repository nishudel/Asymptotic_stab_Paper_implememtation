function y=impact_map(x)
qsd = [x(4);x(5);x(6)];
%findin new values of velocity and the impulse
%Formula 3.23 and 3.24 on the book

m=5;
l=0.5;
r=1.0;
Mt=10;
Mh=15;
g=9.81;

De(1,1)= (r^2*(4*Mh + 4*Mt + 5*m))/4;
De(1,2)=-(m*r^2*cos(x(1) - x(2)))/2;
De(1,3)=  Mt*l*r*cos(x(1)-x(3));
De(1,4)= (r*cos(x(1))*(2*Mh + 2*Mt + 3*m))/2;
De(1,5)= -(r*sin(x(1))*(2*Mh + 2*Mt + 3*m))/2;
De(2,1)=De(1,2);
De(2,2)= (m*r^2)/4;
DE(2,3)=0;
De(2,4)=-(m*r*cos(x(2)))/2;
De(2,5)= (m*r*sin(x(2)))/2;
De(3,1)=De(1,3);
De(3,2)=De(2,3);
De(3,3)= Mt*l^2;
De(3,4)=Mt*l*cos(x(3));
De(3,5)=-Mt*l*sin(x(3));
De(4,1)=De(1,4);
De(4,2)=De(2,4);
De(4,3)=De(3,4);
De(4,4)= Mh + Mt + 2*m;
De(4,5)=0;
De(5,1)=De(1,5);
De(5,2)=De(2,5);
De(5,3)=De(3,5);
De(5,4)=De(4,5);
De(5,5)= Mh + Mt + 2*m;


E2 =[  cos(x(1)), -cos(x(2)), 0, 1, 0;
 -sin(x(1)),  sin(x(2)), 0, 0, 1];

A = [De -E2'; E2 zeros(2,2)];
B = [De*[qsd;zeros(2,1)]; zeros(2,1)];
q_plus = A\B;
y = [x(2); x(1); x(3); q_plus(2); q_plus(1); q_plus(3)];