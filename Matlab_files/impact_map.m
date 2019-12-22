function y=impact_map(x)


%the new angles q1 and q2 will be swapped 
q_n= [x(2);x(1);x(3)];

qsd= [x(5);x(4);x(6)];
%findin new values of velocity and the impulse
%Formula 3.23 and 3.24 on the book
jac_Ye=zeros(2,3);
Ye_abv=eye(3);
Ye=[Ye_abv;jac_Ye];

De =[               125/4, -(5*cos(x(1) - x(2)))/2, 5*cos(x(1) - x(3)), (65*cos(x(1)))/2, -(65*sin(x(1)))/2;
-(5*cos(x(1) - x(2)))/2,                 5/4,              0, -(5*cos(x(2)))/2,   (5*sin(x(2)))/2;
     5*cos(x(1) - x(3)),                   0,            5/2,      5*cos(x(3)),      -5*sin(x(3));
    (65*cos(x(1)))/2,      -(5*cos(x(2)))/2,      5*cos(x(3)),             35,               0;
  -(65*sin(x(1)))/2,       (5*sin(x(2)))/2,     -5*sin(x(3)),              0,              35];

E2 =[  cos(x(1)), -cos(x(2)), 0, 1, 0;
 -sin(x(1)),  sin(x(2)), 0, 0, 1];

delF2= -inv(E2*inv(De)*E2')*E2*Ye;

delqed= (inv(De)*E2'*delF2)+Ye;

delqsd=[delqed(2,:);delqed(1,:);delqed(3,:)]

qd_n=delqsd*qsd;

y=[q_n;qd_n];

