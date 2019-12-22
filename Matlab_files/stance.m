clear 
clc
close all
Yout=[];
tout=[];
n=5;
ts=0;
X0(:,1)=[pi/8;-pi/8;pi/6;1.58;-1.58;0];
for i=2:n
    X0(:,i)=impact_map(X0(:,i-1));
    options=odeset('Events',@impact_check);
    [t,y,te,ye,ie]=ode45(@(t,x)closed_loop(t,x),[0 50],X0(:,i),options); 
    %X0(:,i+1)=impact_map(ye);
    ts=t(length(t));
    Yout=[Yout;y];
    tout=[tout;t];
end

% figure
% plot(tout,Yout(:,1));
% grid on
% figure
% plot(tout,Yout(:,2));
% grid on
% figure
% plot(tout,Yout(:,3));
% grid on
%     
    
figure
plot(Yout(:,1),Yout(:,4));
grid on