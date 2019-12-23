clear 
clc
close all
Yout=[];
tout=[];
n=10;
ts=0;
X0(:,1) = [0.392699081699451;-0.391976695050869;0.525786543993878;1.59948933023117;-1.60439921576902;-0.0148697272114524];
% X0(:,1)=[pi/8;-pi/8;pi/3;1.58;-1.58;0];
options=odeset('Events',@impact_check);

for i=2:n
    if 50-ts<= 0.01
        break;
    else
        X0(:,i)=impact_map(X0(:,i-1));
        [t,y,te,ye,ie]=ode45(@(t,x)closed_loop(t,x),[ts:0.01:500],X0(:,i),options); 
        X0(:,i)= ye;
        ts=t(length(t));
        Yout=[Yout;y];
        tout=[tout;t];
    end
end

figure
plot(tout,Yout(:,1));
grid on
figure
plot(tout,Yout(:,2));
grid on
figure
plot(tout,Yout(:,3));
grid on
    
    
figure
plot(Yout(:,1),Yout(:,4));
grid on

figure
plot(Yout(:,2),Yout(:,5));
grid on