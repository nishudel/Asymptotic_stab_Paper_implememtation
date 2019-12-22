clear
clc

% we are planning to use one coordinates 
% the generalised coordinates and their derivatives:
syms  q1 q2 q3 q1d q2d q3d real
% Masses of the legs, hip and the torso respectively assuming that they are
% lumped
syms  m Mh Mt real
% OTHER CONSTANTS gravity ,length of each leg( assumed to be equal), length
% of torso link
syms g r l real



% we are calculating the kintetic energy and potential enerdy to find the
% lagrangian equation 

q=[q1;q2;q3];
qd=[q1d;q2d;q3d];

% Position of center of mass of each link / mass in the system is
% calculated

%stance leg 
pcm_x1= (r/2)*sin(q1); 
pcm_y1= (r/2)*cos(q1);

%swing leg
pcm_x2= r*sin(q1)-(r/2)*sin(q2); 
pcm_y2= r*cos(q1)-(r/2)*cos(q2);

% hip
pcm_x3= (r)*sin(q1); 
pcm_y3= (r)*cos(q1);

%torso
pcm_x4= r*sin(q1)+l*sin(q3);
pcm_y4= r*cos(q1)+l*cos(q3);

%velocity
pcm_x1d= jacobian(pcm_x1,q)*qd;
pcm_y1d= jacobian(pcm_y1,q)*qd;


pcm_x2d= jacobian(pcm_x2,q)*qd;
pcm_y2d= jacobian(pcm_y2,q)*qd;

pcm_x3d= jacobian(pcm_x3,q)*qd;
pcm_y3d= jacobian(pcm_y3,q)*qd;


pcm_x4d= jacobian(pcm_x4,q)*qd;
pcm_y4d= jacobian(pcm_y4,q)*qd;

%kinetic energy calculation ( we are consdering a lumped mass system so the
%contribution of rotation will be zero)

ke1= (m/2)*((pcm_x1d)'*(pcm_x1d)+ (pcm_y1d)'*(pcm_y1d));
ke2= (m/2)*((pcm_x2d)'*(pcm_x2d)+ (pcm_y2d)'*(pcm_y2d));
ke3= (Mh/2)*((pcm_x3d)'*(pcm_x3d)+ (pcm_y3d)'*(pcm_y3d));
ke4= (Mt/2)*((pcm_x4d)'*(pcm_x4d)+ (pcm_y4d)'*(pcm_y4d));

%Total kinetic energy
K=ke1+ke2+ke3+ke4;

%potential energy
V=(m*pcm_y1+m*pcm_y2+Mh*pcm_y3+Mt*pcm_y4)*g;

%D,G,C matrices
D = jacobian(jacobian(K,qd),qd);
D=simplify(D)
G = jacobian(V,q);
G=simplify(G)
for k=1:3
	for j=1:3
		C(k,j)=0*g;
		for i=1:2
			C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i)) + ...
				diff(D(k,i),q(j)) - ...
				diff(D(i,j),q(k)))*qd(i);
		end
	end
end
C=simplify(C)

% Matrix B 

B=[-1 0;0 -1;1 1]
