
% This is to calculate the impact model of the robot
% very similar to th swing phase, we just introduce two coordinates
% indicating the position of a fixed point till the shift happens at the
% impact

clear
clc

% we are planning to use absolute coordinates :
% the generalised coordinates and their derivatives:q4,q5 are the end
% points of the feet currently in contact with the ground before impact
syms  q1 q2 q3 q4 q5 q1d q2d q3d q4d q5d real
% Masses of the legs, hip and the torso respectively assuming that they are
% lumped
syms  m Mh Mt real
% OTHER CONSTANTS gravity ,length of each leg( assumed to be equal), length
% of torso link
syms g r l real

% we are calculating the kintetic energy to find the D matrix


q=[q1;q2;q3;q4;q5];
qd=[q1d;q2d;q3d;q4d;q5d];

% Position of center of mass of each link / mass in the system is
% calculated. Here, the coordinates q4,q5 also determine the pcm or pcmd 

%stance leg 
pcm_x1= q4+(r/2)*sin(q1); 
pcm_y1= q5+(r/2)*cos(q1);

%swing leg
pcm_x2= q4+r*sin(q1)-(r/2)*sin(q2); 
pcm_y2= q5+r*cos(q1)-(r/2)*cos(q2);

% hip
pcm_x3= q4+(r)*sin(q1); 
pcm_y3= q5+(r)*cos(q1);

%torso
pcm_x4= q4+r*sin(q1)+l*sin(q3);
pcm_y4= q5+r*cos(q1)+l*cos(q3);

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

%D matrix- the only matrix that matters for the impact model
D = jacobian(jacobian(K,qd),qd);
D=simplify(D)

%Now, we calculate the E2 matrix - jacobian(p2,qe)..symbols as used in the
%book



