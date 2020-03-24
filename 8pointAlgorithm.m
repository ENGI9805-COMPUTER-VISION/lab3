clear all
close all
L=300;   %size of the image in pixels
I=zeros(L,L);
 
%define f, cx, cy
f=L;
cx=L/2;
cy=L/2;
 
%create intrinsic camera parameters matrix
Kint=[f 0 cx; 0 f cy; 0 0 1];
 
%create some points on the face of the cube. Homogeneous coord.
P_M=[ 0 0 0 0 0 0 0 0 0 1 2 1 2 1 2;
      2 1 0 2 1 0 2 1 0 0 0 0 0 0 0;
      0 0 0 -1 -1 -1 -2 -2 -2 0 0 -1 -1 -2 -2;
      1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
  NPTS=length(P_M);
  
  %Define pose of cube model with respect to cam1

  ax= 120*pi/180;
  ay=0;
  az=60*pi/180;
  
 Rx=[1 0 0; 0 cos(ax) -sin(ax); 0 sin(ax) cos(ax)];
 Ry=[cos(ay) 0 sin(ay); 0 1 0; -sin(ay) 0 cos(ay)];
 Rz=[cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
 
 R_m_c1=Rx*Ry*Rz;
 t_m_c1=[0;0;5];
 
 M1=[R_m_c1 t_m_c1];

%render image 1
 p1=M1* P_M;
 
 p1(1,:)=p1(1,:) ./ p1(3,:);
 p1(2,:)=p1(2,:) ./ p1(3,:);
 p1(3,:)=p1(3,:) ./ p1(3,:);
 
 figure(1), imshow(I,[]), title('View 1');

%convert image points from normalized to unnormalized
 p1u=Kint * p1;
 for i=1:length(p1u)
     rectangle('Position', [p1u(1,i)-2 p1u(2,i)-2 4 4], 'FaceColor', 'r');
 end

pause;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %set second view 
 %Define rotation of camera2 with respect to camera1
  ax= 0; ay=-25*pi/180; az=0;  
 Rx=[1 0 0; 0 cos(ax) -sin(ax); 0 sin(ax) cos(ax)];
 Ry=[cos(ay) 0 sin(ay); 0 1 0; -sin(ay) 0 cos(ay)];
 Rz=[cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
 
 R_c2_c1=Rx*Ry*Rz;
 %translation of camera2 wrt camera1. cam2 is origin
 t_c2_c1=[3;0;1];
 %calculate the pose of model wrt camera2
 H_m_c1=[R_m_c1 t_m_c1; 0 0 0 1];
 H_c2_c1=[R_c2_c1 t_c2_c1 ; 0 0 0 1];
 H_c1_c2= inv(H_c2_c1);
 H_m_c2=H_c1_c2*H_m_c1;
 
 R_m_c2=H_m_c2(1:3,1:3);
 t_m_c2=H_m_c2(1:3,4);
 %externsic camera 2 matrix
 M2=[R_m_c2 t_m_c2];

 %%%%%%%%%%%%%%%%%%%
 %render image 2
 
 p2=M2*P_M;
  
 p2(1,:)=p2(1,:) ./ p2(3,:);
 p2(2,:)=p2(2,:) ./ p2(3,:);
 p2(3,:)=p2(3,:) ./ p2(3,:);
 
 figure(2), imshow(I,[]), title('View 2');
 
 %convert image points from normalized to unnormalized
 p2u=Kint * p2;
 for i=1:length(p2u)
     rectangle('Position', [p2u(1,i)-2 p2u(2,i)-2 4 4], 'FaceColor', 'b');
 end
 pause;


%calculate the true essential matrix
 
%E=[t]x * R 
 
 t=t_c2_c1;
 E_true=[ 0 -t(3) t(2); t(3) 0 -t(1); -t(2) t(1) 0] * R_c2_c1;
 disp('true E = ');
 disp(E_true);
 pause;

 % Apply preconditioning
 xn=p1(1:2, :);  %xn is a 2xN matrix  
 N=size(xn,2);
 t=(1/N)* sum(xn,2);   %(x,y) centroid of the points
 xnc=xn-t*ones(1,N);  % center the points
%dc is 1xN vector = distance of each new position to (0,0) 
dc=sqrt(sum(xnc.^2)); 
 davg=(1/N)*sum(dc); %average distance to the origin
 s=sqrt(2)/davg; %the scale factor so that the avg dist is sqrt(2)          
 T1=[s*eye(2), -s*t; 0 0 1]; %transformation matrix 
 p1s=T1*p1;
 
 xn=p2(1:2, :);
 N=size(xn,2);
 t=(1/N)* sum(xn,2);   %(x,y) centroid of the points
 xnc=xn-t*ones(1,N);
 dc=sqrt(sum(xnc.^2));
 davg=(1/N)*sum(dc);
 s=sqrt(2)/davg;
 T2=[s*eye(2), -s*t; 0 0 1];
 p2s=T2*p2;
%compute E using 8point linear algorithm
 
 A=[p1s(1,:)'.*p2s(1,:)' p1s(1,:)'.*p2s(2,:)' p1s(1,:)' ...
     p1s(2,:)'.*p2s(1,:)' p1s(2,:)'.*p2s(2,:)' p1s(2,:)' ...
     p2s(1,:)' p2s(2,:)' ones(length(p1s),1)];
 
 [U,D,V]=svd(A);
 x=V(:,size(V,2));
 
 Escale=reshape(x,3,3)';



 %postcondition
 [U,D,V] = svd(Escale);
 %D should be[1 0 0], we enforce that by replacing 
 %it as shown in this line
 Escale=U*diag([1 1 0])*V'; 
%undo scaling done in preconditioning using
 E=T1' * Escale * T2;
 
 disp('Calculated E = ');
 disp(E);
 %Although the above is different but as we said it 
 %is up to a scale, hence divide by E(1,2)
 disp('calculated E after scaling= ');
 disp(E/(-E(1,2)));

