clear all; 
close all;
I = imread('myofficePos1.jpg');
II = imread('myofficePos2.jpg');
imshow(I, [])
imshow(II, [])

%Lab2's intrinsic matrix
Fx = 636;
Fy = 636;
Cx = 317;
Cy = 240;

Mint=[Fx 0 Cx; 0 Fy Cy; 0 0 1];

X_coord=[1.5 3.5 5.5 7.5 4.5 1.5 3.5 5.5 7.5 4.5 1.5 3.5 5.5 7.5 0 0 0 0 0 0 0 0 0 0 0];
Y_coord=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.5 3.5 5.5 2.5 1.5 3.5 5.5 2.5 1.5 3.5 5.5];
Z_coord=[0.5 0.5 0.5 0.5 1.5 2.5 2.5 2.5 2.5 3.5 4.5 4.5 4.5 4.5 0.5 0.5 0.5 1.5 2.5 2.5 2.5 3.5 4.5 4.5 4.5] ;
P_M = [X_coord; Y_coord; Z_coord; ones(1,25)];

NPTS=length(P_M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define pose of model with respect to camera1
% DTR=pi/180;
ax = 1.4012;
ay = -0.5715;
az = 0.1061;

Rx=[1     0             0;
    0    cos(ax)    -sin(ax);
    0    sin(ax)    cos(ax)];
 
Ry=[cos(ay)     0             sin(ay);
    0           1               0;
    -sin(ay)    0             cos(ay)];
 
Rz=[cos(az)     -sin(az)   0;
    sin(az)     cos(az)    0;
    0             0        1];
 
R_m_c1=Rx*Ry*Rz;
Pmorg_c1=[8.4533;17.4056;59.0472];
M=[R_m_c1 Pmorg_c1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Render image 1
p1=M*P_M;
p1(1,:)=p1(1,:)./p1(3,:);
p1(2,:)=p1(2,:)./p1(3,:);
p1(3,:)=p1(3,:)./p1(3,:);
 
figure(1), imshow(I,[]), title('View 1');
%convert image points from normalized to unnormalized
u=Mint*p1;
for i=1:length(u)
    rectangle('Position',[u(1,i)-2 u(2,i)-2 4 4], 'FaceColor', 'r');
end
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup second view
%define rotation of camera 1 wrt camera 2
ax=1.2024;
ay=-1.1887;
az=0.3397;
Rx=[1     0             0;
    0    cos(ax)    -sin(ax);
    0    sin(ax)    cos(ax)];
Ry=[cos(ay)     0             sin(ay);
    0           1               0;
    -sin(ay)    0             cos(ay)];
Rz=[cos(az)     -sin(az)   0;
    sin(az)     cos(az)    0;
    0             0        1];
R_c2_c1=Rx*Ry*Rz;
Pc2org_c1=[-16.3280;16.8507;51.8341];
 
H_m_c1=[R_m_c1 Pmorg_c1 ; 0 0 0 1];
H_c2_c1=[R_c2_c1 Pc2org_c1; 0 0 0 1];
H_c1_c2= inv(H_c2_c1);
H_m_c2=H_c1_c2*H_m_c1;
R_m_c2=H_m_c2(1:3,1:3);
Pmorg_c2=H_m_c2(1:3,4);
%extrinsic parameters matrix
M=[R_m_c2 Pmorg_c2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Render image 2
p2=M*P_M;
p2(1,:)=p2(1,:)./p2(3,:);
p2(2,:)=p2(2,:)./p2(3,:);
p2(3,:)=p2(3,:)./p2(3,:);
 
figure(2), imshow(II,[]), title('View 2');
%convert image points from normalized to unnormalized
u=Mint*p2;
for i=1:length(u)
    rectangle('Position',[u(1,i)-2 u(2,i)-2 4 4], 'FaceColor', 'y');
end
pause


%calculate the essential matrix
t=Pc2org_c1;
E=[0 -t(3) t(2); t(3) 0 -t(1); -t(2) t(1) 0] * R_c2_c1;
disp(E);
 %Draw epipolcar lines
for i=1:length(p2)
    figure(2);
    rectangle('Position',[u(1,i)-2 u(2,i)-2 4 4], 'FaceColor', 'r');
    %compute el=E*p2 where el=[a,b,c] and the equation of the line
    %is ax=by=c=0
    figure(1);
    el=E*p2(:,i);
    %to find two points we assume x =1 and solve for y then assume
    % x=2 and solve for y
    px=1;
    pLine0=[px; (-el(3)-el(1)*px)/el(2); 1];
    px=-2;
    pLine1=[px; (-el(3)-el(1)*px)/el(2); 1];
   %convert to unnormalized
    pLine0=Mint*pLine0; pLine1=Mint*pLine1;
    %draw the line
    line([pLine0(1) pLine1(1)], [pLine0(2) pLine1(2)], 'Color', 'r');
    pause;
end
