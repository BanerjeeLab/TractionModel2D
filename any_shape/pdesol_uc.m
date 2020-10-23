%Shiladitya Banerjee[UCL]
%Edited 01.02.2017

%Compute stresses in any cell geometry
%Unit of length in microns, forces in N

clear
%load contour data
filename='cell_2_pts.txt';%choose appropriate cell contour
data=dlmread(filename);
pixelsize=9E-8;
xpts=pixelsize*data(:,2);
ypts=pixelsize*data(:,1);

tic
%define pde geometry
nm=round(numel(xpts)/2);
Ax=xpts(1:2:length(xpts))';
Ay=ypts(1:2:length(ypts))';
[Xout,Yout]=points2contour(Ax,Ay,1,'ccw');
gdm=[2 nm Xout Yout];
toc
[g1,bt]=decsg(gdm','P1',['P1']');
[g2,bt1]=csgdel(g1,bt);
toc
[p e t]=initmesh(g2);

%use the following commands to refine mesh starting with a custom mesh size
%hmax=1/16;
%[p e t]=initmesh(g1,'Hmax',hmax);
%[p e t]=refinemesh(g2,p,e,t);
%[p e t]=refinemesh(g,p,e,t);



%mechanical parameters
E=0.9E4; % Youngs modulus in Pascals
nu=0.45; % Poisson ratio
G=E/(2*(1+nu)); %shear modulus
mu=2*G*nu/(1-2*nu); %bulk modulus
ht=3E-6;% cell height
Ys=7E9;% substrate stiffness (substrate shear modulus / height) units [N/m^3]
Ya=8E8;%adhesion stiffness [N/m^3]
Yeff=1/((1/Ys)+(1/Ya));% effective stiffness of adhesion
sigma=5E3; % bulk active stress in pascals 
sigma2=0.3E2;% line tension in units [pascals X meters]
ar=815.19E-12;%area

c=[2*G+mu;0;0;G;0;mu;G;0;0;G;mu;0;G;0;0;2*G+mu];

%boundary condition
b=@(p,e,u,time) boundaryFile(p,e,ar,sigma,sigma2,ht);

a=[Yeff/ht;Yeff/ht];
f=[0;0];

%solution
u=assempde(b,p,e,t,c,a,f);

%u = pdenonlin(b,p,e,t,c,a,f,'jacobian','lumped');

%visualize
n=size(p,2);
u1=u(1:n,1);
u2=u(n+1:2*n,1);
st=Yeff*sqrt(u1.^2+u2.^2);
w=0.5*(Yeff^2/Ys)*(u1.^2+u2.^2);
pg=pdegplot(g2);
set(pg,'Color','w');
hold on;
pdeplot(p,e,t,'xydata',st,'colormap',jet,'mesh','off','flowdata',[u1 u2],'xygrid','off','colorbar','on');
set(gca,'DataAspectRatio',[1 1 1]);
%set(gcf,'Color',[0 0 0.5625]);
axis off
set(findobj('Type','line'),'Color','k');
figure(gcf)
%caxis([0 700])
%strain_energy;
