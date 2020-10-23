%units : elastic constant (Pa), length (microns), stress (Pa)

%geometry
filename='fibroblast3_pts_byhand.txt';

% read data and convert distances to microns
data=dlmread(filename);
conv=0.107*1E-6; % conversion factor
xpts=conv*data(:,1);
ypts=conv*data(:,2);

% creating a finite element mesh
nm=numel(xpts);
gdm=[2 nm xpts' ypts'];
[g1,bt]=decsg(gdm','P1',['P1']');
[g2,bt1]=csgdel(g1,bt);
[p e t]=initmesh(g2);
%hmax=1/16;
%[p e t]=initmesh(g1,'Hmax',hmax);
%[p e t]=refinemesh(g2,p,e,t);
%[p e t]=refinemesh(g,p,e,t);



%material properties
E=0.54E4; % Young's modulus
nu=0.43; % Poisson ratio
G=E/(2*(1+nu)); %shear modulus
mu=2*G*nu/(1-2*nu);
ht=3E-6; % cell height
Ys=2.9E9;% substrate rigidity
Ya=1E9;% adhesion
Yeff=1/((1/Ys)+(1/Ya));
sigma=2.4E3; %active stress
sigma2=0.7E-3;%line tension
Pr=235.89E-6; % contour length

c=[2*G+mu;0;0;G;0;mu;G;0;0;G;mu;0;G;0;0;2*G+mu];

b=@(p,e,u,time) boundaryFile(p,e,Pr,sigma,sigma2,ht);

a=[Yeff/ht;Yeff/ht];
f=[0;0];

%solution
u=assempde(b,p,e,t,c,a,f);
%u = pdenonlin(b,p,e,t,c,a,f,'jacobian','lumped');
%visualize
n=size(p,2);
u1=u(1:n,1); % x-component of displacement vector
u2=u(n+1:2*n,1); % y-component of displacement vector
st=Yeff*sqrt(u1.^2+u2.^2); % traction stress magnitude
w=0.5*(Yeff^2/Ys)*(u1.^2+u2.^2); % strain energy density
pg=pdegplot(g2);
pdeplot(p,e,t,'xydata',st,'colormap',jet,'mesh','off','flowdata',[u1 u2],'xygrid','off','colorbar','off');
set(gca,'DataAspectRatio',[1 1 1]);
set(gcf,'Color',[0 0 0.5625]);
axis off
set(findobj('Type','line'),'Color','w');
caxis([0 2500])
%integral; % calculate strain energy and average traction stress


