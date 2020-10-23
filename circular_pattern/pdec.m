%clear
%geometry
r=26E-6;
gdm=[1 0 0 r];
[g,bt]=decsg(gdm','C1',['C1']');
[g1,bt1]=csgdel(g,bt);
hmax=r/25;
[p e t]=initmesh(g1,'Hmax',hmax);
pg=pdegplot(g1);
%[p e t]=initmesh(g);
%[p e t]=refinemesh(g,p,e,t);
%[p e t]=refinemesh(g,p,e,t);
%[p e t]=refinemesh(g,p,2e,t);
ar=pi*r^2;
aspect=1;
Pr=2*pi*r;

%material properties
E=0.54E4; % Young's modulus
nu=0.445; % Poisson ratio
G=E/(2*(1+nu)); %shear modulus
mu=2*G*nu/(1-2*nu);
ht=5E-6;
Ys=1E9;% substrate rigidity
Ya=8E8;
Yeff=1/((1/Ys)+(1/Ya));
sigma=1E3; %bulk active stress
sigma2=0.75E-3;%line tension


c=[2*G+mu;0;0;G;0;mu;G;0;0;G;mu;0;G;0;0;2*G+mu];

%boundary condition

b=@(p,e,u,time) boundaryFile(p,e,Pr,sigma,sigma2,ht);
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
pg=pdegplot(g1);
%set(pg,'Color','k');
%hold on;
pdeplot(p,e,t,'xydata',st,'colormap',parula,'mesh','off','flowdata',[u1 u2],'xygrid','off','colorbar','on');
set(gca,'DataAspectRatio',[1 1 1],'XLim',[-60E-6 60E-6],'YLim',[-30E-6 30E-6]);
%set(gcf,'Color',[0 0 0.5625]);
set(findobj('Type','line'),'Color','w');
axis off
%caxis([0 1000])
strain_energy_circ;
