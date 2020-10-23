%Shiladitya Banerjee[UCL]
%Edited 01.02.2017

%Compute stresses in stadium geometry
%Unit of length in microns, forces in N
Ar = 8000;% area in microns squared

%define geometry
r=37.5E-6;% radius of curvature
L=round(((Ar - pi*(r*10^6)^2)/(2*r)),1)*10^-12;% length

ar=pi*r^2 + 2*r*L; % area
aspect=(L+2*r)/(2*r); % aspect ratio
Pr=2*pi*r + 2*L; % perimeter

%define pde geometry
g1=[1 -L/2 0 r 0 0 0 0 0 0];
g2=[3 4 -L/2 L/2 L/2 -L/2 -r -r r r];
g3=[1 L/2 0 r 0 0 0 0 0 0];
gdm=[g1; g2; g3];
[g,bt]=decsg(gdm','C1+R1+C2',['C1';'R1';'C2']');
[g1,bt1]=csgdel(g,bt);

%define mesh size
hmax=r/20;
%create initial mesh
[p e t]=initmesh(g1,'Hmax',hmax);
pg=pdegplot(g1);

%use the following commands to refine mesh
%[p e t]=initmesh(g);
%[p e t]=refinemesh(g,p,e,t);
%[p e t]=refinemesh(g,p,e,t);
%[p e t]=refinemesh(g,p,2e,t);


%Define mechanical parameters

E=5.8E3; % Youngs modulus in Pascals 
nu=0.43; % Poisson ratio
G=E/(2*(1+nu)); %shear modulus
mu=2*G*nu/(1-2*nu); % bulk modulus 

ht=8E-6; % height of the monolayer
G_sub=2.8E3;
Ys=G_sub/ht;% substrate stiffness (substrate shear modulus / height) units [N/m^3]
Ya=1E9; % adhesion stiffness [N/m^3]
Yeff=1/((1/Ys)+(1/Ya)); % effective stiffness of adhesion
sigma=5.8E2; %bulk active stress in pascals 
sigma2=3E-4; % line tension in units [pascals X meters]


% Elastic modulus matrix
c=[2*G+mu;0;0;G;0;mu;G;0;0;G;mu;0;G;0;0;2*G+mu];

%boundary condition
b=@(p,e,u,time) boundaryFile(p,e,Pr,sigma,sigma2,ht);

%effective driving force
a=[Yeff/ht;Yeff/ht];
f=[0;0];

%compute solution
u=assempde(b,p,e,t,c,a,f);

% use the following for nonlinear problems
% u = pdenonlin(b,p,e,t,c,a,f,'jacobian','lumped');

% output mechanical quantities
n=size(p,2);
u1=u(1:n,1); % x-displacement 
u2=u(n+1:2*n,1); % y-displacement
st=Yeff*sqrt(u1.^2+u2.^2); % traction stress
w=0.5*(Yeff^2/Ys)*(u1.^2+u2.^2); % strain energy density

%visualize
pg=pdegplot(g1);
%set(pg,'Color','k');
%hold on;
pdeplot(p,e,t,'xydata',st,'colormap',jet,'mesh','off','flowdata',[u1 u2],'xygrid','off','colorbar','on');
set(gca,'DataAspectRatio',[1 1 1],'XLim',[-150E-6 150E-6],'YLim',[-100E-6 100E-6]);
%set(gcf,'Color',[0 0 0.5625]);
set(findobj('Type','line'),'Color','w');
axis off
caxis([0 400])


strain_energy_stad; % calculates strain energy
maxstress=max(st) % calculates maximum traction stress
