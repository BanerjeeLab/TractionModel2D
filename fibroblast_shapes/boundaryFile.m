
function[q,gmatrix,h,rmatrix]=boundaryFile(p,e,Pr,sigma,sigma2,ht)
N=2;
ne=size(e,2);
q = zeros(N^2,ne);
gmatrix = zeros(N,ne);
h = zeros(N^2,2*ne);
rmatrix = zeros(N,2*ne);
sml=1E-14;
nvm = zeros(N,ne);%initialize normal vector matrix
nvm2 = zeros(N,ne);
midpts = zeros(N,ne);%initialize edge midpoint matrix
%midpts2=zeros(2,N);
km = zeros(1,N);% initialize curvature matrix
km2 = zeros(1,N);
%


% reorder edges into succesive order - e2 is the new matrix with this
% pieces
e2 = zeros(size(e));
e2(:,1) = e(:,1);
trigger = 1;
for i = 2:ne
    if trigger == 1
        pt = find(e(1,:) == e2(2,i-1));
    else
        pt = find(e(2,:) == e2(1,i-1));
    end

    if isempty(pt)
        
        
        
        if trigger == 1
            pt = find(e(2,:) == e2(2,i-1));
            newcol = e(:,pt(1));
            if newcol == e2(:,i-1)
                newcol = e(:,pt(2));
            end
            e2(:,i) = newcol;
            trigger = 0;
        else
            pt = find(e(1,:) == e2(1,i-1));
            newcol = e(:,pt(1));
            if newcol == e2(:,i-1)
                newcol = e(:,pt(2));
            end
            e2(:,i) = newcol;
            trigger = 1;
        end
    else
        e2(:,i) = e(:,pt);
    end
    
end


for i=1:ne
    %endpoints of the ith edge
    x1=p(1,e2(1,i)); y1=p(2,e2(1,i));
    x2=p(1,e2(2,i)); y2=p(2,e2(2,i)); 
    %midpoints of the ith edge
    xm=(x1+x2)/2;
    ym=(y1+y2)/2;
    
    %[xl;yl] are edge coordinates left of [xm;ym]
    if i==1
        xl=(p(1,e2(2,ne))+p(1,e2(1,ne)))/2;
        yl=(p(2,e2(2,ne))+p(2,e2(1,ne)))/2;
    else
        xl=(p(1,e2(2,i-1))+p(1,e2(1,i-1)))/2;
        yl=(p(2,e2(2,i-1))+p(2,e2(1,i-1)))/2;
    end   
    
    %[xr;yr] are edge coordinates right of [xm;ym]
    if i==ne
        xr=(p(1,e2(1,1))+p(1,e2(2,1)))/2;
        yr=(p(2,e2(1,1))+p(2,e2(2,1)))/2;
    else
        xr=(p(1,e2(1,i+1))+p(1,e2(2,i+1)))/2;
        yr=(p(2,e2(1,i+1))+p(2,e2(2,i+1)))/2;
    end
   
    pm=[xm; ym];pl=[xl; yl];pr=[xr; yr];
    

    
    %fit a parabola through the points pm, pl and pr and compute tangent vectors
    s12=norm(pm-pl)+sml;
    s23=norm(pr-pm)+sml;
    txi=((s23^2)*(xm-xl)+(s12^2)*(xr-xm))/(s12*s23*(s12+s23));
    tyi=((s23^2)*(ym-yl)+(s12)^2*(yr-ym))/(s12*s23*(s12+s23));
    
    tv=[txi;tyi]; % tangent vector
    
    nv=-[tyi;-txi];% unit normal
    tv=tv/norm(tv); % normalize tangent vector
    nv=nv/norm(nv);
    nvm(:,i) = nv; % your normal vector calculation always points outward
    ki=-2*((xm-xl)*(yr-ym)-(ym-yl)*(xr-xm))/sqrt(((xm-xl)^2+(ym-yl)^2)*((xr-xm)^2+(yr-ym)^2)*((xl-xr)^2+(yl-yr)^2));
    km(1,i)=ki;
    midpts(:,i)=pm;
    pts=find(e(1,:) == e2(1,i) & e(2,:) == e2(2,i));
    %midpts2(:,pts)=[(p(1,e(1,pts))+ p(1,e(2,pts)))/2 ; (p(2,e(1,pts))+ p(2,e(2,pts)))/2];
    km2(1,pts)=km(1,i);
    nvm2(:,pts)=nv;
    
    
end





for i=1:ne
    
    gi=zeros(N,1);
    gi(1)=-(sigma+sigma2*(Pr*km2(1,i)/ht))*nvm2(1,i);
    gi(2)=-(sigma+sigma2*(Pr*km2(1,i)/ht))*nvm2(2,i);
    gmatrix(:,i)=gi;
    
end
    
    
    
end