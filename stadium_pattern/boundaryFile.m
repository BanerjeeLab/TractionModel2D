
function[q,gmatrix,h,rmatrix]=boundaryFile(p,e,Pr,sigma,sigma2,ht)
N=2;
ne=size(e,2);
q = zeros(N^2,ne);
gmatrix = zeros(N,ne);
h = zeros(N^2,2*ne);
rmatrix = zeros(N,2*ne);

for i=1:ne
    x1=p(1,e(1,i)); x2=p(1,e(2,i)); 
    y1=p(2,e(1,i)); y2=p(2,e(2,i)); 
    
    if i==1
        x0=p(1,e(1,ne));y0=p(2,e(1,ne));
    else
        x0=p(1,e(1,i-1));y0=p(2,e(1,i-1));
    end   
    if i==ne
        x3=p(1,e(2,ne));y3=p(2,e(2,ne));
    else
        x3=p(1,e(2,i+1));y3=p(2,e(2,i+1));
    end
    
    xm=(x1+x2)/2; xl=(x0+x1)/2; xr=(x2+x3)/2;
    ym=(y1+y2)/2; yl=(y0+y1)/2; yr=(y2+y3)/2;
    si=sqrt((x2-x1)^2+(y2-y1)^2);
    si1=sqrt((x1-x0)^2+(y1-y0)^2);
    s=(si+si1)/2;
    txi=(x2-x1)/si; tyi=(y2-y1)/si;
    txi1=(x1-x0)/si1; tyi1=(y1-y0)/si1;
    nxi=(y2-y1)/sqrt((x2-x1)^2+(y2-y1)^2); nyi=-(x2-x1)/sqrt((x2-x1)^2+(y2-y1)^2); nxi1=tyi1; nyi1=-txi1;
    ki=sqrt((txi-txi1)^2+(tyi-tyi1)^2)/s;
    gi=zeros(N,1);
    gi(1)=-(sigma+sigma2*(Pr*ki)/ht)*nxi;
    gi(2)=-(sigma+sigma2*(Pr*ki)/ht)*nyi;
    gmatrix(:,i)=gi;
end

end