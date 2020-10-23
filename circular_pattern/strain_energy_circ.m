F=0;SE=0;area=0;
for trianglecounter=1:size(t,2)

   % vertex numbers:
   v1n=t(1,trianglecounter);
   v2n=t(2,trianglecounter);
   v3n=t(3,trianglecounter);

   % vertex coordinates:
   v1=p(:,v1n);
   v2=p(:,v2n);
   v3=p(:,v3n);

   A=abs(v1(1,1)*(v2(2,1)-v3(2,1)) + v2(1,1)*(v3(2,1)-v1(2,1)) + v3(1,1)*(v1(2,1)-v2(2,1)) )/2;

   um1=(w(v1n)+w(v2n)+w(v3n))/3;
   um2=(st(v1n)+st(v2n)+st(v3n))/3;
  
   
   SE=SE+um1*A;
   F=F+um2*A;
   area=area+A;
   
  

end
StrainEnergy=SE;
StrainEnergy

%AvgSt=F/area
