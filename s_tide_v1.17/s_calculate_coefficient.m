
function [S] = s_calculate_coefficient(IDP,ntime)
% computes the spline interpolation coefficient(weight)
% written by Haidong Pan
% thanks to Dr. Yuzhe Wang for pointing out a small error in this function
S = zeros(ntime,IDP,'double'); 
if IDP>2
r = (ntime-1)/(IDP-1);           
MM = zeros(IDP,IDP,'double');            
C = zeros(IDP, IDP,'double');
for i=2:IDP-1
    C(i, i-1)  = 0.5;
    C(i, i)     = 2.0;
    C(i, i+1) = 0.5;
end
C(1, 1) = 2;
C(1, 2) = 1;
C(end, end-1) = 1;
C(end, end) = 2.0;   
C_ni = pinv(C);  
 b = zeros(IDP, IDP,'double');
 for i=2:IDP-1
    b(i,i+1) = 1;
    b(i,i-1) = -1;
 end
 b(1,1) = -2;
 b(1,2) = 2;
 b(end, end-1) = -2;
 b(end, end) = 2;

 MM = C_ni*b*(1.5/r);

for t=1:ntime
    k=floor((t-1)/r);
    if k==IDP-1
        k=k-1;
    end
    xi=k*r+1;xi1=(k+1)*r+1;
    k1=k+1;k2=k+2;
    tmp1=(t-xi)*((t-xi1)/r)^2;
    tmp2=(t-xi1)*((t-xi)/r)^2;
    tmp3=(1+2*(t-xi)/r)*((t-xi1)/r)^2;
    tmp4=(1-2*(t-xi1)/r)*((t-xi)/r)^2;
    S(t,:)=tmp1*MM(k1,:)+tmp2*MM(k2,:);
    S(t,k1)=S(t,k1)+tmp3;
    S(t,k2)=S(t,k2)+tmp4;
    
end
else
       S= l_calculate_coefficient(IDP,ntime);
end 

