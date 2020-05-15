function [S] = sinc_calculate_coefficient(IDP,ntime)
% computes the sinc interpolation weights
% Please citing:
%Pan, H., 2018. Application Examples of S_TIDE Toolbox. Technical Report 2018-11. Key Laboratory of Physical Oceanography, Ocean University of China, Qingdao, China, 20pp. 
S = zeros(ntime,IDP,'double'); 
if IDP>2
r = (ntime-1)/(IDP-1);   
xi=zeros(1,IDP,'double');
for i=1:IDP
    xi(i)=(i-1)*r+1;
end
%for i=1:IDP
    %for t=1:ntime
    %S(t,i)=sinc((t-xi(i))/r); 
    %end
%end
t=1:ntime;t=t';
tt=repmat(t,1,IDP);   % the repmat function can greatly accelerate the calculation speed (matlab向量化操作，可以提高50倍以上的运算速度)
xii=repmat(xi,ntime,1);
S=sinc((tt-xii)/r);
else
    S= l_calculate_coefficient(IDP,ntime);
end
