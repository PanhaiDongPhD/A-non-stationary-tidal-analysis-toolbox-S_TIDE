function [S] = l_calculate_coefficient(IDP,ntime)
% computes the linear interpolation weights
% Please citing:
% Guo, Z., H. Pan, A. Cao, and X. Lv (2018), A harmonic analysis method adapted to capturing slow variations of tidal amplitudes and phases, 
%                                             Cont. Shelf Res., 164(June), 37¨C44, doi:10.1016/j.csr.2018.06.005.
S = zeros(ntime,IDP,'double'); 
if IDP>1
r = (ntime-1)/(IDP-1);   
for t=1:ntime
    k=floor((t-1)/r);
    if k==IDP-1
        k=k-1;
    end
    xi=k*r+1;xi1=(k+1)*r+1;
    k1=k+1;k2=k+2;
    S(t,k1)=1-(t-xi)/r;
    S(t,k2)=(t-xi)/r;   
end
else
   S(:,:)=1; 
end