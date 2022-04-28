function [stt,snt,snn]=rotate_stress(sxx,sxy,syy,theta)
% returns stress tensor components in orthogonal (t,n)
% coordinate system at angle theta with respect to 
% original (x,y) coordinate system (t=tangent, n=normal)
  
c = cos(theta); s = sin(theta);
  
stt = sxx.*c.^2 + syy.*s.^2 + 2*sxy.*s.*c;
snn = sxx.*s.^2 + syy.*c.^2 - 2*sxy.*s.*c;
snt  = (syy-sxx).*s.*c + sxy.*(c.^2-s.^2);