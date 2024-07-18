function T = ECI_2_LVLH_T(x)

r = x(1:3,1);
v = x(4:6,1);

ir = r/norm(r);
h =  cross(r,v);
ih = h/norm(h);
iv = cross(ih,ir);
T = [ir'; iv'; ih'];