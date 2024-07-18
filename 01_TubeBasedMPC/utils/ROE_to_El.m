function el = ROE_to_El(ROE, el_ref)

a = el_ref(1);
e = el_ref(2);
i = el_ref(3);
W = el_ref(4);
w = el_ref(5);
M = el_ref(6);

ex = e*cos(w);
ey = e*sin(w);

da = ROE(1);
dl = ROE(2);
dex = ROE(3);
dey = ROE(4);
dix = ROE(5);
diy = ROE(6);

id = dix + i;
Wd = W + diy/(sin(i));
ad = da*a + a;
ud = dl - diy + w + M;
exd = dex + ex;
eyd = dey + ey;

wd = atan2(eyd, exd);
ed = exd/cos(wd);
Md = ud - wd;

el = [ad; ed; id; Wd; wd; Md];


end