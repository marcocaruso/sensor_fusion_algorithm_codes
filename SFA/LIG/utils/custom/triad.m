function q = triad(a, g, b, r)

r1 = g;
r2 = vp(g, r)./norm(vp(g, r));
r3 = vp(r1, r2);

s1 = a;
s2 = vp(a, b)./norm(vp(a, b));
s3 = vp(s1, s2);

Mref = [r1 r2 r3];
Mobs = [s1 s2 s3];

q = dcm2quat(Mref*Mobs');
q = [-q(2:4)'; q(1)];
