function z = vp(x, y)

z = [  0    -x(3)   x(2);
      x(3)    0    -x(1);
     -x(2)   x(1)    0   ];

if nargin > 1, z = z*y; end
