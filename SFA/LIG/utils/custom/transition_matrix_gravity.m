function phi = transition_matrix_gravity(wo, T)

A = -vp(wo).*T;
s = norm(wo)*T;

if (s == 0)
    phi = eye(3);
else
    phi = eye(3) + A.*sinc(s/pi) + (A^2)*(1-cos(s))/(s^2);
end