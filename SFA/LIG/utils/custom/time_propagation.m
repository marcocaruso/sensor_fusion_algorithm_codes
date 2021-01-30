function [xproj, Pproj] = time_propagation(xupdt, Pupdt, o, stdg, T, B, Q)
       
A     = transition_matrix_gravity(o, T);
TP    = blkdiag(A, B);
xproj = TP*xupdt;

sigmaw = eye(3)*(stdg^2*T);
gupdt  = xupdt(1:3);
gamma  = blkdiag(sigmaw, eye(3));
PP     = blkdiag(vp(gupdt), Q);

Pproj  = TP*Pupdt*TP' + PP*gamma*PP';