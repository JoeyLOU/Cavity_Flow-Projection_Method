function f = adv(m,uv,DX,DY)
% advection term, uv need to be specified
% stacked form [f_u; f_v]
u = uv(1:m);
v = uv(m+1:2*m);
f1 = u.*DX*u+v.*DY*u;
f2 = u.*DX*v+v.*DY*v;
f = [f1;f2];
end