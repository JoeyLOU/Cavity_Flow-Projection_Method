function w = uvtow(uv,DX,DY,m)
% calculate omega from u,v
u = uv(1:m);
v = uv(m+1:2*m);
w = DX*v-DY*u;
end