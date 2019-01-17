function p = psolver(temp,DX,DY,L,m,N,dt,t_pts,b_pts,l_pts,r_pts,bd_pts)
% update pressure term
rhs = (DX*temp(1:m)+DY*temp(m+1:2*m))/dt;
L([t_pts; b_pts],:) = DY([t_pts; b_pts],:);
L([l_pts; r_pts],:) = DX([l_pts; r_pts],:);
rhs(bd_pts) = zeros(4*N,1);
p = full(sparse(L)\sparse(rhs));
end