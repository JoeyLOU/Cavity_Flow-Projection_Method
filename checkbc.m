function [ns_errors,bd_errors] = checkbc(x,uv,psi,t_pts,b_pts,l_pts,r_pts,m)
% Error calculations
dy_t_error=norm(uv(t_pts)-16*x.^2.*(1-x).^2,inf)/norm(uv(1:m),inf);
dy_b_error=norm(uv(b_pts),inf)/norm(uv(1:m),inf);
dx_l_error=norm(uv(l_pts+m),inf)/norm(uv(1+m:2*m),inf);
dx_r_error=norm(uv(r_pts+m),inf)/norm(uv(1+m:2*m),inf);
ns_errors=[dy_t_error,dy_b_error,dx_l_error,dx_r_error];

t_psi=norm(psi(t_pts),inf)/norm(psi,inf);
b_psi=norm(psi(b_pts),inf)/norm(psi,inf);
l_psi=norm(psi(l_pts),inf)/norm(psi,inf);
r_psi=norm(psi(r_pts),inf)/norm(psi,inf);
bd_errors=[t_psi,b_psi,l_psi,r_psi];
end