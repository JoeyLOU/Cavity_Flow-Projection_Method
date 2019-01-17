function si = w_to_psi(w,L,N,DY,x,bd_pts,t_pts,cond)
% computes psi from given w, used for startup
% impose bc on L
w = sparse(-w(:,end));
if cond
    L(bd_pts,:) = sparse(4*N,(N+1)^2);
    L(bd_pts,bd_pts) = eye(4*N);
    L(t_pts,:) = DY(t_pts,:);                   % set Newmann bc
    w(bd_pts) = sparse(4*N,1);                  % zeros at x=+-1, y=+-1
    w(t_pts) = 16*x.^2.*(1-x).^2;               % enforce bc at y=1
    % solve laplacian for psi
    si = full(L\w);
    
else
    L(bd_pts,:) = sparse(4*N,(N+1)^2);
    L(bd_pts,bd_pts) = eye(4*N);
    w(bd_pts) = sparse(4*N,1);
    si = full(L\w);
end
end