function LDCF_ND(filename,flag,node,incre,N,Re,dt,tspan)
% LDCF using projection methods
% updated the x range to be [0, 1]
%% setup
% flag = 1;       % 1 for continuation, 0 for initialization
% node = 10;      % last run's node number, 0 for initialization
% incre = 100;
% filename = 'test/test';
% set grid size and Reynolds
if flag == 0     % initialize parameters
%     N = 16;                             % grid size
     m =(1+N)^2;                         % grid number
%     Re = 400;                           % reynolds
%     dt = 1e-2;                           % discretization in time
%     tspan = 500;                        % set time span
else
    load([filename,num2str(node),'.mat']);
end

% Matrix constructions
[D,x] = cheb(N); x = (x+1)/2; y=x;                  % chebyshev D
D = 2*D;
[xx,yy] = meshgrid(x,y); xx=xx(:); yy=yy(:);
xxx=reshape(xx,N+1,N+1); yyy=reshape(yy,N+1,N+1);   % plotting grid
% find boundary points
bd_pts = find(abs(xx)==0 | abs(yy)==0 | abs(xx)==1 | abs(yy)==1);% boundary points
t_pts = find(yy==1);                                % top points
l_pts = find(xx==0);                               % left points
r_pts = find(xx==1);                                % right points
b_pts = find(yy==0);                               % bottom points
% helper matrices
D2 = D^2;
I = eye(N+1);
L = kron(D2,I) + kron(I,D2);                % laplacian
L2 = [L,   zeros(m,m)
    zeros(m,m)   L];                             % big laplacian
DX = kron(D,I);                             % Dx
DY = kron(I,D);                             % Dy
DD = [DX, zeros(m,m)
    zeros(m,m), DY];
% LHS
lhs = [eye(m)/dt-L/2/Re, zeros(m,m)
        zeros(m,m),       eye(m)/dt-L/2/Re];

% dirichlet boundary condition on u and v
% clear rows
lhs([bd_pts;bd_pts+m],:) = zeros(8*N,2*m);
% apply bc on u and v
lhs([bd_pts;bd_pts+m],[bd_pts;bd_pts+m]) = eye(8*N);

LHS = sparse(lhs);

%% startup
if flag == 0
    uv = zeros(2*m,1);      % uv0
    f = adv(m,uv,DX,DY);    % f0
    p = zeros(m,1);         % p0
    % forward euler
    rhs0 = (eye(m*2)/dt+L2/2/Re)*uv-f;
    % boundary conditions
    rhs0([bd_pts;bd_pts+m]) = zeros(8*N,1);
    % u on top and bottom (tangent to boundary)
    rhs0(t_pts) = 16*x.^2.*(1-x).^2 + dt*DX(t_pts,:)*p;
    rhs0(b_pts) = dt*DX(b_pts,:)*p;
    % v on left and right (tangent to boundary)
    rhs0([l_pts;r_pts]+m) = dt*DY([l_pts;r_pts],:)*p;
    
    RHS0 = sparse(rhs0);
    temp = LHS\RHS0;
    p = [p, psolver(temp,DX,DY,L,m,N,dt,t_pts,b_pts,l_pts,r_pts,bd_pts)];
    temp = temp - dt*[DX*p(:,end); DY*p(:,end)];
    uv = [uv, temp];
    f = [f(:,end), adv(m,uv(:,end),DX,DY)];
    
    iter = 1;
    conv = 1;
end

t = iter*dt;
%% iterate
while t< tspan && conv(end)>1e-8
    % RHS
    rhs = (eye(m*2)/dt+L2/2/Re)*uv(:,end) - 3/2*f(:,end) + 1/2*f(:,end-1);
    % boundary conditions
    rhs([bd_pts;bd_pts+m]) = zeros(8*N,1);
    % u on top and bottom (tangent to boundary)
    rhs(t_pts) = 16*x.^2.*(1-x).^2 + ...
        dt*(2*DX(t_pts,:)*p(:,end)-DX(t_pts,:)*p(:,end-1));
    rhs(b_pts) = dt*(2*DX(b_pts,:)*p(:,end)-DX(b_pts,:)*p(:,end-1));
    % v on left and right (tangent to boundary)
    rhs([l_pts+m;r_pts+m]) = dt*(2*DY([l_pts;r_pts],:)*p(:,end)-DY([l_pts;r_pts],:)*p(:,end-1));

    rhs = sparse(rhs);
    temp = full(lhs\rhs);
    % update pressure corrective term
    p = [p(:,end), psolver(temp,DX,DY,L,m,N,dt,t_pts,b_pts,l_pts,r_pts,bd_pts)];
    temp = temp - dt*[DX*p(:,end); DY*p(:,end)];
    uv = [uv(:,end), temp];                                           % update uv
    f = [f(:,end), adv(m,uv(:,end),DX,DY)];
    conv = norm(uv(:,end)-uv(:,end-1),inf)/norm(uv(:,end),inf)/dt;
    iter = iter + 1;
    t = iter*dt;
    if mod(iter,incre)==0
        % save file and display convergence
        node = node + 1;
        fprintf('Iter %i, convergence error = %5.3e\n',iter,conv(end))
        save([filename,num2str(node),'.mat'],'uv','f','p','conv','N','m','Re','dt','tspan','iter','incre');
    end
end

end