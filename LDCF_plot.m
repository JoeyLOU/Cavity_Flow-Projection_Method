function LDCF_plot(filename,node_start,node_final,plotincre,incre)
% plot function for LDCF
clc; close all;
% filename = '49grid/49grid';
% incre = 100;
% node_start = 1200; node_final = 1310;
load([filename,num2str(node_start),'.mat']);
%% helper parameters
% Matrix constructions
[D,x] = cheb(N); x = (x+1)/2; y=x;                  % chebyshev D
dx = zeros(size(x));
for i = 2:length(x)-1
    dx(i) = min(x(i)-x(i+1),x(i-1)-x(i));
end
dx([1,length(x)]) = x(1)-x(2);
[dxx, dyy] = meshgrid(dx,dx); dxx = dxx(:); dyy = dyy(:);
D = D*2;
[xx,yy] = meshgrid(x,y); xx=xx(:); yy=yy(:);
xxx=reshape(xx,N+1,N+1); yyy=reshape(yy,N+1,N+1);   % plotting grid
% find boundary points
bd_pts = find(abs(xx)==0 | abs(yy)==0 | abs(xx)==1 | abs(yy)==1);% boundary points
t_pts = find(yy==1);                               % top points
l_pts = find(xx==0);                               % left points
r_pts = find(xx==1);                               % right points
b_pts = find(yy==0);                               % bottom points
i_pts = find(abs(xx)~=1 & abs(yy)~=1 & abs(xx)~=0 & abs(yy)~=0);
mid_val=x(floor(N/2)+1);
xmid_pts = find(xx==mid_val);
ymid_pts = find(yy==mid_val);

% helper matrices
D2 = D^2;
I = eye(N+1);
L = kron(D2,I) + kron(I,D2);                % laplacian
L2 = [L,   zeros(m,m)
    zeros(m,m)   L];                             % big laplacian
DX = kron(D,I);                             % Dx
DY = kron(I,D);                             % Dy


%% Plot
for node=node_start:plotincre:node_final
    load([filename,num2str(node),'.mat']);
    for i=2:2
        iter = node*incre;
        time = iter*dt;
        uv_temp = uv(:,i);
        % p_temp = p(:,i);
        w = uvtow(uv_temp,DX,DY,m);
        si = w_to_psi(w,L,N,DY,x,bd_pts,t_pts,0);
        figure(1)
        subplot(2,2,1)
        sii=reshape(si,N+1,N+1);
        contourf(xxx,yyy,sii)
        xlabel x, ylabel y, title(['\psi at t= ',num2str(time)])
        colorbar
        
        subplot(2,2,2)
        ww=reshape(w,N+1,N+1);
        contourf(xxx,yyy,ww)
        xlabel x, ylabel y, title(['w at t= ',num2str(time)])
        colorbar
        
        subplot(2,2,3)
        plot(xx(ymid_pts),uv_temp(m+ymid_pts),'-x');
        grid on
        xlabel x, ylabel y, title(['v(x,0.5) at t= ',num2str(time)])
        
        subplot(2,2,4)
        plot(uv_temp(xmid_pts),yy(xmid_pts),'-x');
        grid on
        xlabel x, ylabel y, title(['u(0.5,y) at t= ',num2str(time)])
        
        fprintf('At time %5.3f\n',time)
        [ns_errors,~] = checkbc(x,uv_temp,si,t_pts,b_pts,l_pts,r_pts,m);
        fprintf('u_t Error = %5.3e\tu_b Error = %5.3e\n',ns_errors(1:2))
        fprintf('v_l Error = %5.3e\tv_r Error = %5.3e\n',ns_errors(3:4))
        
        % divergence norm calculation
        div = sum((DX(i_pts,:)*uv_temp(1:m)+DY(i_pts,:)*uv_temp(m+1:2*m)).^2)...
            /(N-1)^2;
        CFL=dt*(uv_temp(1:m)./dxx + uv_temp(m+1:2*m)./dyy);
        fprintf('Div Error = %5.3e\tCFL = %6.4e\n',div,norm(CFL,inf))
        M1 = norm(si(i_pts),inf);
        [a,b] = find(sii == -M1);
        M1_pos = [x(a),y(b)];
        M2 = norm(w(t_pts),inf);
        fprintf('M1 = %5.4e\t\tM2 = %5.4e\n',M1,M2)
        fprintf('s1 = (%5.3f, %5.3f)\n\n',M1_pos)
        
        pause(0.1)
    end
end