clc
clear all
close all

%% Datos

Re = 100;
Nx = 64; % Celdillas en x
Ny = 64; % Celdillas en y
Lx = 1;
Ly = 1;
tf = 20;
CFL = 0.5;
Ulid = 1;

%% Staggered Grid Generation

[grid, u, v, p] = gridGeneration(Lx, Ly, Nx, Ny);
dt = CFL*min(grid.cellMin^2*Re, grid.cellMin);

%% Boundary conditions

bc.uS = zeros(1,Nx-1);
bc.uN = ones(1,Nx-1);
bc.uE = zeros(Ny,1);
bc.uW = zeros(Ny,1);

bc.vS = zeros(1,Nx);
bc.vN = zeros(1,Nx);
bc.vE = zeros(Ny-1,1);
bc.vW = zeros(Ny-1,1);

% figure,
% plot(grid.xc(1,1), grid.yc(1,1), 'kx', grid.xu(1,1), grid.yu(1,1), 'g>',...
%     grid.xv(1,1), grid.yv(1,1), 'r^')
% h = legend('$p$', '$u$', '$v$');
% set(h, 'interpreter', 'latex', 'fontsize', 14)
% %hold on
% plot(grid.x, grid.y, 'k.', grid.x, grid.y, 'k-', grid.y, grid.x, 'k-',...
%     grid.xc, grid.yc, 'kx', grid.xv, grid.yv, 'r^', grid.xu, grid.yu, 'g>')
% set(gca, 'fontname', 'times', 'fontsize', 12);


%% Operators

[D, G, R, M] = DGRM(grid, Nx, Ny);

Lhat = laplacianHat(Nx, Ny, grid);
L.L = M.hat*Lhat.L/R.R;

M.M = M.hat/R.R;
M.inv = inv(M.M);

Ahat = sparse(speye(size(Lhat.L))/dt - 0.5*Lhat.L/Re);
A = M.hat*Ahat/R.R;
dA = decomposition(A);

BN = dt*speye(size(M.M))/M.M + (0.5/Re)*dt*dt*(M.inv*L.L)*M.inv +...
    ((0.5/Re)^2)*(dt^3)*((M.inv*L.L)^2)*M.inv;

LHS = sparse(G.G'*BN*G.G);
dLHS = decomposition(LHS);

%% Simulation

u = reshape(u, [], 1);
v = reshape(v, [], 1);

uOld = u;
vOld = v;

t = 0;
k = 0;

tic
while t<= tf
    
    u = reshape(u, [], 1);
    v = reshape(v, [], 1);

    % Advective terms
    [NhatOld, ~, ~] = convectionHat(grid, uOld, vOld, Nx, Ny, bc);
    [Nhat, ua, va] = convectionHat(grid, u, v, Nx, Ny, bc);
    
    rnHat = explicitTerms(Lhat, Re, dt, Nhat, NhatOld, u, v);  
    rn = M.hat*rnHat;
        
    %% 1. Solve for intermediate velocity
       
    % BC's due to Laplacian
    bc1hat.u = Lhat.ux0*bc.uW + Lhat.ux1*bc.uE + Lhat.uy1*bc.uN' + ...
        Lhat.uy0*bc.uS';
    bc1hat.v = Lhat.vx0*bc.vW + Lhat.vx1*bc.vE + Lhat.vy1*bc.vN' + ...
        Lhat.vy0*bc.vS';

    bc1 = M.hat*[bc1hat.u; bc1hat.v]/Re;
    
    % Flux calculation
    
    q = dA\(rn + bc1);    
    qu = q(1:Ny*(Nx-1));
    qv = q(Ny*(Nx-1)+1:end);
       
  
    %% 2. Solve the Poisson Equation
    
    bc2 = D.uW*(bc.uW*grid.dXp(1)) + D.uE*(bc.uE*grid.dXp(end)) + ...
        D.vS*(bc.vS'*grid.dYp(end)) + D.vN*(bc.vN'*grid.dYp(1));

    RHS = G.G'*q + bc2;
    
    phi = dLHS\RHS;
    
    %% 3. Projection step
    
    q = q - BN*G.G*phi;
       
    vel = R.R\q;

    t = t + dt;
    k = k + 1;

    % Residuals
    %eps(k) = max(abs(C - B*phi));
    epsU(k) = max(abs(u - vel(1:Ny*(Nx-1))));
    epsV(k) = max(abs(v - vel(Ny*(Nx-1)+1:end)));

    % Separation of velocity components
    u = vel(1:Ny*(Nx-1));
    v = vel(Ny*(Nx-1)+1:end);
    
    % Information
    fprintf(['t = ' num2str(t) '. Elapsed time: ' num2str(toc) 's \n']);
    fprintf(['Residuals u = ' num2str(epsU(k)) '.\n'...
        'Residuals v = ' num2str(epsV(k)) '\n \n']);

%     if (mod(k, 0.25) == 0)
%         u = reshape(u, Ny, Nx-1);
% 
%         F(j) = getframe(u);
%         writeVideo, movie
%     
%     end
    
    
    
%     if (mod(k, 50) == 0)
%         u = reshape(u, Ny, Nx-1);
%         figure(2),
%         pcolor(grid.x, grid.y, sqrt(ua.^2 + va.^2)), shading interp
%         colorbar, colormap jet
%         title(['t = ' num2str(t)])
%         drawnow
%     end
    

end
toc

%% Plots
close all

u = reshape(u, Ny, Nx-1);
v = reshape(v, Ny-1, Nx);
phi = reshape(phi, Ny, Nx);
qu = reshape(qu, Ny, Nx-1);
qv = reshape(qv, Ny-1, Nx);

% Contours
figure(2),
pcolor(grid.x, grid.y, hypot(ua,va)), shading interp
colorbar, colormap jet
title(['$u$; $t = $' num2str(t)], 'interpreter', 'latex',...
    'fontsize', 20)
pbaspect([Lx Ly 1])

figure(3),
s = 2;
quiver(grid.x(1:s:end,1:s:end), grid.y(1:s:end,1:s:end),...
    ua(1:s:end,1:s:end), va(1:s:end,1:s:end), 8),
xlim([0 Lx]), ylim([0 Ly])
pbaspect([Lx Ly 1])

figure(4),
xstart = Lx*rand(Nx,1); 
ystart = Ly*rand(Ny,1); 

streamline(grid.x, grid.y, ua, va, xstart, ystart, [0.03, 3000]);
xlim([0 Lx]), ylim([0 Ly])
pbaspect([Lx Ly 1])

% Contours
figure(5),
contourf(grid.xp, grid.yp, phi), %shading interp
colorbar, %colormap jet
title(['$\phi$; $t = $' num2str(t)], 'interpreter', 'latex',...
    'fontsize', 20)
pbaspect([Lx Ly 1])

% Residuals
figure(6),
plot(1:k, epsU, 1:k, epsV)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
h = legend('$u$', '$v$');
set(h, 'interpreter', 'latex', 'fontsize', 16)
xlabel('$N$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\xi$', 'interpreter', 'latex', 'fontsize', 16)
title('Residuals')


%% Validation
close all

uData = load('uData');
vData = load('vData');

% Interpolation
yq = linspace(0, 1, 64);
xq = yq*0 + 0.5;

uInt = interp2(grid.x, grid.y, ua, xq, yq);
vInt = interp2(grid.x, grid.y, va, yq, xq);

% Comparison

figure(4),
plot(uData(:,1), uData(:,2), 'rs', 'markersize', 8), hold on,
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
xlabel('$y$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$u$', 'interpreter', 'latex', 'fontsize', 16)
plot(uInt, yq)
hold off

figure(5),
plot(vData(:,1), vData(:,2), 'rs', 'markersize', 8), hold on,
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$v$', 'interpreter', 'latex', 'fontsize', 16)
plot(yq, vInt)

