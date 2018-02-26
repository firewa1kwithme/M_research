clear;

% улучшенный метод, работает быстрее, без модулей
global T_x T_y bc pw1 pw2
global Nx Ny hx hy
Nx = 20;
Ny = 20;   
hx = 1;   
hy = 1;   
N = 1000;

sle0 = zeros(Nx-2,Ny-2);
sle1 = zeros(Nx-2,Ny-2);
sle_t = zeros(Nx-2,Ny-2);
sle_bc = zeros(Nx-2,Ny-2);
sle_pw1 = zeros(Nx-2,Ny-2);
sle_pw2 = zeros(Nx-2,Ny-2);
    
D = zeros(Nx-2,Ny-2);
DyT = zeros(Nx-2,Ny-2);
DzT = zeros(Nx-2,Ny-2);
DyBC = zeros(Nx-2,Ny-2);
DzBC = zeros(Nx-2,Ny-2);
DyPW1 = zeros(Nx-2,Ny-2);
DzPW1 = zeros(Nx-2,Ny-2);
DyPW2 = zeros(Nx-2,Ny-2);
DzPW2 = zeros(Nx-2,Ny-2);
    
RES = zeros(Nx-2, (Nx-2)*5, N);
 
tic
for i = 1:N
    
 % реализации систем линейных уравнений с различными входными данными
    logNormal();
    upscaling();
        
    Tx_1 = T_x; % transmissibility
    Ty_1 = T_y;
    bc_1 = bc; % boundary conditions
    pw1_1 = pw1; % давление в центральной скважине
    pw2_1 = pw2; % давление в другой скважине
        
    logNormal();
    upscaling();
        
    sle0 = sle(Tx_1, Ty_1, bc_1, pw1_1, pw2_1); % f(x)
    sle1 = sle(T_x, T_y, bc, pw1, pw2); % f(x')
    sle_t = sle(Tx_1, Ty_1, bc, pw1, pw2); % f(y,z') with respect to t
    sle_bc = sle(T_x, T_y, bc_1, pw1, pw2); % f(y,z') with respect to bc
    sle_pw1 = sle(T_x, T_y, bc, pw1_1, pw2); % f(y,z') with respect to pw-s
    sle_pw2 = sle(T_x, T_y, bc, pw1, pw2_1); % f(y,z') with respect to pw2
    
    g = Nx-2;
    RES(:,1:g,i)=sle0;
    RES(:,g+1:2*g,i)=sle1;
    RES(:,2*g+1:3*g,i)=sle_t;
    RES(:,3*g+1:4*g,i)=sle_bc;
    RES(:,4*g+1:5*g,i)=sle_pw1;
    RES(:,5*g+1:6*g,i)=sle_pw2;
    
end
   
f_0=sum(RES,3)./N;
        
for iNA = 1:N
    
    D = D + (RES(:,1:g,iNA)-f_0(:,1:g)).*(RES(:,1:g,iNA)-f_0(:,1:g));
    DyT = DyT+(RES(:,1:g,iNA)-f_0(:,1:g)).*(RES(:,2*g+1:3*g,iNA)-f_0(:,2*g+1:3*g));
    DzT = DzT+(RES(:,g+1:2*g,iNA)-f_0(:,g+1:2*g)).*(RES(:,2*g+1:3*g,iNA)-f_0(:,2*g+1:3*g));
    
    DyBC = DyBC+(RES(:,1:g,iNA)-f_0(:,1:g)).*(RES(:,3*g+1:4*g,iNA)-f_0(:,3*g+1:4*g));
    DzBC = DzBC+(RES(:,g+1:2*g,iNA)-f_0(:,g+1:2*g)).*(RES(:,3*g+1:4*g,iNA)-f_0(:,3*g+1:4*g));
    
    DyPW1 = DyPW1+(RES(:,1:g,iNA)-f_0(:,1:g)).*(RES(:,4*g+1:5*g,iNA)-f_0(:,4*g+1:5*g));
    DzPW1 = DzPW1+(RES(:,g+1:2*g,iNA)-f_0(:,g+1:2*g)).*(RES(:,4*g+1:5*g,iNA)-f_0(:,4*g+1:5*g));
    
    DyPW2 = DyPW2+(RES(:,1:g,iNA)-f_0(:,1:g)).*(RES(:,5*g+1:6*g,iNA)-f_0(:,5*g+1:6*g));
    DzPW2 = DzPW2+(RES(:,g+1:2*g,iNA)-f_0(:,g+1:2*g)).*(RES(:,5*g+1:6*g,iNA)-f_0(:,5*g+1:6*g));
    
end;

D = D./N;
DyT = abs(DyT)./N;
DzT = abs(DzT)./N;
DyTtot=abs(D - DzT); %
% DzTtot=D - DyT; %

DyBC = abs(DyBC)./N;
DzBC = abs(DzBC)./N;
DyBCtot=abs(D - DzBC); %
% DzBCtot=D - DyBC; %

DyPW1 = abs(DyPW1)./N;
DzPW1 = abs(DzPW1)./N;
DyPW1tot=abs(D - DzPW1); %
% DzPW1tot=D - DyPW1; %

DyPW2 = abs(DyPW2)./N;
DzPW2 = abs(DzPW2)./N;
DyPW2tot=abs(D-DzPW2); %
% DzPW2tot=D-DyPW2; %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

SyT = DyT./D;
% SzT = DzT./D;
% SyTtot = 1 - SzT; % %
% SzTtot = 1 - SyT; % %
SyTtot = DyTtot./D; %
% SzTtot = DzTtot./D; %

SyBC = DyBC./D;
% SzBC =DzBC./D;
% SyBCtot = 1-SzBC; % %
% SzBCtot = 1-SyBC; % %
SyBCtot = DyBCtot./D; %
% SzBCtot = DzBCtot./D; %

SyPW1 = DyPW1./D;
% SzPW1 = DzPW1./D;
% SyPW1tot = 1 - SzPW1; % %
% SzPW1tot = 1 - SyPW1; % %
SyPW1tot=DyPW1tot./D; %
% SzPW1tot=DzPW1tot./D; %

SyPW2 = DyPW2./D;
% SzPW2 = DzPW2./D;
% SyPW2tot = 1 - SzPW2; % %
% SzPW2tot = 1 - SyPW2; % %
SyPW2tot=DyPW2tot./D; %
% SzPW2tot=DzPW2tot./D; %

fintime=toc;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

path = 'C:\Users\your_mom\Documents\MATLAB\MR FELL FLOW RATE\results\as p\boundaries 1100\len 10 20 30\ ';
path2 = sprintf('%s%s_%s_%s', path, num2str(Nx), num2str(Ny), num2str(N));

if ~exist(path2, 'dir')
	mkdir(path2);
end

% % % % % отрисовка
set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman');
% axis(axis);
h1 = imagesc(SyT);
axis('xy');
c=colorbar;
c.Limits=[0 1];
% title('S_y for T');
saveas(h1, [path2 '\SyT.png'], 'png');
h5= imagesc(SyTtot);
axis('xy');
c=colorbar;
c.Limits=[0 1];
% title('total S_y for T');
saveas(h5, [path2 '\SyTt.png'], 'png');
% h9 = imagesc(SzT);
% colorbar;
% saveas(h9, [path2 '\SzT.png'], 'png');
% h13= imagesc(SzTtot);
% colorbar;
% saveas(h13, [path2 '\SzTt.png'], 'png');
h2 = imagesc(SyBC);
axis('xy');
c=colorbar;
c.Limits=[0 1];
% title('S_y for bc');
saveas(h2, [path2 '\SyBC.png'], 'png');
h6 = imagesc(SyBCtot);
axis('xy');
c=colorbar;
c.Limits=[0 1];
% title('total S_y for bc');
saveas(h6, [path2 '\SyBCt.png'], 'png');
% h10 = imagesc(SzBC);
% colorbar;
% saveas(h10, [path2 '\SzBC.png'], 'png');
% h14 = imagesc(SzBCtot);
% colorbar;
% saveas(h14, [path2 '\SzBCt.png'], 'png');
h3 = imagesc(SyPW1);
axis('xy');
c=colorbar;
c.Limits=[0 1];
% title('S_y for pressure in wells');
saveas(h3, [path2 '\SyPW1.png'], 'png');
h7 = imagesc(SyPW1tot);
axis('xy');
c=colorbar;
c.Limits=[0 1];
% title('total S_y for pressure in w');
saveas(h7, [path2 '\SyPW1t.png'], 'png');
% h11 = imagesc(SzPW1);
% colorbar;
% saveas(h11, [path2 '\SzPW1.png'], 'png');
% h15 = imagesc(SzPW1tot);
% colorbar;
% saveas(h15, [path2 '\SzPW1t.png'], 'png');
h4 = imagesc(SyPW2);
axis('xy');
c=colorbar;
c.Limits=[0 1];
% title('S_y for pressure in w_2');
saveas(h4, [path2 '\SyPW2.png'], 'png');
h8 = imagesc(SyPW2tot);
axis('xy');
c=colorbar;
c.Limits=[0 1];
% title('total S_y for pressure w_2');
saveas(h8, [path2 '\SyPW2t.png'], 'png');
% h12 = imagesc(SzPW2);
% colorbar;
% saveas(h12, [path2 '\SzPW2.png'], 'png');
% h16 = imagesc(SzPW2tot);
% colorbar;
% saveas(h16, [path2 '\SzPW2t.png'], 'png');
h42 = imagesc(sle0);
axis('xy');
c=colorbar;
% c.Limits=[0 1];
% title('solution');
saveas(h42, [path2 '\Solution.png'], 'png');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% запись в файлик
dlmwrite([path2 '\SyT'], SyT, ';');
dlmwrite([path2 '\SyBC'], SyBC, ';');
dlmwrite([path2 '\SyPW1'], SyPW1, ';');
dlmwrite([path2 '\SyPW2'], SyPW2, ';');
% dlmwrite([path2 '\SzT'], SzT, ';');
% dlmwrite([path2 '\SzBC'], SzBC, ';');
% dlmwrite([path2 '\SzPW1'], SzPW1, ';');
% dlmwrite([path2 '\SzPW2'], SzPW2, ';');
dlmwrite([path2 '\SyTt'], SyTtot, ';');
dlmwrite([path2 '\SyBCt'], SyBCtot, ';');
dlmwrite([path2 '\SyPW1t'], SyPW1tot, ';');
dlmwrite([path2 '\SyPW2t'], SyPW2tot, ';');
% dlmwrite([path2 '\SzTt'], SzTtot, ';');
% dlmwrite([path2 '\SzBCt'], SzBCtot, ';');
% dlmwrite([path2 '\SzPW1t'], SzPW1tot, ';');
% dlmwrite([path2 '\SzPW2t'], SzPW2tot, ';');
dlmwrite([path2 '\time'], fintime);