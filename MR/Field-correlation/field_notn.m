function field_notn(sigma, I, n) 
    clear function;
    global k
    %img and k
%     clear;
    Nx = 20;
    Ny = 20;
    hx = 1;
    hy = 1;
%     n = 100;
    v_0 = zeros(Nx, Ny);
%     I=[10 20 30];
%     sigma = 0.1;
    
    for q = 1:n   
        
        for j = 1:3     
            w(j) = randn;     %для вычисления lambda
        end
            
        ksi_lambda = randn;   %для полей v и v1
        eta_lambda = randn;   %для полей v и v1
        alpha0 = rand;       %для вычисления k
        alpha1 = rand;       %для вычисления k 
        ksi_0 = tan(pi*alpha1/2);	%для вычисления k
        ksi = 1/(1 + ksi_0^2);      %для вычисления k
        eta = alpha0*ksi;  %для вычисления k
     
%вычисление k для вычисления lambda
        while eta >= (ksi_0*ksi)^2
            ksi_0 = tan(pi*rand/2);
            ksi = 1/(1 + ksi_0^2); 
            eta = rand*ksi; 
        end

        k = ksi_0;    %модуль k
%вычисление lambda для скалярного произведения   
        for j = 1:3
            lambda(j) = k*w(j)/(I(j)*(sqrt(w(1)^2 + w(2)^2 + w(3)^2)));%попробовать записать 1 как корень из трех
        end 
    
%задаю x и r
        for iX = 1:Nx            
            for iY = 1:Ny                
                x = [iX*hx iY*hy 0];
                lambda_xr(iX, iY) = dot(lambda, x);
                v(iX, iY) = ksi_lambda*cos(lambda_xr(iX, iY)) + eta_lambda*sin(lambda_xr(iX, iY));
            end            
        end
        
        v_0 = v_0 + v;
        
    end
    
    v = sigma*v_0/sqrt(n);
    meanv = mean(mean(v));
    vfin = meanv + v;
%     k1 = 100*exp(v);
    k = exp(vfin);
%     k3 = n*exp(v);
    
    path = 'C:\Users\your_mom\Documents\MATLAB\MR grafics\figure\FINAL\axis\';
    path2 = sprintf('%s%s_%s_%s', path, num2str(Nx), num2str(Ny), num2str(n));

    if ~exist(path2, 'dir')
        mkdir(path2);
    end
    file_name = strcat('\field_', num2str(sigma), '_', num2str(I));
    set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
    set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman');
    hh = imagesc(k);
    axis('xy');
    c=colorbar;
    c.Limits=[0.5 1.5];
    saveas(hh, [path2 file_name ' k.png'], 'png');
    dlmwrite([path2 file_name ' k'], mean(mean(k)), ';');
%     hh1 = imagesc(k1);
%     colorbar;
%     saveas(hh1, [path2 file_name '_k with 100.png'], 'png');
%     hh2 = imagesc(vfin);
%     colorbar;
%     saveas(hh2, [path2 file_name '_v with mean.png'], 'png');
%     hh3 = imagesc(k2);
%     colorbar;
%     saveas(hh3, [path2 file_name '_k with mean.png'], 'png');
%     hh4 = imagesc(k3);
%     colorbar;
%     saveas(hh4, [path2 file_name 'k with n as dim.png'], 'png');
    
 end
