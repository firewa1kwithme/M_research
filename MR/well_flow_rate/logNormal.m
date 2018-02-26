function logNormal()
clear function;
 
    global v 
    global Nx Ny hx hy

    sigma = 0.1;
    n = 100;
    I = [10 20 30];         %длины корреляции
    v = zeros(Nx, Ny);
    v_0 = zeros(Nx, Ny);
    lambda = [0 0 0];
    
    for q = 1:n    
        
        w = [randn randn randn];
            
        ksi_lambda = randn;   %для полей v и v1
        eta_lambda = randn;   %для полей v и v1
        ksi_0 = tan(pi*rand/2);	%для вычисления k
        ksi = 1/(1 + ksi_0^2);      %для вычисления k
        eta = rand*ksi;  %для вычисления k
     
%вычисление k для вычисления lambda
        while eta >= (ksi_0*ksi)^2
            ksi_0 = tan(pi*rand/2);
            ksi = 1/(1 + ksi_0^2); 
            eta = rand*ksi; 
        end

        k = ksi_0;    %модуль k
%вычисление lambda для скалярного произведения   
        for j = 1:3
            lambda(j) = k*w(j)/(I(j)*(sqrt(w(1)^2 + w(2)^2 + w(3)^2)));
        end 
    
%задаю x и r
        for iX = 1:Nx
            
            for iY = 1:Ny
                x = [iX*hx iY*hy 0];
                lambda_xr(iX, iY) = dot(lambda, x);
                v(iX, iY) = sigma*(ksi_lambda*cos(lambda_xr(iX, iY)) + eta_lambda*sin(lambda_xr(iX, iY)));
            end
            
        end
        
        v_0 = v_0 + v;
    end
    
    vfin = v_0/sqrt(n);
    meanv = mean(mean(v));
    v = meanv+vfin;
    
end