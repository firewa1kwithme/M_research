function upscaling()
clear function;

    global v k T_x T_y bc WI pw1 pw2
    global Nx Ny hx hy 
       
%     mean_k = 100;
    k = exp(v);
    bc = zeros(1, (Nx-2)*(Ny-2));
    k_x = zeros(Nx-1, Ny);
    k_y = zeros(Nx, Ny-1);

%%%%%%%%%%%%%%%%%%%%%%% boundary
    phi_1 = 1;
    phi_2 = 1;
    phi_3 = 0;
    phi_4 = 0;
%%%%%%%%%%%%%%%%%%%%%%%

    WI = zeros((Nx - 2) * (Ny - 2), 1);
    
    pw1 = 1.5*(1+0.1*randn);
    pw2 = 0.5*(1+0.1*randn);
    
    r_0 = 0.28 * sqrt(hx^2 + hy^2)/2;
    r_w = 0.1;
    h=1;

    for iX = floor(Nx/2)
        for iY = floor(Ny/2)
            l = (iY - 2)*(Nx - 2) + (iX - 1);
            WI(l) = WI(l) + 2*pi*h*k(iX, iY)/(log(r_0/r_w));
        end     
    end
    
    for iX = floor(Nx*3/5)
        for iY = floor(Ny*2/5)
            l = (iY - 2)*(Nx - 2) + (iX - 1);
            WI(l) = WI(l) + 2*pi*h*k(iX, iY)/(log(r_0/r_w));
        end     
    end

    for iX = floor(Nx*3/10)
        for iY = floor(Ny*4/5)
            l = (iY - 2)*(Nx - 2) + (iX - 1);
            WI(l) = WI(l) + 2*pi*h*k(iX, iY)/(log(r_0/r_w));
        end     
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    for iX = 2:Nx
        for iY = 1:Ny    
                
            k_x(iX - 1,iY) = 2*k(iX - 1,iY)*k(iX,iY)/(k(iX - 1,iY) + k(iX,iY));
                
        end
    end

    for iX = 1:Nx
        for iY = 2:Ny
            
            k_y(iX,iY - 1) = 2*k(iX,iY - 1)*k(iX,iY)/(k(iX,iY - 1) + k(iX,iY));
            
        end
    end

    T_x = hy * k_x / hx;
    T_y = hx * k_y / hy;

    for iX = 2:Nx-1
    
        iY = 2;
        l = (Nx - 2)*(iY-2) + (iX-1);
        bc(l) = bc(l)+(phi_1*(Nx-iX)+phi_2*(iX-1))*(1+0.1*randn)/(Nx-1);    
    
        iY = Ny-1;
        l = (Nx - 2)*(iY-2) + (iX-1);
        bc(l) = bc(l)+(phi_4*(Nx-iX)+phi_3*(iX-1))*(1+0.1*randn)/(Nx-1);
    
    end

    for iY = 2:Ny-1
    
        iX = 2;
        l = (Nx - 2)*(iY-2) + (iX-1);
        bc(l) = bc(l)+(phi_1*(Ny-iY)+phi_4*(iY-1))*(1+0.1*randn)/(Ny-1);    
    
        iX = Nx-1;
        l = (Nx - 2)*(iY-2) + (iX-1);
        bc(l) = bc(l)+(phi_2*(Ny-iY)+phi_3*(iY-1))*(1+0.1*randn)/(Ny-1);
    
    end
    
end