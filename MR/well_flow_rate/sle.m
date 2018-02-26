function pp = sle(T_x, T_y, bc, pw1, pw2)
clear function;
    
    %global T_x T_y bc pw1 pw2
    global Nx Ny WI
   
    A = zeros((Nx - 2) * (Ny - 2));
    b = zeros((Nx - 2) * (Ny - 2), 1);
    q = zeros((Nx - 2) * (Ny - 2), 1);
%%%%%%%%%%%%%   matrix A
    for iX = 2:Nx-1
        for iY = 2:Ny-1
            
            l = (iY - 2)*(Nx - 2) + (iX - 1);
            A(l,l) =-T_x(iX-1,iY)-T_x(iX,iY)-T_y(iX,iY)-T_y(iX,iY-1)-WI(l);
            
        end
    end
    
     for iX = 3:Nx - 1
        for iY = 2:Ny - 1
            
            l = (iY - 2)*(Nx - 2) + (iX - 1);
            A(l-1,l) = T_x(iX-1,iY);
            
        end
    end

    for iX = 2:Nx - 2
        for iY = 2:Ny - 1
            
            l = (iY - 2)*(Nx - 2) + (iX - 1);
            A(l+1,l) = T_x(iX,iY);
            
        end
    end

    for iX = 2:Nx - 1
        for iY = 3:Ny - 1
            
            l = (iY - 2)*(Nx - 2) + (iX - 1);
            A(l-Nx+2,l) = T_y(iX,iY-1);
%             A(l,l-1) = T_y(iX,iY-1);
            
        end
    end

    for iX = 2:Nx-1
        for iY = 2:Ny-2
            
            l = (iY - 2)*(Nx - 2) + (iX - 1);
            A(l+Nx-2,l) = T_y(iX,iY);
%             A(l,l+1) = T_y(iX,iY);
            
        end
    end

%%%%%%%%%%%%%%%%%%%%%% b
    for iX = 2:Nx-1
    
        iY=2;
        l = (Nx - 2)*(iY-2) + (iX-1);
        b(l)= -T_y(iX,iY-1)*bc(l);    
    
        iY=Ny-1;
        l = (Nx - 2)*(iY-2) + (iX-1);
        b(l)= -T_y(iX,iY)*bc(l);
    
    end

    for iY = 2:Ny-1
    
        iX = 2;
        l = (Nx - 2)*(iY-2) + (iX-1);
        b(l)=  -T_x(iX-1,iY)*bc(l);    
    
        iX = Nx-1;
        l = (Nx - 2)*(iY-2) + (iX-1);   
        b(l)=  -T_x(iX,iY)*bc(l);
    
    end
    
    % b либо сделать pw как вектор, либо как WI задать отдельно не 0 b
    % то есть нигде нет wi потом отдельно прибать значения в индексах wi
    
    for iX = floor(Nx/2)
        for iY = floor(Ny/2)
            l = (iY - 2)*(Nx - 2) + (iX - 1);
            b(l) = b(l) - WI(l)*pw1;
        end     
    end
    
    for iX = floor(Nx*3/5)
        for iY = floor(Ny*2/5)
            l = (iY - 2)*(Nx - 2) + (iX - 1);
            b(l) = b(l) - WI(l)*pw2;
        end     
    end

    for iX = floor(Nx*3/10)
        for iY = floor(Ny*4/5)
            l = (iY - 2)*(Nx - 2) + (iX - 1);
            b(l) = b(l) - WI(l)*pw2;
        end     
    end    
    
    A_sp=sparse(A);
    p = A_sp\b;
    
    pp = reshape(p,Nx-2,Ny-2);
    
%    % думаю, стоит pw задать как вектор.     
%     for iX = floor(Nx*4/5)
%         for iY = floor(Ny*4/5)
%             l = (iY - 2)*(Nx - 2) + (iX - 1);
%             q(l) = q(l) + WI(l)*(p(l)-pw1);
%         end     
%     end 
%     for iX = floor(Nx*3/5)
%         for iY = floor(Ny*2/5)
%             l = (iY - 2)*(Nx - 2) + (iX - 1);
%             q(l) = q(l)+WI(l)*(p(l)-pw2);
%         end     
%     end
%     
%     for iX = floor(Nx*3/10)
%         for iY = floor(Ny*4/5)
%             l = (iY - 2)*(Nx - 2) + (iX - 1);
%             q(l) = q(l)+WI(l)*(p(l)-pw2);
%         end     
%     end
%     
%     qq = reshape(p,Nx-2,Ny-2);
%     %imagesc(pp);
end