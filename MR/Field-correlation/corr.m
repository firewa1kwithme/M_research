 function corr(sigma, I, R)
%grafics
    clear function;
%     R = 40000; %100000, ���� ������ ����� - 40000
    N = 50;
    n = 50; %��������� �� 50 ���� 10
    %I = [1 1 1];         %����� ����������
    cor = zeros(1, N);

    for l = 1:R
    
        v_0 = zeros(1, N);
    
        for q = 1:n
        
            for j = 1:3
            
                w(j) = randn;     %��� ���������� lambda
            end
    
            ksi_lambda = randn;   %��� ����� v � v1
            eta_lambda = randn;   %��� ����� v � v1
            ksi_0 = tan(pi*rand/2);	%��� ���������� k
            ksi = 1/(1 + ksi_0^2);      %��� ���������� k
            eta = rand*ksi;  %��� ���������� k
     
%���������� k ��� ���������� lambda
            while eta >= (ksi_0*ksi)^2
                ksi_0 = tan(pi*rand/2);
                ksi = 1/(1 + ksi_0^2); 
                eta = rand*ksi; 
            end

            k = ksi_0;    %������ k
    
%���������� lambda ��� ���������� ������������   
            for j = 1:3
                
                lambda(j) = k*w(j)/(I(j)*(sqrt(w(1)^2 + w(2)^2 + w(3)^2)));%����������� �������� 1 ��� ������ �� ����
 
            end 
    
%����� x � r
            for i = 1:N
                x(i,1) = 0 + 5*(i-1)/(N);
                x(i,2) = 0 + 10*(i-1)/(N);
                x(i,3) = 0 + 15*(i-1)/(N);
                r(i) = sqrt((x(i,1)/I(1))^2 + (x(i,2)/I(2))^2 + (x(i,3)/I(3))^2);
    
%���������� ���������� ������������  
                lambda_xr(i) = dot(lambda, x(i,:));
                v(i) = sigma*(ksi_lambda*cos(lambda_xr(i)) + eta_lambda*sin(lambda_xr(i)));

            end
            v_0 = v_0 + v;

        end

        v = v_0/sqrt(n);
%�������������� �������
        for i = 1:N
            cor(i) = cor(i) + v(1)*v(i);
        end
    end
    
        path = 'C:\Users\your_mom\Documents\MATLAB\MR grafics\figure\testofQuite\';
    if ~exist(path, 'dir')
        mkdir(path);
    end
    
    set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
    set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman');
    dr = sqrt((5/N)^2 + (10/N)^2 + (15/N)^2);
    cor_exact = sigma^2*exp(-r);
    x1 = 0:dr:(N-1)*dr;
    hh = figure;
    plot(x1, (cor/R),'r', x1, cor_exact, 'g', 'LineWidth', 3);
    xlabel('�');
    ylabel('correlation function');
    legend('theoretical correlation function','exact correlation function');
    grid on;
    axis on;
    file_name = strcat('corr_', num2str(sigma), '_', num2str(I), '.png');
    saveas(hh, [path file_name], 'png');  
    
 end

