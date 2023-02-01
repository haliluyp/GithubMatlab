function [exactF] = exact_functions(x,y,Ftype,i)
% dimension = 2

if Ftype == 'A' % 2*2对称向量的系数A,3个未知量
    if i == 1 % A11
        exactF = 2*exp(x);
    end
    if i == 2 % A12，A21
        exactF = ones(size(x,1),size(x,2));
    end
    if i == 3 % A12，A22
        exactF = exp(y);
    end
end

if Ftype == 'B' % 2*1向量系数B,2个未知量
    if i == 1 % B1
        exactF = x;
    end
    if i == 2 % B2
        exactF = y;
    end
end

if Ftype == 'c' % 标量系数c
    if i == 1
        exactF = sin(x);
    end
end

if Ftype == 'u' % 精确解u,及其偏导
    if i == 1 % u
        exactF = sin(pi*x).*sin(pi*y);
    end
    if i == 2 % Dux
        exactF = pi*cos(pi*x).*sin(pi*y);
    end
    if i == 3 % Duy
        exactF = pi*sin(pi*x).*cos(pi*y);
    end    
end

if Ftype == 'f' % 右端函数f
    if i == 1
        exactF = 2*pi^2*exp(x).*sin(pi*x).*sin(pi*y)-2*pi*exp(x).*cos(pi*x).*sin(pi*y)-pi^2*cos(pi*x).*cos(pi*y) + ...
                 pi^2*exp(y).*sin(pi*x).*sin(pi*y)-pi*exp(y).*sin(pi*x).*cos(pi*y)-pi^2*cos(pi*x).*cos(pi*y) +...
                 pi*x.*cos(pi*x).*sin(pi*y)+ pi*y.*sin(pi*x).*cos(pi*y) + sin(x).*sin(pi*x).*sin(pi*y);
    end
end

end

