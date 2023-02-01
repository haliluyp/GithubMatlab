function [exactF] = exact_functions(x,y,z,Ftype,i)
% dimension = 3
ModleProblem = 1;

if ModleProblem == 1
if Ftype == 'A' % 3*3对称向量的系数A,6个未知量
    if i == 1 % A11
        exactF = ones(size(x,1),size(x,2));
    end
    if i == 2 % A12，A21
        exactF = zeros(size(x,1),size(x,2));
    end
    if i == 3 % A13，A31
        exactF = zeros(size(x,1),size(x,2));
    end
    if i == 4 % A22
        exactF = ones(size(x,1),size(x,2));
    end
    if i == 5 % A23，A32
        exactF = zeros(size(x,1),size(x,2));
    end
    if i == 6 % A33
        exactF = ones(size(x,1),size(x,2));
    end
end

if Ftype == 'B' % 3*1向量系数B,3个未知量
    if i == 1 % B1
        exactF = zeros(size(x,1),size(x,2));
    end
    if i == 2 % B2
        exactF = zeros(size(x,1),size(x,2));
    end
    if i == 3 % B3
        exactF = zeros(size(x,1),size(x,2));
    end
end

if Ftype == 'c' % 标量系数c
    if i == 1 % c
        exactF = zeros(size(x,1),size(x,2));
    end
end

if Ftype == 'u' % 精确解u,及其偏导
    if i == 1 % u
        exactF = x.*(1-x).*y.*(1-y).*z.*(1-z);
    end
    if i == 2 % Dux
        exactF = (1-2*x).*y.*(1-y).*z.*(1-z);
    end
    if i == 3 % Duy
        exactF = x.*(1-x).*(1-2*y).*z.*(1-z);
    end   
    if i == 4 % Duz
        exactF = x.*(1-x).*y.*(1-y).*(1-2*z);
    end   
end

if Ftype == 'f' % 右端函数f
    if i == 1
        exactF = 2*y.*(1-y).*z.*(1-z)+2*x.*(1-x).*z.*(1-z)+2*y.*(1-y).*x.*(1-x);
    end
end
end

if ModleProblem == 2
if Ftype == 'A' % 3*3对称向量的系数A,6个未知量
    if i == 1 % A11
        exactF = ones(size(x,1),size(x,2));
    end
    if i == 2 % A12，A21
        exactF = zeros(size(x,1),size(x,2));
    end
    if i == 3 % A13，A31
        exactF = zeros(size(x,1),size(x,2));
    end
    if i == 4 % A22
        exactF = 2*ones(size(x,1),size(x,2));
    end
    if i == 5 % A23，A32
        exactF = zeros(size(x,1),size(x,2));
    end
    if i == 6 % A33
        exactF = 3*ones(size(x,1),size(x,2));
    end
end

if Ftype == 'B' % 3*1向量系数B,3个未知量
    if i == 1 % B1
        exactF = 2*ones(size(x,1),size(x,2));
    end
    if i == 2 % B2
        exactF = 3*ones(size(x,1),size(x,2));
    end
    if i == 3 % B3
        exactF = 4*ones(size(x,1),size(x,2));
    end
end

if Ftype == 'c' % 标量系数c
    if i == 1 % c
        exactF = 10*ones(size(x,1),size(x,2));
    end
end

if Ftype == 'u' % 精确解u,及其偏导
    if i == 1 % u
        exactF =  sin(pi*x).*sin(pi*y).*sin(pi*z);
    end
    if i == 2 % Dux
        exactF = pi*cos(pi*x).*sin(pi*y).*sin(pi*z);
    end
    if i == 3 % Duy
        exactF = pi*sin(pi*x).*cos(pi*y).*sin(pi*z);
    end   
    if i == 4 % Duz
        exactF = pi*sin(pi*x).*sin(pi*y).*cos(pi*z);
    end   
end

if Ftype == 'f' % 右端函数f
    if i == 1
        exactF = (6*pi^2+10)*sin(pi*x).*sin(pi*y).*sin(pi*z)+2*pi*cos(pi*x).*sin(pi*y).*sin(pi*z)+...
                  3*pi*sin(pi*x).*cos(pi*y).*sin(pi*z)+4*pi*sin(pi*x).*sin(pi*y).*cos(pi*z);
    end
end
end

if ModleProblem == 3
if Ftype == 'A' % 3*3对称向量的系数A,6个未知量
    if i == 1 % A11
        exactF = exp(x)+3;
    end
    if i == 2 % A12，A21
        exactF = zeros(size(x,1),size(x,2));
    end
    if i == 3 % A13，A31
        exactF = zeros(size(x,1),size(x,2));
    end
    if i == 4 % A22
        exactF = 6+sin(y);
    end
    if i == 5 % A23，A32
        exactF = zeros(size(x,1),size(x,2));
    end
    if i == 6 % A33
        exactF = 9+z.^6;
    end
end

if Ftype == 'B' % 3*1向量系数B,3个未知量
    if i == 1 % B1
        exactF = x;
    end
    if i == 2 % B2
        exactF = y.^2;
    end
    if i == 3 % B3
        exactF = z.^3;
    end
end

if Ftype == 'c' % 标量系数c
    if i == 1 % c
        exactF = cos(x);
    end
end

if Ftype == 'u' % 精确解u,及其偏导
    if i == 1 % u
        exactF = (exp(x)-exp(1)).*(exp(x)-1).*(exp(y)-exp(1)).*(exp(y)-1).*(exp(z)-exp(1)).*(exp(z)-1);
    end
    if i == 2 % Dux
        exactF = (2*exp(2*x)-(exp(1)+1)*exp(x)).*(exp(y)-exp(1)).*(exp(y)-1).*(exp(z)-exp(1)).*(exp(z)-1);
    end
    if i == 3 % Duy
        exactF = (exp(x)-exp(1)).*(exp(x)-1).*(2*exp(2*y)-(exp(1)+1)*exp(y)).*(exp(z)-exp(1)).*(exp(z)-1);
    end   
    if i == 4 % Duz
        exactF = (exp(x)-exp(1)).*(exp(x)-1).*(exp(y)-exp(1)).*(exp(y)-1).*(2*exp(2*z)-(exp(1)+1)*exp(z));
    end   
end

if Ftype == 'f' % 右端函数f
    if i == 1
        exactF = -(exp(x).*(2*exp(2*x)-(exp(1)+1)*exp(x))+(exp(x)+3).*(4*exp(2*x)-(exp(1)+1)*exp(x))).*(exp(y)-exp(1)).*(exp(y)-1).*(exp(z)-exp(1)).*(exp(z)-1)...
                 -(exp(x)-exp(1)).*(exp(x)-1).*(cos(y).*(2*exp(2*y)-(exp(1)+1)*exp(y))+(6+sin(y)).*(4*exp(2*y)-(exp(1)+1)*exp(y))).*(exp(z)-exp(1)).*(exp(z)-1)...
                 -(exp(x)-exp(1)).*(exp(x)-1).*(exp(y)-exp(1)).*(exp(y)-1).*(6*z.^5.*(2*exp(2*z)-(exp(1)+1)*exp(z))+(9+z.^6).*(4*exp(2*z)-(exp(1)+1)*exp(z)))+...
                 x.*(2*exp(2*x)-(exp(1)+1)*exp(x)).*(exp(y)-exp(1)).*(exp(y)-1).*(exp(z)-exp(1)).*(exp(z)-1)+...
                 y.^2.*(exp(x)-exp(1)).*(exp(x)-1).*(2*exp(2*y)-(exp(1)+1)*exp(y)).*(exp(z)-exp(1)).*(exp(z)-1)+...
                 z.^3.*(exp(x)-exp(1)).*(exp(x)-1).*(exp(y)-exp(1)).*(exp(y)-1).*(2*exp(2*z)-(exp(1)+1)*exp(z))+...
                 cos(x).*(exp(x)-exp(1)).*(exp(x)-1).*(exp(y)-exp(1)).*(exp(y)-1).*(exp(z)-exp(1)).*(exp(z)-1);
    end
end
end

end

