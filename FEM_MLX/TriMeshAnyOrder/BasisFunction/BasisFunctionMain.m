function [value_phi,value_dphi_j] = BasisFunctionMain(g2,r)
%% r次lagrange计算节点的面积坐标矩阵
%   ac(3,dof)
%   行：三角形的三个面积坐标，和为1
%   列：dof个计算节点（杨辉三角排序）
[ac] = area_coord(r);

%% r次lagrange的阶乘矩阵（辅助生成基函数及其导在高斯点的取值）
% fac_basis(length(g2),r+1,3)
[fac_basis] = factorial_basis(g2,r);

%% r次lagrange基函数在高斯点出取值
% value_phi(length(g2),dof)
[value_phi] = basis_phi_value(ac,fac_basis,r);

%% r次lagrange基函数偏导在高斯点出取值
% value_phi(length(g2),dof，3)
[value_dphi_j] = partial_basis_phi_value(ac,fac_basis,g2,r);
end