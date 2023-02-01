function [value_phi,value_dphi_j] = BasisFunctionMain(g3,r)
%% r次lagrange计算节点的体积坐标矩阵
%   ac(4,dof)
%   行：四面体的四个体积坐标，和为1
%   列：dof个计算节点（杨辉三角排序）
[ac] = area_coord(r);

%% r次lagrange的阶乘矩阵（辅助生成基函数及其导在高斯点的取值）
% fac_basis(length(g2),r+1,4)
[fac_basis] = factorial_basis(g3,r);

%% r次lagrange基函数在高斯点处取值
% value_phi(length(g2),dof)
[value_phi] = basis_phi_value(ac,fac_basis,r);

%% r次lagrange基函数偏导在高斯点处取值
% value_phi(length(g2),dof，4)
[value_dphi_j] = partial_basis_phi_value(ac,fac_basis,g3,r);
end