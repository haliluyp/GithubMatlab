function [value_phi,value_dphi_j] = BasisFunctionMain(g3,r)
%% r��lagrange����ڵ������������
%   ac(4,dof)
%   �У���������ĸ�������꣬��Ϊ1
%   �У�dof������ڵ㣨�����������
[ac] = area_coord(r);

%% r��lagrange�Ľ׳˾��󣨸������ɻ��������䵼�ڸ�˹���ȡֵ��
% fac_basis(length(g2),r+1,4)
[fac_basis] = factorial_basis(g3,r);

%% r��lagrange�������ڸ�˹�㴦ȡֵ
% value_phi(length(g2),dof)
[value_phi] = basis_phi_value(ac,fac_basis,r);

%% r��lagrange������ƫ���ڸ�˹�㴦ȡֵ
% value_phi(length(g2),dof��4)
[value_dphi_j] = partial_basis_phi_value(ac,fac_basis,g3,r);
end