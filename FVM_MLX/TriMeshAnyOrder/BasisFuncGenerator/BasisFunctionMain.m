function [value_phi,value_dphi_j] = BasisFunctionMain(gv,r)
%% r��lagrange����ڵ������������
%   ac(3,dof)
%   �У������ε�����������꣬��Ϊ1
%   �У�dof������ڵ㣨�����������
[ac] = area_coord(r);

%% r��lagrange�Ľ׳˾��󣨸������ɻ��������䵼�ڸ�˹���ȡֵ��
% fac_basis(length(gv),r+1,3)
[fac_basis] = factorial_basis(gv,r);

%% r��lagrange�������ڸ�˹���ȡֵ
% value_phi(length(gv),dof)
[value_phi] = basis_phi_value(ac,fac_basis,r);

%% r��lagrange������ƫ���ڸ�˹���ȡֵ
% value_phi(length(gv),dof��3)
[value_dphi_j] = partial_basis_phi_value(ac,fac_basis,gv,r);
end