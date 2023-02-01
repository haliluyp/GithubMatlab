function [Uh] = FVMsolver(r,p,t,e,inner,P,Ln,S,dual_para)
dof = size(t,1);
num_t = size(t,2);
num_p = size(p,2);
K = zeros(dof^2,num_t);
F = zeros(dof,num_t);

[ac] = area_coord(r);
[ReLoc] = relocatenodes(ac,r); %������������������ĶԳƷ�ʽ�����������ڵ�

for d2 = 1:3 % �����ε�ԪK���ڱ������ķ�Ϊ'3'���Ӳ���
    for d1 = 1:2 % �Ӳ��� ����K �ı߷�Ϊ '2'��С������
        [Dual] = dual_info(r,ReLoc,d2,d1,dual_para);
        [K0] = get_K(r,P,Ln,S,Dual);
        K = K + K0;
        [F0] = get_F(r,P,Dual);
        F = F + F0;
    end
end

K = reshape(K',1,dof^2*num_t);
F = reshape(F',1,dof*num_t);

[K] = assemble_K(K,t,num_p);
[F] = assemble_F(F,t,num_p);

X(:,:,1) = p(1,e);X(:,:,2) = p(2,e);
Ue = exact_functions(X,'u',1);
[K,F] = with_boundary_condition(K,F,e,inner,Ue); 
Uh = zeros(num_p,1);
Uh(inner,1) = K\F;
Uh(e,1) = Ue;
%%%%%%%
end