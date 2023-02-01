function [F] = get_F(value_phi,V,f,w3)
dof = size(value_phi,2);
num_t = size(V,2);

% ��һ���� aij*dphi*dphj�ڵ�Ԫ�ϻ��� 
F = zeros(dof,num_t);
for i = 1:length(w3)
    for j = 1:dof
        F(j,:) = w3(i)*f(i,:).*value_phi(i,j) + F(j,:);
    end
end

F = F.*repmat(V,dof,1);

F = reshape(F',1,dof*num_t);
end