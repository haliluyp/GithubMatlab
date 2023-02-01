function [F] = get_F(value_phi,S,f,w2)
dof = size(value_phi,2);
num_t = size(S,2);

% ��һ���� aij*dphi*dphj�ڵ�Ԫ�ϻ��� 
F = zeros(dof,num_t);
for i = 1:length(w2)
    for j = 1:dof
        F(j,:) = w2(i)*f(i,:).*value_phi(i,j) + F(j,:);
    end
end

F = F.*repmat(S,dof,1);

F = reshape(F',1,dof*num_t);

end