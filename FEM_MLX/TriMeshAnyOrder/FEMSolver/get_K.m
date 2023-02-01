function [K] = get_K(value_phi,value_dphi_j,Ln,S,Cof,w2)  
dof = size(value_phi,2);
num_t = size(Ln,2);

% ��һ���� aij*dphi*dphj�ڵ�Ԫ�ϻ��� 
K1 = zeros(dof^2,num_t);
for i = 1:length(w2)
    for j = 1:dof
        value_dphij = [(value_dphi_j(i,j,1)*Ln(1,:,1)+value_dphi_j(i,j,2)*Ln(1,:,2)+value_dphi_j(i,j,3)*Ln(1,:,3));
                       (value_dphi_j(i,j,1)*Ln(2,:,1)+value_dphi_j(i,j,2)*Ln(2,:,2)+value_dphi_j(i,j,3)*Ln(2,:,3))];
        for k = j:dof
            value_dphik = [(value_dphi_j(i,k,1)*Ln(1,:,1)+value_dphi_j(i,k,2)*Ln(1,:,2)+value_dphi_j(i,k,3)*Ln(1,:,3));
                           (value_dphi_j(i,k,1)*Ln(2,:,1)+value_dphi_j(i,k,2)*Ln(2,:,2)+value_dphi_j(i,k,3)*Ln(2,:,3))];
            K1((j-1)*dof+k,:) = w2(i)*sum([Cof(i,:,1).*value_dphij(1,:)+Cof(i,:,2).*value_dphij(2,:);
                                           Cof(i,:,2).*value_dphij(1,:)+Cof(i,:,3).*value_dphij(2,:)].*value_dphik,1)...
                                + K1((j-1)*dof+k,:);
        end
    end
end

% �Գ���
for j = 2:dof
    for k = 1:(j-1)
        K1((j-1)*dof+k,:) = K1((k-1)*dof+j,:);
    end
end

K1 = K1./repmat(S,dof^2,1);

% �ڶ����� bj*dphj*phi�ڵ�Ԫ�ϻ��� 
K2 = zeros(dof^2,num_t);
for i = 1:length(w2)
    for j = 1:dof
        value_dphij = [(value_dphi_j(i,j,1)*Ln(1,:,1)+value_dphi_j(i,j,2)*Ln(1,:,2)+value_dphi_j(i,j,3)*Ln(1,:,3));
                        (value_dphi_j(i,j,1)*Ln(2,:,1)+value_dphi_j(i,j,2)*Ln(2,:,2)+value_dphi_j(i,j,3)*Ln(2,:,3))];
        for k = 1:dof
            K2((j-1)*dof+k,:) = w2(i)*sum([Cof(i,:,4).*value_dphij(1,:);
                                           Cof(i,:,5).*value_dphij(2,:)]).*value_phi(i,k)...
                                + K2((j-1)*dof+k,:);
        end
    end
end

% �������� c*phi*phj�ڵ�Ԫ�ϻ��� 
K3 = zeros(dof^2,num_t);
for i = 1:length(w2)
    for j = 1:dof
        for k = j:dof
            K3((j-1)*dof+k,:) = w2(i)*Cof(i,:,6).*value_phi(i,j).*value_phi(i,k)...
                                + K3((j-1)*dof+k,:);
        end
    end
end

% �Գ���
for j = 2:dof
    for k = 1:(j-1)
        K3((j-1)*dof+k,:) = K3((k-1)*dof+j,:);
    end
end

K3 = K3.*repmat(S,dof^2,1);

K = K1 + K2 + K3;

K = reshape(K',1,dof^2*num_t);
end