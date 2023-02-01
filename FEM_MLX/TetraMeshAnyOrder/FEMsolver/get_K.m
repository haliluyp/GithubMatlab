function [K] = get_K(value_phi,value_dphi_j,Ln,V,Cof,w3) 
dof = size(value_phi,2);
num_t = size(Ln,2);

% 第一部分 aij*dphi*dphj在单元上积分 
K1 = zeros(dof^2,num_t);
for i = 1:length(w3)
    for j = 1:dof
        value_dphij = [value_dphi_j(i,j,1)*Ln(1,:,1)+value_dphi_j(i,j,2)*Ln(1,:,2)+value_dphi_j(i,j,3)*Ln(1,:,3)+value_dphi_j(i,j,4)*Ln(1,:,4);
                       value_dphi_j(i,j,1)*Ln(2,:,1)+value_dphi_j(i,j,2)*Ln(2,:,2)+value_dphi_j(i,j,3)*Ln(2,:,3)+value_dphi_j(i,j,4)*Ln(2,:,4);
                       value_dphi_j(i,j,1)*Ln(3,:,1)+value_dphi_j(i,j,2)*Ln(3,:,2)+value_dphi_j(i,j,3)*Ln(3,:,3)+value_dphi_j(i,j,4)*Ln(3,:,4)];
        for k = j:dof
            value_dphik = [value_dphi_j(i,k,1)*Ln(1,:,1)+value_dphi_j(i,k,2)*Ln(1,:,2)+value_dphi_j(i,k,3)*Ln(1,:,3)+value_dphi_j(i,k,4)*Ln(1,:,4);
                            value_dphi_j(i,k,1)*Ln(2,:,1)+value_dphi_j(i,k,2)*Ln(2,:,2)+value_dphi_j(i,k,3)*Ln(2,:,3)+value_dphi_j(i,k,4)*Ln(2,:,4);
                            value_dphi_j(i,k,1)*Ln(3,:,1)+value_dphi_j(i,k,2)*Ln(3,:,2)+value_dphi_j(i,k,3)*Ln(3,:,3)+value_dphi_j(i,k,4)*Ln(3,:,4)];
            K1((j-1)*dof+k,:) = w3(i)*sum([Cof(i,:,1).*value_dphij(1,:)+Cof(i,:,2).*value_dphij(2,:)+Cof(i,:,3).*value_dphij(3,:);
                                           Cof(i,:,2).*value_dphij(1,:)+Cof(i,:,4).*value_dphij(2,:)+Cof(i,:,5).*value_dphij(3,:);
                                           Cof(i,:,3).*value_dphij(1,:)+Cof(i,:,5).*value_dphij(2,:)+Cof(i,:,6).*value_dphij(3,:)].*value_dphik,1)...
                                + K1((j-1)*dof+k,:);
        end
    end
end

% 对称性
for j = 2:dof
    for k = 1:(j-1)
        K1((j-1)*dof+k,:) = K1((k-1)*dof+j,:);
    end
end

K1 = K1./repmat(V,dof^2,1);

% 第二部分 bj*dphj*phi在单元上积分 
K2 = zeros(dof^2,num_t);
for i = 1:length(w3)
    for j = 1:dof
        value_dphij = [value_dphi_j(i,j,1)*Ln(1,:,1)+value_dphi_j(i,j,2)*Ln(1,:,2)+value_dphi_j(i,j,3)*Ln(1,:,3)+value_dphi_j(i,j,4)*Ln(1,:,4);
                       value_dphi_j(i,j,1)*Ln(2,:,1)+value_dphi_j(i,j,2)*Ln(2,:,2)+value_dphi_j(i,j,3)*Ln(2,:,3)+value_dphi_j(i,j,4)*Ln(2,:,4);
                       value_dphi_j(i,j,1)*Ln(3,:,1)+value_dphi_j(i,j,2)*Ln(3,:,2)+value_dphi_j(i,j,3)*Ln(3,:,3)+value_dphi_j(i,j,4)*Ln(3,:,4)];
        for k = 1:dof
            K2((j-1)*dof+k,:) = w3(i)*sum([Cof(i,:,7).*value_dphij(1,:);
                                           Cof(i,:,8).*value_dphij(2,:);
                                           Cof(i,:,9).*value_dphij(3,:)]).*value_phi(i,k)...
                                + K2((j-1)*dof+k,:);
        end
    end
end


% 第三部分 c*phi*phj在单元上积分 
K3 = zeros(dof^2,num_t);
for i = 1:length(w3)
    for j = 1:dof
        for k = j:dof
            K3((j-1)*dof+k,:) = w3(i)*Cof(i,:,10).*value_phi(i,j).*value_phi(i,k)...
                                + K3((j-1)*dof+k,:);
        end
    end
end

% 对称性
for j = 2:dof
    for k = 1:(j-1)
        K3((j-1)*dof+k,:) = K3((k-1)*dof+j,:);
    end
end

K3 = K3.*repmat(V,dof^2,1);


K = K1 + K2 + K3;

K = reshape(K',1,dof^2*num_t);
end