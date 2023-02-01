function [value_dphi_j] = partial_basis_phi_value(ac,fac_basis,g2,r)
dof = (r+1)*(r+2)/2; % degree of freedom
length_g = size(g2,1);

A = zeros(length_g,r,3);
C1 = zeros(length_g,r-1,3);
C = zeros(length_g,r+1,3);
for k = 1:3
    A(:,:,k) = r*g2(:,k)-repmat(0:1:r-1,length_g,1); 
    for i = 1:(r-1)
        for j = 1:(i+1)
            flag = 1:(i+1);
            flag(j) = []; 
            C1(:,i,k) = C1(:,i,k) + prod(A(:,flag,k),2); %  prod(A,2), A的所有列相乘
        end
    end
    C(:,:,k) = [repmat([0,r],length_g,1)  r*C1(:,:,k)];
end

B = 1./factorial(0:1:r);  % [1/0!, 1/1!, 1/2!, ... 1/r!]

value_dphi_j = zeros(length_g,dof,3);

for i = 1:dof
    value_dphi_j(:,i,1) = B(ac(1,i)+1)*B(ac(2,i)+1)*B(ac(3,i)+1)*...
        C(:,ac(1,i)+1,1).*fac_basis(:,ac(2,i)+1,2).*fac_basis(:,ac(3,i)+1,3);
    value_dphi_j(:,i,2) = B(ac(1,i)+1)*B(ac(2,i)+1)*B(ac(3,i)+1)*...
        fac_basis(:,ac(1,i)+1,1).*C(:,ac(2,i)+1,2).*fac_basis(:,ac(3,i)+1,3);
    value_dphi_j(:,i,3) = B(ac(1,i)+1)*B(ac(2,i)+1)*B(ac(3,i)+1)*...
        fac_basis(:,ac(1,i)+1,1).*fac_basis(:,ac(2,i)+1,2).*C(:,ac(3,i)+1,3);
end
end