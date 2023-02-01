function [value_dphi_j] = partial_basis_phi_value(ac,fac_basis,g3,r)
dof = (r+1)*(r+2)*(r+3)/6; % degree of freedom
length_g = size(g3,1);

A = zeros(length_g,r,4);
C1 = zeros(length_g,r-1,4);
C = zeros(length_g,r+1,4);
for k = 1:4
    A(:,:,k) = r*g3(:,k)-repmat(0:1:r-1,length_g,1); 
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

value_dphi_j = zeros(length_g,dof,4);
for i = 1:dof
    value_dphi_j(:,i,1) = B(ac(1,i)+1)*B(ac(2,i)+1)*B(ac(3,i)+1)*B(ac(4,i)+1)*...
        C(:,ac(1,i)+1,1).*fac_basis(:,ac(2,i)+1,2).*fac_basis(:,ac(3,i)+1,3).*fac_basis(:,ac(4,i)+1,4);
    value_dphi_j(:,i,2) = B(ac(1,i)+1)*B(ac(2,i)+1)*B(ac(3,i)+1)*B(ac(4,i)+1)*...
        fac_basis(:,ac(1,i)+1,1).*C(:,ac(2,i)+1,2).*fac_basis(:,ac(3,i)+1,3).*fac_basis(:,ac(4,i)+1,4);
    value_dphi_j(:,i,3) = B(ac(1,i)+1)*B(ac(2,i)+1)*B(ac(3,i)+1)*B(ac(4,i)+1)*...
        fac_basis(:,ac(1,i)+1,1).*fac_basis(:,ac(2,i)+1,2).*C(:,ac(3,i)+1,3).*fac_basis(:,ac(4,i)+1,4);
    value_dphi_j(:,i,4) = B(ac(1,i)+1)*B(ac(2,i)+1)*B(ac(3,i)+1)*B(ac(4,i)+1)*...
        fac_basis(:,ac(1,i)+1,1).*fac_basis(:,ac(2,i)+1,2).*fac_basis(:,ac(3,i)+1,3).*C(:,ac(4,i)+1,4);
end

end