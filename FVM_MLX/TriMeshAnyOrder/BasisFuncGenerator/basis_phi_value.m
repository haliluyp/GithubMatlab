function [value_phi] = basis_phi_value(ac,fac_basis,r)
dof = (r+1)*(r+2)/2; % degree of freedom

B = 1./factorial(0:1:r);  % [1/0!, 1/1!, 1/2!, ... 1/r!]


value_phi = zeros(size(fac_basis,1),dof);
for i = 1:dof
    value_phi(:,i) = B(ac(1,i)+1)*B(ac(2,i)+1)*B(ac(3,i)+1)*...
                     fac_basis(:,ac(1,i)+1,1).*fac_basis(:,ac(2,i)+1,2).*fac_basis(:,ac(3,i)+1,3);
end

end