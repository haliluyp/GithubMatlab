function [fac_basis] = factorial_basis(g3,r)
%%%  生成阶乘矩阵: fac_basisi = [1,  r*Li, (r*Li)*(r*Li-1), 
%%%                 (r*Li)*(r*Li-1)*(r*Li-2), ...,  (r*Li)*...*(r*Li-(r-1))]

length_g = size(g3,1);
fac_basis = zeros(length_g,r+1,4);

for k = 1:4
    fac_basis(:,:,k) = [ones(length_g,1), r*g3(:,k)-repmat(0:1:r-1,length_g,1)];
    if r>1
        for i = (r+1):-1:3
            fac_basis(:,i,k) = prod(fac_basis(:,1:i,k),2); % fac_basis的1:i列连乘
        end
    end
end

end