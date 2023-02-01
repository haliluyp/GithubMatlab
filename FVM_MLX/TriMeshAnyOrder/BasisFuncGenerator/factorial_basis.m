function [fac_basis] = factorial_basis(gv,r)
%%%  ���ɽ׳˾���: fac_basis = [1,  r*lamda, (r*lamda)*(r*lamda-1), 
%%%    (r*lamda)*(r*lamda-1)*(r*lamda-2), ...,  (r*lamda)*...*(r*lamda-(r-1))]

length_g = size(gv,1);
fac_basis = zeros(length_g,r+1,3);

for k = 1:3
    fac_basis(:,:,k) = [ones(length_g,1), r*gv(:,k)-repmat(0:1:r-1,length_g,1)];
    if r>1
        for i = (r+1):-1:3
            fac_basis(:,i,k) = prod(fac_basis(:,1:i,k),2);
        end
    end
end

end