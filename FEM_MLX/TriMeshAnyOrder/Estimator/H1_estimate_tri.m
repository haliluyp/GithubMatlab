function  [H1_error] = H1_estimate_tri(Ln,S,value_dphi_j,exact_dUg,Uh,w2)
dof = size(value_dphi_j,2);

H1_error=0;
Ln(:,:,1) = Ln(:,:,1)./repmat(S,2,1);
Ln(:,:,2) = Ln(:,:,2)./repmat(S,2,1);
Ln(:,:,3) = Ln(:,:,3)./repmat(S,2,1);
for i=1:length(w2)
    value_dphi_x = (value_dphi_j(i,:,1))'.*repmat(Ln(1,:,1),dof,1)+(value_dphi_j(i,:,2))'.*repmat(Ln(1,:,2),dof,1)+(value_dphi_j(i,:,3))'.*repmat(Ln(1,:,3),dof,1);
    value_dphi_y = (value_dphi_j(i,:,1))'.*repmat(Ln(2,:,1),dof,1)+(value_dphi_j(i,:,2))'.*repmat(Ln(2,:,2),dof,1)+(value_dphi_j(i,:,3))'.*repmat(Ln(2,:,3),dof,1);
    
    H1_error=H1_error+w2(i)*((exact_dUg(i,:,1)-sum(value_dphi_x.*Uh)).^2+(exact_dUg(i,:,2)-sum(value_dphi_y.*Uh)).^2);
end
H1_error=sqrt(sum(S.*H1_error));
end