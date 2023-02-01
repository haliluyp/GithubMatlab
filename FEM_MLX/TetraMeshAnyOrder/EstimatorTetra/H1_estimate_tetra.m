function  [H1_error] = H1_estimate_tetra(Ln,V,value_dphi_j,exact_Ug,Uh,w3)
dof = size(value_dphi_j,2);

H1_error=0;
for i = 1:4
    Ln(:,:,i) = Ln(:,:,i)./repmat(V,3,1);
end
for i=1:length(w3)
    value_dphi_x = (value_dphi_j(i,:,1))'.*repmat(Ln(1,:,1),dof,1)+(value_dphi_j(i,:,2))'.*repmat(Ln(1,:,2),dof,1)+(value_dphi_j(i,:,3))'.*repmat(Ln(1,:,3),dof,1)+(value_dphi_j(i,:,4))'.*repmat(Ln(1,:,4),dof,1);
    value_dphi_y = (value_dphi_j(i,:,1))'.*repmat(Ln(2,:,1),dof,1)+(value_dphi_j(i,:,2))'.*repmat(Ln(2,:,2),dof,1)+(value_dphi_j(i,:,3))'.*repmat(Ln(2,:,3),dof,1)+(value_dphi_j(i,:,4))'.*repmat(Ln(2,:,4),dof,1);
    value_dphi_z = (value_dphi_j(i,:,1))'.*repmat(Ln(3,:,1),dof,1)+(value_dphi_j(i,:,2))'.*repmat(Ln(3,:,2),dof,1)+(value_dphi_j(i,:,3))'.*repmat(Ln(3,:,3),dof,1)+(value_dphi_j(i,:,4))'.*repmat(Ln(3,:,4),dof,1);
    
    H1_error=H1_error+w3(i)*((exact_Ug(i,:,1)-sum(value_dphi_x.*Uh)).^2+(exact_Ug(i,:,2)-sum(value_dphi_y.*Uh)).^2+(exact_Ug(i,:,3)-sum(value_dphi_z.*Uh)).^2);
end
H1_error=sqrt(sum(V.*H1_error));
end