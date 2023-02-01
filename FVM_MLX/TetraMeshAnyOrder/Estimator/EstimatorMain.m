function [Uh,H1_error,H1_super_error,L2_error,L2_super_error] = EstimatorMain(r,p,t,e,inner,dual_para)
%% 得到数值解相关的误差结果

ver_label  = [1,r*(r+1)*(r+2)/6+1, r*(r+1)*(r+2)/6+r*(r+1)/2+1,(r+1)*(r+2)*(r+3)/6];
P = zeros(3,size(t,2),4);
for i = 1:4
    P(:,:,i) = p(:,t(ver_label(i),:));
end

[Ln,V] = get_mesh_info(P);

%%%%%%%%%%%%%%%%%%%%%% FEMsolver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Uh] = FVMsolver(r,p,t,e,inner,P,Ln,V,dual_para);

U_h =  Uh(t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[w3,g3]=tetra_gausspoints(15);
[value_phi,value_dphi_j] = BasisFunctionMain(g3,r); % 基函数及其导在高斯点g3上取值

X = zeros(size(g3,1),size(P,2),size(P,1));
for i = 1:size(P,1)
    X(:,:,i) = g3(:,1)*P(i,:,1) + g3(:,2)*P(i,:,2) + g3(:,3)*P(i,:,3) + g3(:,4)*P(i,:,4);
end
exact_Ug = zeros(size(X,1),size(X,2),3);
for i = 1:4
    exact_Ug(:,:,i) = exact_functions(X,'u',i);
end

X1(:,:,1) = p(1,:);X1(:,:,2) = p(2,:);X1(:,:,3) = p(3,:);
U = exact_functions(X1,'u',1);
U_I = U(t);

H1_error = H1_estimate_tetra(Ln,V,value_dphi_j,exact_Ug(:,:,2:4),U_h,w3);

H1_super_error = H1_super_estimate_tetra(Ln,V,value_dphi_j,U_h,U_I,w3);

L2_error = L2_estimate_tetra(value_phi,exact_Ug(:,:,1),U_h,V,w3);

L2_super_error = L2_super_estimate_tetra(value_phi,U_h,U_I,V,w3);

end