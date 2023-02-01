function [Uh,H1_error,H1_super_error,L2_error,L2_super_error] = EstimatorMain(r,p,t,e,inner,dual_para)
% 得到数值解相关的误差结果

ver_label =  [1, r*(r+1)/2+1, (r+1)*(r+2)/2]; % 三个顶点的序号
P = zeros(2,size(t,2),3);
for i = 1:3
    P(:,:,i) = p(:,t(ver_label(i),:));
end
[Ln,S] = get_mesh_info(P);

%%%%%%%%%%%%%%%%%%%%%% FVMsolver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Uh] = FVMsolver(r,p,t,e,inner,P,Ln,S,dual_para);

U_h =  Uh(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[w2,g2] = tri_gausspoints(13);
[value_phi,value_dphi_j] = BasisFunctionMain(g2,r);

X = zeros(size(g2,1),size(P,2),size(P,1));
for i = 1:size(P,1)
    X(:,:,i) = g2(:,1)*P(i,:,1) + g2(:,2)*P(i,:,2) + g2(:,3)*P(i,:,3);
end
exact_Ug = zeros(size(X,1),size(X,2),3);
for i = 1:3
    exact_Ug(:,:,i) = exact_functions(X,'u',i);
end

X1(:,:,1) = p(1,:);X1(:,:,2) = p(2,:);
U = exact_functions(X1,'u',1);
U_I = U(t);
H1_error = H1_estimate_tri(Ln,S,value_dphi_j,exact_Ug(:,:,2:3),U_h,w2);

H1_super_error = H1_super_estimate_tri(Ln,S,value_dphi_j,U_h,U_I,w2);

L2_error = L2_estimate_tri(value_phi,exact_Ug(:,:,1),U_h,S,w2);

L2_super_error = L2_super_estimate_tri(value_phi,U_h,U_I,S,w2);

end