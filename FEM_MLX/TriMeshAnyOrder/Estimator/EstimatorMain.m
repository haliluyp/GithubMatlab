function [Uh,H1_error,H1_super_error,L2_error,L2_super_error] = EstimatorMain(r,p,t,e,inner,ver_label)
%% 得到数值解相关的误差结果
[w2,g2] = tri_gauss_points(13);
[value_phi,value_dphi_j] = BasisFunctionMain(g2,r); % 基函数及其导在高斯点g2上取值

P1 = p(:,t(ver_label(1),:));P2 = p(:,t(ver_label(2),:));P3 = p(:,t(ver_label(3),:));
[Ln,S] = get_mesh_info(P1,P2,P3);

%%%%%%%%%%%%%%%%%%%%%% FEMsolver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Uh,exact_Ug] = FEMsolver(w2,g2,p,t,e,inner,P1,P2,P3,value_phi,value_dphi_j,Ln,S);

U_h =  Uh(t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = exact_functions(p(1,:),p(2,:),'u',1);
U_I = U(t);
H1_error = H1_estimate_tri(Ln,S,value_dphi_j,exact_Ug(:,:,2:3),U_h,w2);

H1_super_error = H1_super_estimate_tri(Ln,S,value_dphi_j,U_h,U_I,w2);

L2_error = L2_estimate_tri(value_phi,exact_Ug(:,:,1),U_h,S,w2);

L2_super_error = L2_super_estimate_tri(value_phi,U_h,U_I,S,w2);

end