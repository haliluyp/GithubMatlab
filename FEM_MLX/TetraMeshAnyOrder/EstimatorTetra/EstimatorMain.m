function [Uh,H1_error,H1_super_error,L2_error,L2_super_error] = EstimatorMain(r,p,t,e,inner)
%% 得到数值解相关的误差结果
[w3,g3]=tetra_gauss_points(15);
[value_phi,value_dphi_j] = BasisFunctionMain(g3,r); % 基函数及其导在高斯点g2上取值

vert2 = r*(r+1)*(r+2)/6+1; vert3 = r*(r+1)*(r+2)/6+r*(r+1)/2+1; vert4 = (r+1)*(r+2)*(r+3)/6;
P1 = p(:,t(1,:));P2 = p(:,t(vert2,:));P3 = p(:,t(vert3,:));P4 = p(:,t(vert4,:));
[Ln,V] = get_mesh_info(P1,P2,P3,P4);

%%%%%%%%%%%%%%%%%%%%%% FEMsolver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Uh,exact_Ug] = FEMsolver(w3,g3,p,t,e,inner,P1,P2,P3,P4,value_phi,value_dphi_j,Ln,V);

U_h =  Uh(t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = exact_functions(p(1,:),p(2,:),p(3,:),'u',1);
U_I = U(t);

H1_error = H1_estimate_tetra(Ln,V,value_dphi_j,exact_Ug(:,:,2:4),U_h,w3);

H1_super_error = H1_super_estimate_tetra(Ln,V,value_dphi_j,U_h,U_I,w3);

L2_error = L2_estimate_tetra(value_phi,exact_Ug(:,:,1),U_h,V,w3);

L2_super_error = L2_super_estimate_tetra(value_phi,U_h,U_I,V,w3);

end