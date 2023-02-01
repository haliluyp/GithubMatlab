function [K] = get_K(r,P,Ln,V,Dual)
dof = (r+1)*(r+2)*(r+3)/6;
num_t = size(Ln,2);

[w2,g2] = tri_gausspoints(13);
[w3,g3] = tetra_gausspoints(15);

K1 = zeros(dof^2,num_t);
Ns = zeros(1,3);
CofL = zeros(size(g2,1),size(P,2),6);
for s = 1:size(Dual.S,1)
    Tmatrix(1,:) = [1 Dual.N(Dual.S(s,1),1:3)];
    Tmatrix(2,:) = [1 Dual.N(Dual.S(s,2),1:3)];
    Tmatrix(3,:) = [1 Dual.N(Dual.S(s,3),1:3)];
    for i = 1:3
        Tm = Tmatrix; Tm(:,i+1) = [];
        Ns(i) = (-1)^(i+1)*det(Tm);
    end
    Dual_Sn = Ns(1)*Ln(:,:,1) + Ns(2)*Ln(:,:,2) + Ns(3)*Ln(:,:,3);
    gv2 = g2(:,1).*Dual.N(Dual.S(s,1),:) + g2(:,2).*Dual.N(Dual.S(s,2),:) + g2(:,3).*Dual.N(Dual.S(s,3),:);
    [~,v_dphi_j] = BasisFunctionMain(gv2,r);
    X1(:,:,1) = gv2(:,1)*P(1,:,1) + gv2(:,2)*P(1,:,2) + gv2(:,3)*P(1,:,3) + gv2(:,4)*P(1,:,4);
    X1(:,:,2) = gv2(:,1)*P(2,:,1) + gv2(:,2)*P(2,:,2) + gv2(:,3)*P(2,:,3) + gv2(:,4)*P(2,:,4);
    X1(:,:,3) = gv2(:,1)*P(3,:,1) + gv2(:,2)*P(3,:,2) + gv2(:,3)*P(3,:,3) + gv2(:,4)*P(3,:,4);
    for i = 1:6
        CofL(:,:,i) = exact_functions(X1,'A',i);
    end
    Sd = Dual.Sd{1,s};
    for j = 1:dof
        value_dphij1 = v_dphi_j(:,j,1).*Ln(1,:,1)+v_dphi_j(:,j,2).*Ln(1,:,2)+v_dphi_j(:,j,3).*Ln(1,:,3)+v_dphi_j(:,j,4).*Ln(1,:,4);
        value_dphij2 = v_dphi_j(:,j,1).*Ln(2,:,1)+v_dphi_j(:,j,2).*Ln(2,:,2)+v_dphi_j(:,j,3).*Ln(2,:,3)+v_dphi_j(:,j,4).*Ln(2,:,4);
        value_dphij3 = v_dphi_j(:,j,1).*Ln(3,:,1)+v_dphi_j(:,j,2).*Ln(3,:,2)+v_dphi_j(:,j,3).*Ln(3,:,3)+v_dphi_j(:,j,4).*Ln(3,:,4);
        Int_S = sum(w2.*((CofL(:,:,1).*value_dphij1 + CofL(:,:,2).*value_dphij2 + CofL(:,:,3).*value_dphij3).*repmat(Dual_Sn(1,:),length(w2),1) + ...
                         (CofL(:,:,2).*value_dphij1 + CofL(:,:,4).*value_dphij2 + CofL(:,:,5).*value_dphij3).*repmat(Dual_Sn(2,:),length(w2),1) + ...
                         (CofL(:,:,3).*value_dphij1 + CofL(:,:,5).*value_dphij2 + CofL(:,:,6).*value_dphij3).*repmat(Dual_Sn(3,:),length(w2),1)),1);
        for m = 1:size(Sd,2)
            K1((j-1)*dof+Sd(1,m),:) = K1((j-1)*dof+Sd(1,m),:) + Sd(2,m)*Int_S;
        end
    end
end
K1 = K1./repmat(V,dof^2,1);

v4 = zeros(size(P,1),size(P,2),3);
K2 = zeros(dof^2,num_t);
K3 = zeros(dof^2,num_t);
for v = 1:size(Dual.V,1)
    for k1 = 1:4
            v4(:,:,k1) = Dual.N(Dual.V(v,k1),1)*P(:,:,1) + Dual.N(Dual.V(v,k1),2)*P(:,:,2) + ...
                         Dual.N(Dual.V(v,k1),3)*P(:,:,3) + Dual.N(Dual.V(v,k1),4)*P(:,:,4);
    end
    V_dual(1,:) = tetra_volume(v4);% 
    
    X2(:,:,1) = g3(:,1).*v4(1,:,1) + g3(:,2).*v4(1,:,2) + g3(:,3).*v4(1,:,3) + g3(:,4).*v4(1,:,4); 
    X2(:,:,2) = g3(:,1).*v4(2,:,1) + g3(:,2).*v4(2,:,2) + g3(:,3).*v4(2,:,3) + g3(:,4).*v4(2,:,4); 
    X2(:,:,3) = g3(:,1).*v4(3,:,1) + g3(:,2).*v4(3,:,2) + g3(:,3).*v4(3,:,3) + g3(:,4).*v4(3,:,4); 
    
    CofV(:,:,1) = exact_functions(X2,'B',1);
    CofV(:,:,2) = exact_functions(X2,'B',2);
    CofV(:,:,3) = exact_functions(X2,'B',3);
    CofV(:,:,4) = exact_functions(X2,'c',1);
    gv3 = g3(:,1).*Dual.N(Dual.V(v,1),:) + g3(:,2).*Dual.N(Dual.V(v,2),:) + g3(:,3).*Dual.N(Dual.V(v,3),:) + g3(:,4).*Dual.N(Dual.V(v,4),:);
    [v_phi,v_dphi_j] = BasisFunctionMain(gv3,r);
    Vd = Dual.Vd{1,v};
    for j = 1:dof
        value_dphij1 = v_dphi_j(:,j,1).*Ln(1,:,1)+v_dphi_j(:,j,2).*Ln(1,:,2)+v_dphi_j(:,j,3).*Ln(1,:,3)+v_dphi_j(:,j,4).*Ln(1,:,4);
        value_dphij2 = v_dphi_j(:,j,1).*Ln(2,:,1)+v_dphi_j(:,j,2).*Ln(2,:,2)+v_dphi_j(:,j,3).*Ln(2,:,3)+v_dphi_j(:,j,4).*Ln(2,:,4);
        value_dphij3 = v_dphi_j(:,j,1).*Ln(3,:,1)+v_dphi_j(:,j,2).*Ln(3,:,2)+v_dphi_j(:,j,3).*Ln(3,:,3)+v_dphi_j(:,j,4).*Ln(3,:,4);
        Int_V_B =  V_dual.*sum(w3.*(CofV(:,:,1).*value_dphij1 + CofV(:,:,2).*value_dphij2 + CofV(:,:,3).*value_dphij3),1);
        
        Int_V_c =  V_dual.*sum(w3.*CofV(:,:,4).*v_phi(:,j),1);
        
        for m = 1:size(Vd,2)
            K2((j-1)*dof+Vd(1,m),:) = K2((j-1)*dof+Vd(1,m),:) + Vd(2,m)*Int_V_B;
            K3((j-1)*dof+Vd(1,m),:) = K3((j-1)*dof+Vd(1,m),:) + Vd(2,m)*Int_V_c;
        end
    end
end
K2 = K2./repmat(V,dof^2,1);

K = - K1 + K2 + K3;
end