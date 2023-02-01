function [K] = get_K(r,P,Ln,S,Dual)
dof = (r+1)*(r+2)/2;
num_t = size(Ln,2);


[w1,g1] = line_gausspoints(5);
[w2,g2] = tri_gausspoints(13);

K1 = zeros(dof^2,num_t);
for l = 1:size(Dual.L,1)
    Dual_Ln = -(Dual.N(Dual.L(l,2),2)-Dual.N(Dual.L(l,1),2))*Ln(:,:,1) + (Dual.N(Dual.L(l,2),1)-Dual.N(Dual.L(l,1),1))*Ln(:,:,2);
    gv1 = g1(:,1).*Dual.N(Dual.L(l,1),:) + g1(:,2).*Dual.N(Dual.L(l,2),:);
    [~,v_dphi_j] = BasisFunctionMain(gv1,r);
    X1(:,:,1) = gv1(:,1)*P(1,:,1) + gv1(:,2)*P(1,:,2) + gv1(:,3)*P(1,:,3);
    X1(:,:,2) = gv1(:,1)*P(2,:,1) + gv1(:,2)*P(2,:,2) + gv1(:,3)*P(2,:,3);
    CofL(:,:,1) = exact_functions(X1,'A',1);
    CofL(:,:,2) = exact_functions(X1,'A',2);
    CofL(:,:,3) = exact_functions(X1,'A',3);
    Ld = Dual.Ld{1,l};
    for j = 1:dof
        value_dphij1 = (v_dphi_j(:,j,1).*Ln(1,:,1)+v_dphi_j(:,j,2).*Ln(1,:,2)+v_dphi_j(:,j,3).*Ln(1,:,3));
        value_dphij2 = (v_dphi_j(:,j,1).*Ln(2,:,1)+v_dphi_j(:,j,2).*Ln(2,:,2)+v_dphi_j(:,j,3).*Ln(2,:,3));
        Int_L = sum(w1.*((CofL(:,:,1).*value_dphij1 + CofL(:,:,2).*value_dphij2).*repmat(Dual_Ln(1,:),length(w1),1) + ...
                         (CofL(:,:,2).*value_dphij1 + CofL(:,:,3).*value_dphij2).*repmat(Dual_Ln(2,:),length(w1),1)),1);
        for m = 1:size(Ld,2)
            K1((j-1)*dof+Ld(1,m),:) = K1((j-1)*dof+Ld(1,m),:) + Ld(2,m)*Int_L;
        end
    end
end
K1 = K1./repmat(S,dof^2,1);

v = zeros(size(P,1),size(P,2),3);
K2 = zeros(dof^2,num_t);
K3 = zeros(dof^2,num_t);
for s = 1:size(Dual.S,1)
    for k1 = 1:3
            v(:,:,k1) = Dual.N(Dual.S(s,k1),1)*P(:,:,1) + Dual.N(Dual.S(s,k1),2)*P(:,:,2) + ...
                        Dual.N(Dual.S(s,k1),3)*P(:,:,3);
    end
    X2(:,:,1) = g2(:,1).*v(1,:,1) + g2(:,2).*v(1,:,2) + g2(:,3).*v(1,:,3);
    X2(:,:,2) = g2(:,1).*v(2,:,1) + g2(:,2).*v(2,:,2) + g2(:,3).*v(2,:,3);
    S_dual(1,:) = abs((v(1,:,2)-v(1,:,1)).*(v(2,:,3)-v(2,:,1))-(v(2,:,2)-v(2,:,1)).*(v(1,:,3)-v(1,:,1))); % 两倍的对偶面积
    CofS(:,:,1) = exact_functions(X2,'B',1);
    CofS(:,:,2) = exact_functions(X2,'B',2);
    CofS(:,:,3) = exact_functions(X2,'c',1);
    gv2 = g2(:,1)*Dual.N(Dual.S(s,1),:) + g2(:,2)*Dual.N(Dual.S(s,2),:) + g2(:,3)*Dual.N(Dual.S(s,3),:);
    [v_phi,v_dphi_j] = BasisFunctionMain(gv2,r);
    Sd = Dual.Sd{1,s};
    for j = 1:dof
        value_dphij1 = (v_dphi_j(:,j,1).*Ln(1,:,1)+v_dphi_j(:,j,2).*Ln(1,:,2)+v_dphi_j(:,j,3).*Ln(1,:,3));
        value_dphij2 = (v_dphi_j(:,j,1).*Ln(2,:,1)+v_dphi_j(:,j,2).*Ln(2,:,2)+v_dphi_j(:,j,3).*Ln(2,:,3));
        Int_S_B =  S_dual.*sum(w2.*(CofS(:,:,1).*value_dphij1 + CofS(:,:,2).*value_dphij2),1);
        
        Int_S_c =  S_dual.*sum(w2.*CofS(:,:,3).*v_phi(:,j),1);
        
        for m = 1:size(Sd,2)
            K2((j-1)*dof+Sd(1,m),:) = K2((j-1)*dof+Sd(1,m),:) + Sd(2,m)*Int_S_B;
            K3((j-1)*dof+Sd(1,m),:) = K3((j-1)*dof+Sd(1,m),:) + Sd(2,m)*Int_S_c;
        end
    end
end
K2 = K2./repmat(S,dof^2,1);

K = K1 + K2 + K3;
end