function [F] = get_F(r,P,Dual)
dof = (r+1)*(r+2)*(r+3)/6;
num_t = size(P,2);
F = zeros(dof,num_t);

v4 = zeros(size(P,1),size(P,2),4);
[w3,g3] = tetra_gausspoints(15);
for v = 1:size(Dual.V,1)
    for k1 = 1:4
            v4(:,:,k1) = Dual.N(Dual.V(v,k1),1)*P(:,:,1) + Dual.N(Dual.V(v,k1),2)*P(:,:,2) + ...
                        Dual.N(Dual.V(v,k1),3)*P(:,:,3) + Dual.N(Dual.V(v,k1),4)*P(:,:,4);
    end
    V_dual(1,:) = tetra_volume(v4);% 
    
    X2(:,:,1) = g3(:,1).*v4(1,:,1) + g3(:,2).*v4(1,:,2) + g3(:,3).*v4(1,:,3) + g3(:,4).*v4(1,:,4); 
    X2(:,:,2) = g3(:,1).*v4(2,:,1) + g3(:,2).*v4(2,:,2) + g3(:,3).*v4(2,:,3) + g3(:,4).*v4(2,:,4); 
    X2(:,:,3) = g3(:,1).*v4(3,:,1) + g3(:,2).*v4(3,:,2) + g3(:,3).*v4(3,:,3) + g3(:,4).*v4(3,:,4); 
    CofV_f = exact_functions(X2,'f',1);
    Int_V = V_dual.*sum(w3.*CofV_f,1);
    
    Vd = Dual.Vd{1,v};
    for m = 1:size(Vd,2)
            F(Vd(1,m),:) = F(Vd(1,m),:) + Vd(2,m)*Int_V;
    end
end

end