function [F] = get_F(r,P,Dual)
dof = (r+1)*(r+2)/2;
num_t = size(P,2);
F = zeros(dof,num_t);

v = zeros(size(P,1),size(P,2),3);
[w2,g2] = tri_gausspoints(13);
for s = 1:size(Dual.S,1)
    for k1 = 1:3
            v(:,:,k1) = Dual.N(Dual.S(s,k1),1)*P(:,:,1) + Dual.N(Dual.S(s,k1),2)*P(:,:,2) + ...
                        Dual.N(Dual.S(s,k1),3)*P(:,:,3);
    end
    X2(:,:,1) = g2(:,1).*v(1,:,1) + g2(:,2).*v(1,:,2) + g2(:,3).*v(1,:,3);
    X2(:,:,2) = g2(:,1).*v(2,:,1) + g2(:,2).*v(2,:,2) + g2(:,3).*v(2,:,3);
    S_dual(1,:) = abs((v(1,:,2)-v(1,:,1)).*(v(2,:,3)-v(2,:,1))-(v(2,:,2)-v(2,:,1)).*(v(1,:,3)-v(1,:,1))); % 两倍的对偶面积
    CofS_f = exact_functions(X2,'f',1);
    Sd = Dual.Sd{1,s};
    
    Int_S = S_dual.*sum(w2.*CofS_f,1);
    for m = 1:size(Sd,2)
            F(Sd(1,m),:) = F(Sd(1,m),:) + Sd(2,m)*Int_S;
    end
end

end