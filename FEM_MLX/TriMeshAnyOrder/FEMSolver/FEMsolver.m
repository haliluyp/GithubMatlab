function [Uh,exact_Ug] = FEMsolver(w2,g2,p,t,e,inner,P1,P2,P3,value_phi,value_dphi_j,Ln,S)
num_i = size(g2,1);
x1 = g2(:,1).*repmat(P1(1,:),num_i,1) + g2(:,2).*repmat(P2(1,:),num_i,1) + ...
     g2(:,3).*repmat(P3(1,:),num_i,1);
y1 = g2(:,1).*repmat(P1(2,:),num_i,1) + g2(:,2).*repmat(P2(2,:),num_i,1) + ...
     g2(:,3).*repmat(P3(2,:),num_i,1);
Cof(:,:,1) = exact_functions(x1,y1,'A',1);
Cof(:,:,2) = exact_functions(x1,y1,'A',2);
Cof(:,:,3) = exact_functions(x1,y1,'A',3);
Cof(:,:,4) = exact_functions(x1,y1,'B',1);
Cof(:,:,5) = exact_functions(x1,y1,'B',2);
Cof(:,:,6) = exact_functions(x1,y1,'c',1);
f = exact_functions(x1,y1,'f',1);
for i = 1:3
    exact_Ug(:,:,i) = exact_functions(x1,y1,'u',i);
end
[K] = get_K(value_phi,value_dphi_j,Ln,S,Cof,w2);
num_p = size(p,2);
[K] = assemble_K(K,t,num_p);
[F] = get_F(value_phi,S,f,w2);
[F] = assemble_F(F,t,num_p);

Ue = exact_functions(p(1,e),p(2,e),'u',1);
[K,F] = with_boundary_condition(K,F,e,inner,Ue); 
Uh = zeros(num_p,1);
Uh(inner,1) = K\F;
Uh(e,1) = Ue;
end