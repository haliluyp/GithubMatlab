function [Uh,exact_Ug] = FEMsolver(w3,g3,p,t,e,inner,P1,P2,P3,P4,value_phi,value_dphi_j,Ln,V)
num_i = size(g3,1);
x1 = g3(:,1).*repmat(P1(1,:),num_i,1) + g3(:,2).*repmat(P2(1,:),num_i,1) + ...
     g3(:,3).*repmat(P3(1,:),num_i,1) + g3(:,4).*repmat(P4(1,:),num_i,1);
y1 = g3(:,1).*repmat(P1(2,:),num_i,1) + g3(:,2).*repmat(P2(2,:),num_i,1) + ...
     g3(:,3).*repmat(P3(2,:),num_i,1) + g3(:,4).*repmat(P4(2,:),num_i,1);
z1 = g3(:,1).*repmat(P1(3,:),num_i,1) + g3(:,2).*repmat(P2(3,:),num_i,1) + ...
     g3(:,3).*repmat(P3(3,:),num_i,1) + g3(:,4).*repmat(P4(3,:),num_i,1);

for i = 1:6
    Cof(:,:,i) = exact_functions(x1,y1,z1,'A',i);
end
for i = 1:3
    Cof(:,:,i+6) = exact_functions(x1,y1,z1,'B',i);
end
Cof(:,:,10) = exact_functions(x1,y1,z1,'c',1);

f = exact_functions(x1,y1,z1,'f',1);
for i = 1:4
    exact_Ug(:,:,i) = exact_functions(x1,y1,z1,'u',i);
end
[K] = get_K(value_phi,value_dphi_j,Ln,V,Cof,w3);
num_p=size(p,2);
[K] = assemble_K(K,t,num_p);
[F] = get_F(value_phi,V,f,w3);
[F] = assemble_F(F,t,num_p);

Ue = exact_functions(p(1,e),p(2,e),p(3,e),'u',1);
[K,F] = with_boundary_condition(K,F,e,inner,Ue); 
Uh = zeros(num_p,1);
Uh(inner,1) = K\F;
Uh(e,1) = Ue;

end