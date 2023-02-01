function  [L2_error] = L2_super_estimate_tetra(value_phi,Uh1,Uh2,V,w3)

L2_error=0;
for i=1:length(w3)
    L2_error=L2_error+w3(i)*(value_phi(i,:)*(Uh1-Uh2)).^2;
end
L2_error=sqrt(sum(V.*L2_error));
end