function   [L2_error] = L2_estimate_tetra(value_phi,exact_Ug,Uh,V,w3)

L2_error=0;
for i=1:length(w3)
    L2_error=L2_error+w3(i)*(exact_Ug(i,:)-value_phi(i,:)*Uh).^2;
end
L2_error=sqrt(sum(V.*L2_error));
end