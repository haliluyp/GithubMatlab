function  [L2_error] = L2_estimate_tri(value_phi,exact_U,Uh,S,w2)

L2_error=0;
for i=1:length(w2)
    L2_error=L2_error+w2(i)*(exact_U(i,:)-value_phi(i,:)*Uh).^2;
end
L2_error=sqrt(sum(S.*L2_error));
end