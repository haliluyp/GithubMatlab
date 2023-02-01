function  [L2_error] = L2_super_estimate_tri(value_phi,Uh1,Uh2,S,w2)

L2_error=0;
for i=1:length(w2)
    L2_error=L2_error+w2(i)*(value_phi(i,:)*(Uh1-Uh2)).^2;
end
L2_error=sqrt(sum(S.*L2_error));
end