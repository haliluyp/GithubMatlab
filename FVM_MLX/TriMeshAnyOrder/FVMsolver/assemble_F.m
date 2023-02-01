function [F]=assemble_F(F,t,num_p)
dof = size(t,1);
num_t = size(t,2);
ii=reshape(t',dof*num_t,1);
F=accumarray(ii,F,[num_p,1]);
end