function [K]=assemble_K(K,t,num_p)
dof = size(t,1);
num_t = size(t,2);
ii = reshape((repmat(t,dof,1))',1,dof^2*num_t);
jj = reshape((repmat(t,1,dof))',1,dof^2*num_t);
K = sparse(ii,jj,K,num_p,num_p);
end
