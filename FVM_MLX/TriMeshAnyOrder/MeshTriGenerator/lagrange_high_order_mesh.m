function [p_r,t_r] = lagrange_high_order_mesh(ac,p,t,r)
num_t=size(t,2);
dof = (r+1)*(r+2)/2; % degree of freedom

ac = ac/r;

p_r1 = zeros(2*dof,num_t);
for i=1:dof
    p_r1(2*i-1:2*i,:) = p(:,t(1,:))* ac(1,i) + ...
        p(:,t(2,:))* ac(2,i) + p(:,t(3,:))* ac(3,i);
end

p_r2 = reshape(p_r1,2,dof*num_t);

[p_r3,~,ic] = unique(p_r2','rows');
p_r = p_r3';

t_r=reshape(ic,dof,num_t);

end