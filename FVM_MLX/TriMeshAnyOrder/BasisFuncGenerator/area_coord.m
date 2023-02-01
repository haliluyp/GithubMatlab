function area_coord_p = area_coord(r)
dof = (r+1)*(r+2)/2;
area_coord_p = zeros(3,dof);
%%% r��Lagrange��ֵ�ڵ�������꣺P_1 |  P_2    P_3 |..| P_dof-r .. P_dof
%%%                  lamda1      r  | r-1    r-1  |  |  0            0
%%%                  lamda2      0  |   1      0  |  |  r            0 
%%%                  lamda3      0  |   0      1  |  |  0            r

%%%   ����r��Lagrange��ֵ�ڵ�һ����r+1�㣬��i����i���ڵ�����ʼ�ڵ����Ϊ (i-1)*i/2+1

for i = 1:r+1  % r��Lagrange��ֵ�ĵ�i��
    for j = 1:i  % ÿ����i���ڵ� 
        area_coord_p(1,(i)*(i-1)/2+j) = r-i+1;
        area_coord_p(2,(i)*(i-1)/2+j) = i-j;
        area_coord_p(3,(i)*(i-1)/2+j) = j-1;
    end
end

end