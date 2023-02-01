function area_coord_p = area_coord(r)
dof = (r+1)*(r+2)*(r+3)/6;
area_coord_p = zeros(4,dof);
%%% r��Lagrange��ֵ�ڵ�������꣺P_1 |  P_2   P_3  P_4|..| P_dof_1 .. P_dof
%%%                      L1      r  |  r-1   r-1  r-1|  |  0            0
%%%                      L2      0  |   1     0    0 |  |  r            0 
%%%                      L3      0  |   0     1    0 |  |  0            0
%%%                      L4      0  |   0     0    1 |  |  0            r

%%%   ����r��Lagrange��ֵ�ڵ�һ����r+1�㣬��i����i���ڵ�����ʼ�ڵ����Ϊ (i-1)*i/2+1

for i = 1:r+1  % r��Lagrange��ֵ�ĵ�i��
    for j = 1:i % ÿ��Ķ�ά�����εĵ�j��
        for k = 1:j % ��ά�����εĵ�j����j���ڵ�
            area_coord_p(1,(i-1)*(i)*(i+1)/6+(j-1)*j/2+k) = r-i+1;
            area_coord_p(2,(i-1)*(i)*(i+1)/6+(j-1)*j/2+k) = i-j;
            area_coord_p(3,(i-1)*(i)*(i+1)/6+(j-1)*j/2+k) = j-k;
            area_coord_p(4,(i-1)*(i)*(i+1)/6+(j-1)*j/2+k) = k-1;
        end
    end
end

end