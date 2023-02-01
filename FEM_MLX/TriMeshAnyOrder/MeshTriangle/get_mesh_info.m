function [Ln,S] = get_mesh_info(P1,P2,P3)
% 三角形单元各边长度乘对应单位外法向量：li*ni
Ln(:,:,1)=-[P3(2,:)-P2(2,:);P2(1,:)-P3(1,:)];
Ln(:,:,2)=-[P1(2,:)-P3(2,:);P3(1,:)-P1(1,:)];
Ln(:,:,3)=-[P2(2,:)-P1(2,:);P1(1,:)-P2(1,:)]; 


S=abs((P2(1,:)-P1(1,:)).*(P3(2,:)-P1(2,:))-...
      (P2(2,:)-P1(2,:)).*(P3(1,:)-P1(1,:))); % 三角形单元面积的两倍
end