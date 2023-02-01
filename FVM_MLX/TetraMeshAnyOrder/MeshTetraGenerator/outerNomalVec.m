function [OutVec] = outerNomalVec(Vec1,Vec2 )
%  根据两个向量Vec1,Vec2（右手Vec1->Vec2），求外法向量OutVec

% Veci(3,num_vec); 

OutVec(1,:) = Vec1(2,:).*Vec2(3,:) - Vec1(3,:).*Vec2(2,:);
OutVec(2,:) = -(Vec1(1,:).*Vec2(3,:) - Vec1(3,:).*Vec2(1,:));
OutVec(3,:) = Vec1(1,:).*Vec2(2,:) - Vec1(2,:).*Vec2(1,:);
end

