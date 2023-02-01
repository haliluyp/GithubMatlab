function [Ln,V] = get_mesh_info(P)
% 四面体单元各面2倍面积乘对应单位外法向量：Si*ni

Vec24 = P(:,:,4)-P(:,:,2);Vec23 = P(:,:,3)-P(:,:,2);Vec21 = P(:,:,1)-P(:,:,2); 
Vec13 = P(:,:,3)-P(:,:,1);Vec14 = P(:,:,4)-P(:,:,1); 

Ln(:,:,1)= outerNomalVec(Vec24,Vec23);
Ln(:,:,2)= outerNomalVec(Vec13,Vec14);
Ln(:,:,3)= outerNomalVec(Vec21,Vec24);
Ln(:,:,4)= outerNomalVec(Vec23,Vec21); %  四面体四个面的单位外法向量*两倍对应面面积
             
V = tetra_volume(P); % 四面体体积的六倍 6V

end