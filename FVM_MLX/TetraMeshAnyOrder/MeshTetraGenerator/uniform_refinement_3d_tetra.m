function [p1,t1]=uniform_refinement_3d_tetra(p,t)
nump=size(p,2);numt=size(t,2);

edge1=[t([2 3],:);t([3 4],:);t([4 2],:);t([1 2],:);t([1 3],:);t([1 4],:)]; % 12 * num_t : 每两行是一条有向边的两个端点的节点标号

edge2=reshape(edge1,2,6*numt); %重排六行边为两行边，增加6倍列数

d_tag=edge2(2,:)-edge2(1,:); %输出edge2边 末端点标签减起端点标签的值（判断是否从小标签端点到大标签端点）

w_edge=find(d_tag<0);  % 输出edge2中从大标签端点到小标签端点的边

% edge3=edge2';
% edge3(w_edge,1)=edge3(w_edge,1)+(d_tag(w_edge))';
% edge3(w_edge,2)=edge3(w_edge,2)-(d_tag(w_edge))'; % edge3: 将edge2 全部重新输出为 从小标签端点到大标签端点 的边

edge2(:,w_edge) = flip(edge2(:,w_edge),1);

[edge4,~,ic]=unique(edge2','rows'); % edge4: 排除重复边后得到的边，edge3(ia)=edge4 , edge4(ic)=edge3

midp=(p(:,edge4(:,1))+p(:,edge4(:,2)))/2; % edge4中的边中点

p1=[p,midp]; % 一致加密一次后得到的节点

midp_in_t=reshape(ic,6,numt)+nump; % 第i行是所有单元的第i条边的中点标签，标签需要加上前面的顶点标签总数

% t1=reshape([midp_in_t(4,:);midp_in_t(5,:);midp_in_t(6,:);t(4,:);
%             t(1,:);midp_in_t(1,:);midp_in_t(3,:);midp_in_t(4,:);
%             t(2,:);midp_in_t(2,:);midp_in_t(1,:);midp_in_t(5,:);
%             t(3,:);midp_in_t(3,:);midp_in_t(2,:);midp_in_t(6,:);
%             midp_in_t(2,:);midp_in_t(5,:);midp_in_t(4,:);midp_in_t(1,:);
%             midp_in_t(2,:);midp_in_t(4,:);midp_in_t(3,:);midp_in_t(1,:);
%             midp_in_t(2,:);midp_in_t(4,:);midp_in_t(5,:);midp_in_t(6,:);
%             midp_in_t(2,:);midp_in_t(3,:);midp_in_t(4,:);midp_in_t(6,:)],4,8*numt);
t1=reshape([t(1,:);midp_in_t(4,:);midp_in_t(5,:);midp_in_t(6,:);
            t(2,:);midp_in_t(1,:);midp_in_t(4,:);midp_in_t(3,:);
            t(3,:);midp_in_t(1,:);midp_in_t(2,:);midp_in_t(5,:);
            t(4,:);midp_in_t(2,:);midp_in_t(3,:);midp_in_t(6,:);
            midp_in_t(1,:);midp_in_t(2,:);midp_in_t(5,:);midp_in_t(4,:);
            midp_in_t(1,:);midp_in_t(2,:);midp_in_t(4,:);midp_in_t(3,:);
            midp_in_t(6,:);midp_in_t(2,:);midp_in_t(4,:);midp_in_t(5,:);
            midp_in_t(6,:);midp_in_t(2,:);midp_in_t(3,:);midp_in_t(4,:)],4,8*numt);
end