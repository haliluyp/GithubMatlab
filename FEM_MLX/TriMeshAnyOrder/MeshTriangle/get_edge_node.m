function [e_n,ic] = get_edge_node(t)
numt=size(t,2);

edge1=[t([1 3 6],:);t([6 5 4],:);t([4 2 1],:)]; % 9 * num_t : 每三行是一条有向边所有节点标号

edge2=reshape(edge1,3,3*numt); %重排六行边为两行边，增加3倍列数

d_tag=edge2(3,:)-edge2(1,:); %输出edge2边 末端点标签减起端点标签的值（判断是否从小标签端点到大标签端点）

w_edge=find(d_tag<0);  % 输出edge2中从大标签端点到小标签端点的边

edge2(:,w_edge) = flip(edge2(:,w_edge),1);

[e_n,~,ic]=unique(edge2','rows'); % edge4: 排除重复边后得到的边，edge3(ia)=edge4 , edge4(ic)=edge3

e_n=e_n';
end