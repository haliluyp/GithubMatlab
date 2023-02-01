function [e_n,ic] = get_edge_node(t)
numt=size(t,2);

edge1=[t([1 3 6],:);t([6 5 4],:);t([4 2 1],:)]; % 9 * num_t : ÿ������һ����������нڵ���

edge2=reshape(edge1,3,3*numt); %�������б�Ϊ���бߣ�����3������

d_tag=edge2(3,:)-edge2(1,:); %���edge2�� ĩ�˵��ǩ����˵��ǩ��ֵ���ж��Ƿ��С��ǩ�˵㵽���ǩ�˵㣩

w_edge=find(d_tag<0);  % ���edge2�дӴ��ǩ�˵㵽С��ǩ�˵�ı�

edge2(:,w_edge) = flip(edge2(:,w_edge),1);

[e_n,~,ic]=unique(edge2','rows'); % edge4: �ų��ظ��ߺ�õ��ıߣ�edge3(ia)=edge4 , edge4(ic)=edge3

e_n=e_n';
end