function [p1,t1]=uniform_refinement_3d_tetra(p,t)
nump=size(p,2);numt=size(t,2);

edge1=[t([2 3],:);t([3 4],:);t([4 2],:);t([1 2],:);t([1 3],:);t([1 4],:)]; % 12 * num_t : ÿ������һ������ߵ������˵�Ľڵ���

edge2=reshape(edge1,2,6*numt); %�������б�Ϊ���бߣ�����6������

d_tag=edge2(2,:)-edge2(1,:); %���edge2�� ĩ�˵��ǩ����˵��ǩ��ֵ���ж��Ƿ��С��ǩ�˵㵽���ǩ�˵㣩

w_edge=find(d_tag<0);  % ���edge2�дӴ��ǩ�˵㵽С��ǩ�˵�ı�

% edge3=edge2';
% edge3(w_edge,1)=edge3(w_edge,1)+(d_tag(w_edge))';
% edge3(w_edge,2)=edge3(w_edge,2)-(d_tag(w_edge))'; % edge3: ��edge2 ȫ���������Ϊ ��С��ǩ�˵㵽���ǩ�˵� �ı�

edge2(:,w_edge) = flip(edge2(:,w_edge),1);

[edge4,~,ic]=unique(edge2','rows'); % edge4: �ų��ظ��ߺ�õ��ıߣ�edge3(ia)=edge4 , edge4(ic)=edge3

midp=(p(:,edge4(:,1))+p(:,edge4(:,2)))/2; % edge4�еı��е�

p1=[p,midp]; % һ�¼���һ�κ�õ��Ľڵ�

midp_in_t=reshape(ic,6,numt)+nump; % ��i�������е�Ԫ�ĵ�i���ߵ��е��ǩ����ǩ��Ҫ����ǰ��Ķ����ǩ����

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