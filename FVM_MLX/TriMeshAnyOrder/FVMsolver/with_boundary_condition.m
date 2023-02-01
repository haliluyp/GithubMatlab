function [K,F] = with_boundary_condition(K,F,e,inner,Ue)

F(inner,1)=F(inner,1)-K(inner,e)*Ue'; % ���ñ�ֵ�������¸�ֵ�Ҷ˾���B

K(e,:)=[];K(:,e)=[];
F(e,:)=[];
end