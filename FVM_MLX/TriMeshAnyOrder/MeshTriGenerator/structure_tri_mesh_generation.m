function [p,t]=structure_tri_mesh_generation(interval,meshsize,type)
if type ==1
[X,Y] = meshgrid(linspace(interval(1,1),interval(1,2),meshsize(1)+1),linspace(interval(2,1),interval(2,2),meshsize(2)+1));
num_p=(meshsize(1)+1)*(meshsize(2)+1); % .num_p: 剖分得到的总节点数
x=reshape(X,1,num_p);y=reshape(Y,1,num_p);
p=[x;y];% p(2,length): 节点标号及坐标

t=delaunay(x,y)';
% num_t=2*meshsize(1)*meshsize(2);
% for i=1:meshsize(2)
%     T_1st_point(1,1+meshsize(1)*(i-1):meshsize(1)+meshsize(1)*(i-1))=1+(meshsize(1)+1)*(i-1):1+(meshsize(1)+1)*(i-1)+(meshsize(1)-1);
% end
% t(:,1:num_t/2)=[T_1st_point(1,:);T_1st_point(1,:)+meshsize(1)+1;T_1st_point(1,:)+meshsize(1)+2];
% t(:,num_t/2+1:num_t)=[T_1st_point(1,:);T_1st_point(1,:)+meshsize(1)+2;T_1st_point(1,:)+1];
end


%% 第二类结构网格，米字型网格（看作四个矩形形成的大矩形剖分是一致的）
if type ==2
[X,Y] = meshgrid(linspace(interval(1,1),interval(1,2),2*meshsize(1)+1),...
                 linspace(interval(2,1),interval(2,2),2*meshsize(2)+1));
             
num_p = (2*meshsize(1)+1)*(2*meshsize(2)+1); % .num_p: 剖分得到的总节点数
x = reshape(X,1,num_p);y = reshape(Y,1,num_p);
p=[x;y];

% 三角形网格在一个大矩形上的编号
T_point = zeros(9,1);
for i = 1:3
    T_point([1;2;3]+3*(i-1)) = [1;2;3]+(2*meshsize(2)+1)*(i-1);
end

t1(:,1) = [1;4;5];t1(:,2) = [1;5;2];
t2(:,1) = [2;5;3];t2(:,2) = [3;5;6];

t0(:,1:2) = T_point(t1);t0(:,3:4) = T_point(t2);
t0(:,5:6) = T_point(t2+2);t0(:,7:8) = T_point(t1+4);

d_num = size(t0,2);
t_line = reshape(repmat(t0,1,1,meshsize(2))+reshape(0:2:2*meshsize(2)-2,1,1,meshsize(2)),...
    3,d_num*meshsize(2));

t = reshape(repmat(t_line,1,1,meshsize(1))+reshape(0:2*(2*meshsize(2)+1):2*(meshsize(1)-1)*(2*meshsize(2)+1),1,1,meshsize(1)),...
    3,d_num*meshsize(1)*meshsize(2));
end

%% 第三类结构网格，鱼骨型网格（山脉型）
if type ==3
[X,Y] = meshgrid(linspace(interval(1,1),interval(1,2),2*meshsize(1)+1),...
                 linspace(interval(2,1),interval(2,2),2*meshsize(2)+1));
             
num_p = (2*meshsize(1)+1)*(2*meshsize(2)+1); % .num_p: 剖分得到的总节点数
x = reshape(X,1,num_p);y = reshape(Y,1,num_p);
p=[x;y];

% 三角形网格在一个大矩形上的编号
T_point = zeros(9,1);
for i = 1:3
    T_point([1;2;3]+3*(i-1)) = [1;2;3]+(2*meshsize(2)+1)*(i-1);
end

t1(:,1) = [1;4;5];t1(:,2) = [1;5;2];
t2(:,1) = [4;7;5];t2(:,2) = [7;8;5];

t0(:,1:2) = T_point(t1);t0(:,3:4) = T_point(t1+1);
t0(:,5:6) = T_point(t2);t0(:,7:8) = T_point(t2+1);

d_num = size(t0,2);
t_line = reshape(repmat(t0,1,1,meshsize(2))+reshape(0:2:2*meshsize(2)-2,1,1,meshsize(2)),...
    3,d_num*meshsize(2));

t = reshape(repmat(t_line,1,1,meshsize(1))+reshape(0:2*(2*meshsize(2)+1):2*(meshsize(1)-1)*(2*meshsize(2)+1),1,1,meshsize(1)),...
    3,d_num*meshsize(1)*meshsize(2));
end
end