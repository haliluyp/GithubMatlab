function [p,t] = structure_tetra_mesh_generation(interval,meshsize,mesh_type)
if mesh_type == 1  %  Kuhn's partition
[X,Y] = meshgrid(linspace(interval(1,1),interval(1,2),meshsize(1)+1),...
                 linspace(interval(2,1),interval(2,2),meshsize(2)+1));
num_p1 = (meshsize(1)+1)*(meshsize(2)+1); % .num_p1: 每个xy平面剖分得到的总节点数
x = repmat(reshape(X,1,num_p1),1,meshsize(3)+1);
y = repmat(reshape(Y,1,num_p1),1,meshsize(3)+1);
z = reshape(repmat(linspace(interval(3,1),interval(3,2),meshsize(3)+1),num_p1,1),...
            1,num_p1*(meshsize(3)+1));
p = [x;y;z];% p(3,:): 节点坐标

% 四面体剖分点编号
T_point([1;2]) = [1;2];
T_point([3;4]) = [1;2]+meshsize(1)+1;
T_point([5;6]) = [1;2]+num_p1;
T_point([7;8]) = T_point([5;6])+meshsize(1)+1;

% Type1 
t1(:,1) = T_point([1;3;2;7]);t1(:,2) = T_point([1;2;5;7]);t1(:,3) = T_point([2;6;5;7]);
t1(:,4) = T_point([2;3;4;7]);t1(:,5) = T_point([2;4;8;7]);t1(:,6) = T_point([2;8;6;7]);

d_num=6;
t_line = reshape(repmat(t1,1,1,meshsize(2))+reshape(0:1:meshsize(2)-1,1,1,meshsize(2)),...
    4,d_num*meshsize(2));

t_plane = reshape(repmat(t_line,1,1,meshsize(1))+reshape(0:(meshsize(2)+1):(meshsize(1)-1)*(meshsize(2)+1),1,1,meshsize(1)),...
    4,d_num*meshsize(1)*meshsize(2));

t =  reshape(repmat(t_plane,1,1,meshsize(3))+reshape(0:num_p1:(meshsize(3)-1)*num_p1,1,1,meshsize(3)),...
    4,d_num*meshsize(1)*meshsize(2)*meshsize(3));
end

if mesh_type == 2  % 一对角单独分离处理的立方体六分
[X,Y] = meshgrid(linspace(interval(1,1),interval(1,2),meshsize(1)+1),...
                 linspace(interval(2,1),interval(2,2),meshsize(2)+1));
num_p1 = (meshsize(1)+1)*(meshsize(2)+1); % .num_p1: 每个xy平面剖分得到的总节点数
x = repmat(reshape(X,1,num_p1),1,meshsize(3)+1);
y = repmat(reshape(Y,1,num_p1),1,meshsize(3)+1);
z = reshape(repmat(linspace(interval(3,1),interval(3,2),meshsize(3)+1),num_p1,1),...
            1,num_p1*(meshsize(3)+1));
p = [x;y;z];% p(3,:): 节点坐标

% 四面体剖分点编号
T_point([1;2]) = [1;2];
T_point([3;4]) = [1;2]+meshsize(1)+1;
T_point([5;6]) = [1;2]+num_p1;
T_point([7;8]) = T_point([5;6])+meshsize(1)+1;

% Type2(单独角)
t1(:,1) = T_point([1;3;2;5]);t1(:,2) = T_point([2;5;3;7]);t1(:,3) = T_point([2;6;5;7]);
t1(:,4) = T_point([2;3;4;7]);t1(:,5) = T_point([2;4;6;7]);t1(:,6) = T_point([4;6;7;8]);

d_num=6;
t_line = reshape(repmat(t1,1,1,meshsize(2))+reshape(0:1:meshsize(2)-1,1,1,meshsize(2)),...
    4,d_num*meshsize(2));

t_plane = reshape(repmat(t_line,1,1,meshsize(1))+reshape(0:(meshsize(2)+1):(meshsize(1)-1)*(meshsize(2)+1),1,1,meshsize(1)),...
    4,d_num*meshsize(1)*meshsize(2));

t =  reshape(repmat(t_plane,1,1,meshsize(3))+reshape(0:num_p1:(meshsize(3)-1)*num_p1,1,1,meshsize(3)),...
    4,d_num*meshsize(1)*meshsize(2)*meshsize(3));
end


%% 第二类结构网格，一个正方体分为5个四面体，8个正方体看作一块的一致剖分（当作二次元）
if mesh_type == 3
[X,Y] = meshgrid(linspace(interval(1,1),interval(1,2),2*meshsize(1)+1),...
                 linspace(interval(2,1),interval(2,2),2*meshsize(2)+1));
num_p1 = (2*meshsize(1)+1)*(2*meshsize(2)+1); % .num_p1: 每个xy平面剖分得到的总节点数
x = repmat(reshape(X,1,num_p1),1,2*meshsize(3)+1);
y = repmat(reshape(Y,1,num_p1),1,2*meshsize(3)+1);
z = reshape(repmat(linspace(interval(3,1),interval(3,2),2*meshsize(3)+1),num_p1,1),...
            1,num_p1*(2*meshsize(3)+1));
p = [x;y;z];% p(3,:): 节点坐标

% 四面体剖分点编号
T_point = zeros(27,1);
for i = 1:3
    T_point([1;2;3]+3*(i-1)) = [1;2;3]+(2*meshsize(1)+1)*(i-1);
end
for i = 1:2
    T_point((1:9)'+9*i) = T_point((1:9)')+num_p1*i;
end

t1(:,1) = [1;4;2;10];t1(:,2) = [13;10;14;4];t1(:,3) = [11;14;10;2];
t1(:,4) = [5;2;4;14];t1(:,5) = [10;14;4;2];

t2(:,1) = [13;14;16;4];t2(:,2) = [7;8;4;16];t2(:,3) = [17;16;14;8];
t2(:,4) = [5;4;8;14];t2(:,5) = [4;16;8;14];

t0(:,1:5) = T_point(t1); t0(:,6:10) = T_point(t2-2); t0(:,11:15) = T_point(t2);t0(:,16:20) = T_point(t1+4);
t0(:,21:25) = T_point(t2+6); t0(:,26:30) = T_point(t1+10); t0(:,31:35) = T_point(t1+12);t0(:,36:40) = T_point(t2+10);

d_num = size(t0,2);
t_line = reshape(repmat(t0,1,1,meshsize(2))+reshape(0:2:2*meshsize(2)-2,1,1,meshsize(2)),...
    4,d_num*meshsize(2));

t_plane = reshape(repmat(t_line,1,1,meshsize(1))+reshape(0:2*(2*meshsize(2)+1):2*(meshsize(1)-1)*(2*meshsize(2)+1),1,1,meshsize(1)),...
    4,d_num*meshsize(1)*meshsize(2));

t =  reshape(repmat(t_plane,1,1,meshsize(3))+reshape(0:2*num_p1:2*(meshsize(3)-1)*num_p1,1,1,meshsize(3)),...
    4,d_num*meshsize(1)*meshsize(2)*meshsize(3));
end
end