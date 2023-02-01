function [Dual] = dual_info(r,ReLoc,d2,d1,dual_para)
% 考虑1/6的三角形

ver_label =  [1, r*(r+1)/2+1, (r+1)*(r+2)/2]; % 三个顶点的序号

% 三角形三个顶点:  对于任意r不变
Dual.N(1,:) = (ReLoc.ac(:,ver_label(d1+1),d2))'; 
Dual.N(2,:) = [1/2,1/2,1/2]*r; Dual.N(2,d2) = 0; 
Dual.N(3,:) = [1/3,1/3,1/3]*r; 

%%%% 
if r ==1
    Dual.L = [2 3];
    
    Dual.S = [1 2 3];
    
    tN = ReLoc.ver(d2,ver_label(d1+1)); % 小三角形中有哪些计算节点（序号）
    
    
    Dual.Ld = cell(1,size(Dual.L,1));
    Dual.Ld{1,1} = [tN(1);
                     (-1)^(d1+1)];
                 
                 
    Dual.Sd = cell(1,size(Dual.S,1));
    Dual.Sd{1,1} = [tN(1);
                     1];
end

if r == 2
    Dual.N(4,:) = 2*dual_para.a*Dual.N(2,:) + (1-2*dual_para.a)*Dual.N(1,:);
    Dual.N(5,:) = 3/2*dual_para.b*Dual.N(3,:) + (1-3/2*dual_para.b)*Dual.N(1,:);
    
    
    Dual.L = [4 5; 5 3];
    
    Dual.S = [1 4 5; 1 2 3];
    
    
    tN = [ReLoc.ver(d2,ver_label(d1+1)),ReLoc.ver(d2,ver_label(2)+1)]; % 小三角形中有哪些计算节点（序号）
    
    Dual.Ld = cell(1,size(Dual.L,1));
    Dual.Ld{1,1} = [tN(1) tN(2);
                   (-1)^(d1+1) -(-1)^(d1+1)];
    Dual.Ld{1,2} = [tN(2);
                   -(-1)^(d1+1)];
                 
    Dual.Sd = cell(1,size(Dual.S,1));
    Dual.Sd{1,1} = [tN(1) tN(2);
                     1    -1];
    Dual.Sd{1,2} = [tN(2);
                     1   ];
end

if r == 3
    Dual.N(4,:) = dual_para.a*Dual.N(2,:) + (1-dual_para.a)*Dual.N(1,:);
    Dual.N(5,:) = dual_para.b*Dual.N(3,:) + (1-dual_para.b)*Dual.N(1,:);
    Dual.N(6,:) = dual_para.c*Dual.N(3,:) + (1-dual_para.c)*Dual.N(2,:);
    
    Dual.L = [4 5; 5 6; 6 2];
    
    Dual.S = [1 4 5; 3 5 6; 1 2 3];
    
    
    tN = [ReLoc.ver(d2,ver_label(d1+1)),ReLoc.ver(d2,ver_label(d1+1)+(-1)^(d1+1)), 5]; % 小三角形中有哪些计算节点（序号）
    
    Dual.Ld = cell(1,size(Dual.L,1));
    Dual.Ld{1,1} = [tN(1) tN(2);
                   (-1)^(d1+1) -(-1)^(d1+1)];
    Dual.Ld{1,2} = [tN(3)    tN(2);
                   (-1)^(d1+1) -(-1)^(d1+1)];
    Dual.Ld{1,3} = [tN(2);
                   -(-1)^(d1+1)];
                 
    Dual.Sd = cell(1,size(Dual.S,1));
    Dual.Sd{1,1} = [tN(1) tN(2);
                     1    -1];
    Dual.Sd{1,2} = [tN(3) tN(2);
                     1    -1];
    Dual.Sd{1,3} = [tN(2);
                     1   ];
end

Dual.N = Dual.N/r;

end