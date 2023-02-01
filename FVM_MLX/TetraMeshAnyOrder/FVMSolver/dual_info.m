function [Dual] = dual_info(r,ReLoc,d3,d2,d1,dual_para)
% 考虑1/24的四面体
Reloc_v = [1 2 3 4; 2 3 1 4; 3 4 1 2; 4 1 3 2];
Reloc_l = [2 3; 3 1; 1 2];

ver_label = [1, r*(r+1)*(r+2)/6+1, r*(r+1)*(r+2)/6+r*(r+1)/2+1, (r+1)*(r+2)*(r+3)/6]; % 四面体的四个顶点序号

tri_v = ReLoc.ver(d3,ver_label(2):ver_label(4));


% 四面体四个顶点: 对于任意r不变
Dual.N(1,:) = [0,0,0,0]; Dual.N(1,Reloc_v(d3,Reloc_l(d2,d1)+1)) = r;
Dual.N(2,:) = [1/2,1/2,1/2,1/2]*r; Dual.N(2,d3) = 0; Dual.N(2,Reloc_v(d3,d2+1)) = 0; 
Dual.N(3,:) = [1/3,1/3,1/3,1/3]*r; Dual.N(3,d3) = 0;
Dual.N(4,:) = [1/4,1/4,1/4,1/4]*r; 

if r ==1
    Dual.S = [2 3 4]; %有向
    
    Dual.V = [1 2 3 4];
    
    tN = tri_v(Reloc_l(d2,d1));
    
    Dual.Sd = cell(1,size(Dual.S,1));
    Dual.Sd{1,1} = [tN(1);
                    (-1)^(d1+1)];
          
    Dual.Vd = cell(1,size(Dual.V,1));
    Dual.Vd{1,1} = [tN(1);
                     1];
end

if r == 2
    tri_r2 = [tri_v([4 6 5]); tri_v([6 1 3]); tri_v([1 4 2])];
    
    Dual.N(5,:) = 2*dual_para.a*Dual.N(2,:) + (1-2*dual_para.a)*Dual.N(1,:); 
    Dual.N(6,:) = 3/2*dual_para.b*Dual.N(3,:) + (1-3/2*dual_para.b)*Dual.N(1,:);
    Dual.N(7,:) = 4/3*dual_para.c*Dual.N(4,:) + (1-4/3*dual_para.c)*Dual.N(1,:);
    
    Dual.S = [5 6 7;3 4 6; 4 7 6]; %有向
    
    Dual.V = [1 2 3 4; 1 5 6 7];
    
    tN(1) = tri_r2(d2,d1);
    tN(2) = tri_r2(d2,3);
    
    Dual.Sd = cell(1,size(Dual.S,1));
    Dual.Sd{1,1} = [tN(1)  tN(2);
                    (-1)^(d1+1) -(-1)^(d1+1)];
    Dual.Sd{1,2} = [tN(2);
                    -(-1)^(d1+1)];
    Dual.Sd{1,3} = [tN(2);
                    -(-1)^(d1+1)];
          
    Dual.Vd = cell(1,size(Dual.V,1));
    Dual.Vd{1,1} = [tN(2);
                     1];
    Dual.Vd{1,2} = [tN(1) tN(2);
                     1   -1];
    
end

if r == 3
    tri_r3 = [tri_v([7 8 10 9 5]);  tri_v([10 6 1 3 5]); tri_v([1 2 7 4 5])];
    
    Dual.N(5,:) = dual_para.a1*Dual.N(2,:) + (1-dual_para.a1)*Dual.N(1,:); 
    Dual.N(6,:) = dual_para.a2*Dual.N(3,:) + (1-dual_para.a2)*Dual.N(1,:);
    Dual.N(7,:) = dual_para.a3*Dual.N(4,:) + (1-dual_para.a3)*Dual.N(1,:);
    
    Dual.N(8,:) = dual_para.b1*Dual.N(3,:) + (1-dual_para.b1)*Dual.N(1,:); 
    Dual.N(9,:) = dual_para.b2*Dual.N(3,:) + (1-dual_para.b2)*Dual.N(2,:);
    Dual.N(10,:) = dual_para.b3*Dual.N(4,:) + (1-dual_para.b3)*Dual.N(1,:);
    Dual.N(11,:) = dual_para.b4*Dual.N(4,:) + (1-dual_para.b4)*Dual.N(2,:);
    
    Dual.S = [5 6 7; 8 11 10; 8 9 11;
              6 8 7; 8 10 7; 2 11 9]; %有向
    
    Dual.V = [1 2 3 4; 1 5 6 7;
              3 10 11 4; 3 8 11 10; 3 11 8 9];
    
    tN(1) = tri_r3(d2,2*d1-1); tN(2) = tri_r3(d2,2*d1);
    tN(3) = tri_r3(d2,5);
    
    Dual.Sd = cell(1,size(Dual.S,1));
    Dual.Sd{1,1} = [tN(1)  tN(2);
                    (-1)^(d1+1) -(-1)^(d1+1)];
    Dual.Sd{1,2} = [tN(3) tN(2);
                    (-1)^(d1+1) -(-1)^(d1+1)];
    Dual.Sd{1,3} = [tN(3) tN(2);
                    (-1)^(d1+1) -(-1)^(d1+1)];
    Dual.Sd{1,4} = [tN(2);
                    -(-1)^(d1+1)];
    Dual.Sd{1,5} = [tN(2);
                    -(-1)^(d1+1)];
    Dual.Sd{1,6} = [tN(2);
                    -(-1)^(d1+1)];
          
    Dual.Vd = cell(1,size(Dual.V,1));
    Dual.Vd{1,1} = [tN(2);
                     1];
    Dual.Vd{1,2} = [tN(1) tN(2);
                     1   -1];
    Dual.Vd{1,3} = [tN(3) tN(2);
                     1   -1];
    Dual.Vd{1,4} = [tN(3) tN(2);
                     1   -1];
    Dual.Vd{1,5} = [tN(3) tN(2);
                     1   -1];
    
end

Dual.N = Dual.N/r;

end

