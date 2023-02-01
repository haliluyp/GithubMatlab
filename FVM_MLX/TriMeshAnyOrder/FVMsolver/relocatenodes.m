function [ReLoc] = relocatenodes(ac,r)
dof = (r+1)*(r+2)/2;
ReLoc.ver = zeros(3,dof);
ReLoc.ver(1,:) = 1:dof;

Reloc_v = [1 2 3; 2 3 1; 3 1 2];

ReLoc.ac = zeros(size(ac,1),size(ac,2),3);
ReLoc.ac(:,:,1) = ac;

ac1(Reloc_v(2,:),:) = ac(Reloc_v(1,:),:);
ac2(Reloc_v(3,:),:) = ac(Reloc_v(1,:),:);

ReLoc.ac(:,:,2) = ac1;
ReLoc.ac(:,:,3) = ac2;

for i = 1:dof
    [~,loc1] = ismember((ac1(:,i))',ac','rows');
    [~,loc2] = ismember((ac2(:,i))',ac','rows');
    ReLoc.ver(2,i) = loc1;
    ReLoc.ver(3,i) = loc2;
end

end

