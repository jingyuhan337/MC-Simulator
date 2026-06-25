%%
clc
clear
%%
dim=[0.2 0.2 0.2];%%millimeter size of voxel
% dim1=[1 1 1];%%millimeter
CenterPoint=[0 0 0];
%%

r=0.01;%diameter
r_mean=0.01;
r_std=0.001;
alpha=100;
beta=10^-4;

dis_myocyte1=1*10^-3;%distance between mycytes
dis_myocyte2=1*10^-3;

ele=85:-5:-85;%%orientation
azi=0;

%%

for kkk=1
for kk=1
 kk
ele11=ele(kk);
n_step=1;%1
n_layer=100;
ele1=(ele11-n_layer*n_step:n_step:ele11+n_layer*n_step)*pi/180;
azi1=azi*pi/180;


for k=1:length(ele1)
dir1(k,:)=[cos(ele1(k)).*cos(azi1) cos(ele1(k)).*sin(azi1) sin(ele1(k))];
end
dir2=cross(dir1(1,:),dir1(2,:))/norm(cross(dir1(1,:),dir1(2,:)));


% dir1=[dir111;dir111;dir111;dir111;dir111];
for k=1:size(dir1,1)
dir3(k,:)=cross(dir1(k,:),dir2);
end


MiddlePoint=zeros(1,3);
MiddlePoint1=cell(1,1);
m=1;

for ii=-n_layer:n_layer
MiddlePoint(m,:)=CenterPoint+(dir2*sqrt((2*r+dis_myocyte2)^2-(r+dis_myocyte1/2)^2)+dir3(ii+n_layer+1,:)*(r+dis_myocyte1/2))*ii;
for j=-100:100
    MiddlePoint1{m}(j+101,:)=MiddlePoint(m,:)+dir3(ii+n_layer+1,:)*(2*r+dis_myocyte1)*j;
end
m=m+1;
end

StartPoint=zeros(1,3);
EndPoint=zeros(1,3);
MiddlePoint2=zeros(1,3);
h=zeros(1,1);
n=1;
for i=1:size(MiddlePoint,1)

for k=1:size(MiddlePoint1{i},1)
if dir1(i,1)~=0&&dir1(i,2)==0&&dir1(i,3)~=0

InterPoint=zeros(4,3);
InterPoint(1,:)=PlaneLineInterPoint([1,0,0],dir1(i,:),[CenterPoint(1)-dim(1)/2 0 0],MiddlePoint1{i}(k,:));
InterPoint(2,:)=PlaneLineInterPoint([1,0,0],dir1(i,:),[CenterPoint(1)+dim(1)/2 0 0],MiddlePoint1{i}(k,:));
InterPoint(3,:)=PlaneLineInterPoint([0,0,1],dir1(i,:),[0 0 CenterPoint(3)-dim(3)/2],MiddlePoint1{i}(k,:));
InterPoint(4,:)=PlaneLineInterPoint([0,0,1],dir1(i,:),[0 0 CenterPoint(3)+dim(3)/2],MiddlePoint1{i}(k,:));

InterPoint(1,1)=CenterPoint(1)-dim(1)/2;
InterPoint(2,1)=CenterPoint(1)+dim(1)/2;
InterPoint(3,3)=CenterPoint(3)-dim(3)/2;
InterPoint(4,3)=CenterPoint(3)+dim(3)/2;
InterPoint1=zeros(2,3);
m=1;
for j=1:4
    if sum([InterPoint(j,:)<CenterPoint-dim/2,InterPoint(j,:)>CenterPoint+dim/2])
           continue;
    end
    InterPoint1(m,:)=InterPoint(j,:);
    m=m+1;
end
if m==3

if dot((InterPoint1(2,:)-InterPoint1(1,:))/norm(InterPoint1(2,:)-InterPoint1(1,:)),dir1(i,:))>0.5
    StartPoint(n,:)=InterPoint1(1,:);
    EndPoint(n,:)=InterPoint1(2,:);
    MiddlePoint2(n,:)=MiddlePoint1{i}(k,:);

    h(n)=norm(StartPoint(n,:)-EndPoint(n,:));
%      if h(n)>r
    n=n+1;
%     else
%         continue;
%     end

else
    StartPoint(n,:)=InterPoint1(2,:);
    EndPoint(n,:)=InterPoint1(1,:);
   MiddlePoint2(n,:)=MiddlePoint1{i}(k,:);
    h(n)=norm(StartPoint(n,:)-EndPoint(n,:));
%     if h(n)>r
    n=n+1;
%     else
%         continue;
%     end
end

end
end
if dir1(i,1)~=0&&dir1(i,2)==0&&dir1(i,3)==0

InterPoint=zeros(2,3);
InterPoint(1,:)=PlaneLineInterPoint([1,0,0],dir1(i,:),[CenterPoint(1)-dim(1)/2 0 0],MiddlePoint1{i}(k,:));
InterPoint(2,:)=PlaneLineInterPoint([1,0,0],dir1(i,:),[CenterPoint(1)+dim(1)/2 0 0],MiddlePoint1{i}(k,:));

InterPoint(1,1)=CenterPoint(1)-dim(1)/2;
InterPoint(2,1)=CenterPoint(1)+dim(1)/2;

InterPoint1=zeros(2,3);
m=1;
for j=1:2
    if sum([InterPoint(j,:)<CenterPoint-dim/2,InterPoint(j,:)>CenterPoint+dim/2])
           continue;
    end
    InterPoint1(m,:)=InterPoint(j,:);
    m=m+1;
end
if m==3

if dot((InterPoint1(2,:)-InterPoint1(1,:))/norm(InterPoint1(2,:)-InterPoint1(1,:)),dir1(i,:))>0.5
    StartPoint(n,:)=InterPoint1(1,:);
    EndPoint(n,:)=InterPoint1(2,:);
    MiddlePoint2(n,:)=MiddlePoint1{i}(k,:);

    h(n)=norm(StartPoint(n,:)-EndPoint(n,:));
%      if h(n)>r
    n=n+1;
%     else
%         continue;
%     end

else
    StartPoint(n,:)=InterPoint1(2,:);
    EndPoint(n,:)=InterPoint1(1,:);
   MiddlePoint2(n,:)=MiddlePoint1{i}(k,:);
    h(n)=norm(StartPoint(n,:)-EndPoint(n,:));
%     if h(n)>r
    n=n+1;
%     else
%         continue;
%     end
end

end
end
if dir1(i,1)~=0&&dir1(i,2)~=0&&dir1(i,3)~=0
InterPoint=zeros(6,3);
InterPoint(1,:)=PlaneLineInterPoint([1,0,0],dir1(i,:),[CenterPoint(1)-dim(1)/2 0 0],MiddlePoint1{i}(k,:));
InterPoint(2,:)=PlaneLineInterPoint([1,0,0],dir1(i,:),[CenterPoint(1)+dim(1)/2 0 0],MiddlePoint1{i}(k,:));
  InterPoint(3,:)=PlaneLineInterPoint([0,1,0],dir1(i,:),[0 CenterPoint(2)-dim(2)/2 0],MiddlePoint1{i}(k,:));
  InterPoint(4,:)=PlaneLineInterPoint([0,1,0],dir1(i,:),[0 CenterPoint(2)+dim(2)/2 0],MiddlePoint1{i}(k,:));
InterPoint(5,:)=PlaneLineInterPoint([0,0,1],dir1(i,:),[0 0 CenterPoint(3)-dim(3)/2],MiddlePoint1{i}(k,:));
InterPoint(6,:)=PlaneLineInterPoint([0,0,1],dir1(i,:),[0 0 CenterPoint(3)+dim(3)/2],MiddlePoint1{i}(k,:));

InterPoint(1,1)=CenterPoint(1)-dim(1)/2;
InterPoint(2,1)=CenterPoint(1)+dim(1)/2;
 InterPoint(3,2)=CenterPoint(2)-dim(2)/2;
 InterPoint(4,2)=CenterPoint(2)+dim(2)/2;
InterPoint(5,3)=CenterPoint(3)-dim(3)/2;
InterPoint(6,3)=CenterPoint(3)+dim(3)/2;
InterPoint1=zeros(2,3);
m=1;
for j=1:6
    if sum([InterPoint(j,:)<CenterPoint-dim/2,InterPoint(j,:)>CenterPoint+dim/2])
           continue;
    end
    InterPoint1(m,:)=InterPoint(j,:);
    m=m+1;
end
if m==3

if dot((InterPoint1(2,:)-InterPoint1(1,:))/norm(InterPoint1(2,:)-InterPoint1(1,:)),dir1(i,:))>0.5
    StartPoint(n,:)=InterPoint1(1,:);
    EndPoint(n,:)=InterPoint1(2,:);
    MiddlePoint2(n,:)=MiddlePoint1{i}(k,:);

    h(n)=norm(StartPoint(n,:)-EndPoint(n,:));
%      if h(n)>r
    n=n+1;
%     else
%         continue;
%     end

else
    StartPoint(n,:)=InterPoint1(2,:);
    EndPoint(n,:)=InterPoint1(1,:);
   MiddlePoint2(n,:)=MiddlePoint1{i}(k,:);
    h(n)=norm(StartPoint(n,:)-EndPoint(n,:));
%     if h(n)>r
    n=n+1;
%     else
%         continue;
%     end
end

end

end

end

end

%%%%Consideration of myocyte length
h1=0.06;
size_StartPoint1=0;
StartPoint1=zeros(1,3);
EndPoint1=zeros(1,3);
for j=1:size(StartPoint,1)
    dir11=(EndPoint(j,:)-StartPoint(j,:))/norm(EndPoint(j,:)-StartPoint(j,:));
    num=floor(norm(StartPoint(j,:)-EndPoint(j,:))/h1);
%     dd=(num*h1-norm(StartPoint{k1,k2,k3}(j,:)-EndPoint{k1,k2,k3}(j,:)))/2;
if num==0
    StartPoint1(size_StartPoint1+1,:)=StartPoint(j,:);
EndPoint1(size_StartPoint1+1,:)=EndPoint(j,:);
else
for i=1:num
   
    StartPoint1(size_StartPoint1+i,:)=StartPoint(j,:)+h1*dir11*(i-1);
    EndPoint1(size_StartPoint1+i,:)=StartPoint1(size_StartPoint1+i,:)+h1*dir11;
end
StartPoint1(size_StartPoint1+num+1,:)=StartPoint(j,:)+h1*dir11*num;
EndPoint1(size_StartPoint1+num+1,:)=EndPoint(j,:);
end
size_StartPoint1=size(StartPoint1,1);
end



 %%
%%%%%Consideration of diameter heterogeneity
rt=zeros(1,length(StartPoint));
StartPoint1=zeros(length(StartPoint),3);
EndPoint1=zeros(length(StartPoint),3);

for i=1:length(StartPoint)
    rt(i)=gamrnd(alpha,beta);
end
rt1=sort(rt,'descend');

m=1;
pt=zeros(1,length(StartPoint));
pt(m)=ceil(rand(1)*length(StartPoint));

StartPoint1(m,:)=StartPoint(pt(m),:);
EndPoint1(m,:)=EndPoint(pt(m),:);
 m=2;
for i=1:1000
    while m<=length(StartPoint)
       pt(m) =ceil(rand(1)*length(StartPoint));
          for k=1:m-1
              if pt(m)==pt(k)
                  break;
              end
               dis1=norm(StartPoint(pt(m),:)-StartPoint(pt(k),:));
               dis2=norm(EndPoint(pt(m),:)-EndPoint(pt(k),:));
              if dis1<rt1(m)+rt1(k)||dis2<rt1(m)+rt1(k)
                  break;
              end
              if k==m-1

                 StartPoint1(m,:)=StartPoint(pt(m),:);
                  EndPoint1(m,:)=EndPoint(pt(m),:);
                 m=m+1;
              end

          end
    end

    if m==length(MiddlePoint2)+1
        break;
    end
end


%%
%output
 fid=fopen('local\startpoint_mean_r=0.01_std_r=0.001_step=1_DM=50_0.6_0.6_0.6.txt','w');
fprintf(fid,'%f\n',size(StartPoint1,1));%%myocyte number
 for i=1:size(StartPoint1,1)

  fprintf(fid,'%f\n',rt1(i));

 end
 for i=1:size(StartPoint1,1)
  for j=1:3
  fprintf(fid,'%f\n',StartPoint1(i,j));
  end
 end
 for i=1:size(EndPoint1,1)
  for j=1:3
  fprintf(fid,'%f\n',EndPoint1(i,j));
  end
 end
fclose(fid);
% end
% end
