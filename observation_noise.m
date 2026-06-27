clc
clear
%%

%creat signal in 6/15/21 direction
% u6=[1/sqrt(2) 0 1/sqrt(2);
%     -1/sqrt(2) 0 1/sqrt(2);
%     0 1/sqrt(2) 1/sqrt(2);
%     0 1/sqrt(2) -1/sqrt(2);
%     1/sqrt(2) 1/sqrt(2) 0;
%     -1/sqrt(2) 1/sqrt(2) 0];
s=2*cos(pi/5);
A=1/sqrt(1+s^2);
B=s/sqrt(1+s^2);
C=(s-1)/2;
D=s/2;
u6=[0 A -B;0 -A -B;A -B 0;A B 0;-B 0 A; -B 0 -A];
u15=[0 0 1;0 1 0;C -D -0.5;C D -0.5;-C -D -0.5;C -D 0.5;-0.5 C D;0.5 -C D;0.5 C -D;0.5 C D;-D 0.5 C;D -0.5 C;D 0.5 -C;-D -0.5 -C;-1 0 0];
u21=[0 0 1;0 A -B;0 -A -B;0 1 0;C -D -0.5;C D -0.5;-C -D -0.5;C -D 0.5;-0.5 C D;0.5 -C D;0.5 C -D;0.5 C D;A -B 0;A B 0;-D 0.5 C;D -0.5 C;D 0.5 -C;-D -0.5 -C;-B 0 A; -B 0 -A;-1 0 0];
% [u162,fcs] = icosahedron(2);
%%
B1=[100 300 500 1000 2000 3000];%s/mm^2
thegma1=[2 5 10 20 30];%%ms

Delta=[20 40 60 80 100];%ms
for t1=1:length(B1)
    for t2=1:length(thegma1)
        for t3=1:length(Delta)
qq1(t1,t2,t3)=sqrt(B1(t1)*10^3./(Delta(t3)-thegma1(t2)/3));
q1(t1,t2,t3)=qq1(t1,t2,t3)/(2*pi);
        end
    end
end

%%
 mean_MD=cell(length(B1),length(thegma1),length(Delta));
mean_FA=cell(length(B1),length(thegma1),length(Delta));
mean_azi=cell(length(B1),length(thegma1),length(Delta));
mean_ele=cell(length(B1),length(thegma1),length(Delta));
mean_EigenValues1=cell(length(B1),length(thegma1),length(Delta));

mean_EigenValues2=cell(length(B1),length(thegma1),length(Delta));
mean_EigenValues3=cell(length(B1),length(thegma1),length(Delta));
% Dif=cell(length(B1),length(thegma1),length(Delta));
% sim_vector1=cell(length(B1),length(thegma1),length(Delta));
S0_reference=1;
% The constant phase is optional. It has no impact on the Rician magnitude.
s0=exp(1i*pi/4);
% SNR=20;
% sigma = 1 / SNR;
SNR_b0_dB=20;
sigma = S0_reference / (10^(SNR_b0_dB/20));
u=u21;
for k=1
for kkk=7
    a=k
    aa=kkk
fileID=fopen(['D:\postdoc\paper\A Robust Multi-Scale Computational Framework\data\startpoint_mean_r=0.01_std_r=0.001_step=1_DM=70_azi=0_ele=-81.25_12.5_81.25\',num2str(k),'_',num2str(kkk),'_startpoint_mean_r=0.01_std_r=0.001_step=1_DM=70_0.4_0.4_0.4_p=0.02_molecules=25000_DELTA=0.001_tau=0.001_D=1_2.5_cycle=130.txt']);
% fileID=fopen('D:\doctorat\programme\C\inmyocyte_outmyocyte_permeability_heterogeneity2\inmyocyte_outmyocyte_permeability_heterogeneity2\r_mean=10_std=1_varied angle_0.7_0.7_0.7_new\8_r_mean=0.01_r_std=0.001_DM=50_0.7_0.7_0.7_per=0.002_molecules=50000_DELTA=0.001_tau=0.1_D=1_2.5_cycle=120_newnew.txt');
 aaa=textscan(fileID,'%f');
  fclose(fileID);  %%
n_molecule=25000;
n_step=131;%%
step_location=cell(1,n_molecule);
for i=1:n_molecule
    for j=1:n_step

       step_location{i}(j,:)=aaa{1}(3*n_step*(i-1)+3*j-2:3*n_step*(i-1)+3*j);
    end
end
for t1=1:length(B1)
for t2=1:length(thegma1)
for t3=1:length(Delta)
   
signal=zeros(1,length(u));
displacement=zeros(1,n_molecule);
for j=1:length(u)
        for i=1:n_molecule

        displacement(i)=dot((sum(step_location{i}(Delta(t3)+1:Delta(t3)+thegma1(t2),:))-sum(step_location{i}(1:thegma1(t2),:)))/(thegma1(t2)),u(j,:));
        end
        signal(j)=(1/n_molecule)*sum(cos(qq1(t1,t2,t3)*displacement));
end

signal_noisy=zeros(100,length(u));
for rep=1:100



    N1 = sigma .* randn(1,length(u));
    N2 = sigma .* randn(1,length(u));
    Sstar = s0 .* signal(:).';
    % Rician magnitude signal
    signal_noisy(rep,:) = abs(Sstar + N1 + 1i .* N2);




DTIdata=struct();
n_voxel=1;

for j=1:length(u)

    DTIdata(j).VoxelData = abs(signal_noisy(rep,j));
    DTIdata(j).Gradient = u(j,:);
    DTIdata(j).Bvalue=B1(t1);%%%(s/mm^2)

end

 [MDv,FAv,VectorF,Y,DifT]=DTI(DTIdata);
 % Dif{t1,t2,t3}(k,kkk,:,:)=DifT;
 MD{t1,t2,t3}(rep)=MDv;
 FA{t1,t2,t3}(rep)=FAv;
[EigenVectors,D]=eig(DifT);
EigenValues=diag(D);
EigenValues1{t1,t2,t3}(rep)=EigenValues(1);
EigenValues2{t1,t2,t3}(rep)=EigenValues(2);
EigenValues3{t1,t2,t3}(rep)=EigenValues(3);
VectorF1=reshape(VectorF,3,3);
Fiber_dir=VectorF1(:,end);%%fiber direction
sim_vector=squeeze(Fiber_dir)';
% sim_vector1(rep)=sim_vector;
%ele~[-1/2*pi,1/2*pi]
sim_ele{t1,t2,t3}(rep)=asin(sim_vector(3))*180/pi;

%azi~[0,181]
if sim_vector(1)>=0
  sim_azi{t1,t2,t3}(rep)=atan(sim_vector(2)/sim_vector(1))*180/pi;
else
  sim_azi{t1,t2,t3}(rep)=(atan(sim_vector(2)/sim_vector(1))+pi)*180/pi;
end
if  sim_azi{t1,t2,t3}(rep)>90
    sim_azi{t1,t2,t3}(rep)=sim_azi{t1,t2,t3}(rep)-180;
end

end
mean_MD{t1,t2,t3}(k,kkk)=mean(MD{t1,t2,t3});
mean_FA{t1,t2,t3}(k,kkk)=mean(FA{t1,t2,t3});
mean_EigenValues1{t1,t2,t3}(k,kkk)=mean(EigenValues1{t1,t2,t3});
mean_EigenValues2{t1,t2,t3}(k,kkk)=mean(EigenValues2{t1,t2,t3});
mean_EigenValues3{t1,t2,t3}(k,kkk)=mean(EigenValues3{t1,t2,t3});
mean_azi{t1,t2,t3}(k,kkk)=mean(sim_azi{t1,t2,t3});
mean_ele{t1,t2,t3}(k,kkk)=mean(sim_ele{t1,t2,t3});

end
end
end
end
end
%%
x=[20 40 60 80 100];
lineColor = {
    [0.74,0.22,0.30],  
    [0.82,0.40,0.20],  
    [0.18,0.64,0.60],   
    [0.12,0.58,0.70],   
    [0.10,0.36,0.66],   
    [0.28,0.30,0.62],  
    [0.48,0.26,0.52]}; 
fillColor = {
    [0.95,0.74,0.76],
    [0.97,0.82,0.74],
    [0.70,0.88,0.86],
    [0.68,0.86,0.92],
    [0.66,0.78,0.92],
    [0.76,0.78,0.90],
    [0.84,0.76,0.88]};

figure('Color','w','Position',[100,100,1200,700]);

hold on;grid on;grid minor;
set(gca,'GridAlpha',0.2);

% 循环绘制6条曲线+阴影带
for k = 1:6
    % 生成模拟数据
    for kk=2
        for kkk=1:5
    mu(k,kkk) = mean(reshape(FA{k,kk,kkk},1,100));
    stdVal(k,kkk) = std(reshape(FA{k,kk,kkk},1,100));
        end
    end
    % 填充阴影
    fill([x,flip(x)],[mu(k,:)-stdVal(k,:),flip(mu(k,:)+stdVal(k,:))],...
        fillColor{k},'FaceAlpha',0.8,'EdgeColor','none');
    % 均值曲线
    plot(x,mu(k,:),'Color',lineColor{k},'LineWidth',2.2);
end

 % ylim([0.1,0.6]);
 xticks([20, 40, 60, 80, 100]);
% yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]);
% xlabel('X');ylabel('Value');
% title('7 Groups Mean Curves with Shaded Std Band');
set(gca,'FontSize',25,'FontWeight','bold','Box','off');
hold off;

%%
error_mu_MD=abs(mu_MD-mu_MD_noise)./mu_MD;
error_mu_FA=abs(mu_FA-mu_FA_noise)./mu_FA;
 %%
 for k = 1:6
    % 生成模拟数据
    for kk=1:5
        for kkk=1:5
    
    std_FA(k,kk,kkk) = std(reshape(FA{k,kk,kkk},1,100));
        end
    end
 end
    %%
% mean_bias_MD=mean_bias_MD_model4_vb;
        % data=squeeze(std_FA(3,:,:))';
     % data=squeeze(std_FA(1:6,:,2))';
  data=squeeze(std_FA(1:6,2,:));
data1=round(data,1);
% figure;heatmap(data1,'CellLabelColor','none','ColorLimits',[0.006 0.012]);
% figure;heatmap(data1,'FontSize',18,'Colormap',parula,'ColorLimits',[0 1.5]);
figure;heatmap(data1,'FontSize',18,'ColorLimits',[0 0.3],'Colormap',parula);