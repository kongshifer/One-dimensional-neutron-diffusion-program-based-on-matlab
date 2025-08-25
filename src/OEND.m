function ONED()
clc;clear;
format long
input=readmatrix('input1.xlsx');

N_max=input(1,3);%最大迭代次数
yip_k = input(1,4);
yip_f = input(1,5);
left  = input(1,6);%left boundary
right = input(1,7);%right boundary
n=input(1,1);%分段数
n_material=input(1,2);
Lx(1:n)=zeros;
Nx(1:n)=zeros;
for i=1:1:n
    Lx(1,i)=input(2,i);% x方向长度
end
for i=1:1:n
    Nx(1,i)=input(3,i);%x方向网格数
end
sumNx = sum(Nx);
Delta(1:sum(Nx)) =zeros;% 网格大小/网格区间长度
tep=1;
for j=1:1:n
    for i=tep:1:(tep+Nx(j)-1)
        Delta(i)= Lx(j)/Nx(j);
    end
    tep=i+1;
end
% 截面参数(假设常数)
Sigma_t_tep(1,n)=zeros;
fission_tep(1,n)=zeros; 
Sigma_s_tep(1,n)=zeros;
for i=1:1:n
    index=4+input(4,i);
    Sigma_t_tep(1,i)=input(index,2);%总截面
    fission_tep(1,i)=input(index,3); %裂变截面乘有效裂变中子数
    Sigma_s_tep(1,i)=input(index,4);%散射截面
end


Sigma_a_tep= Sigma_t_tep-Sigma_s_tep;%吸收截面
Sigma_f_tep= Sigma_a_tep;% 裂变截面

phi_n=ones((sum(Nx)-1),N_max+1);%中子通量初值
M((sumNx-1):(sumNx-1)) = zeros; % 消失项矩阵
F((sumNx-1):(sumNx-1)) = zeros; % 裂变源项矩阵
S= zeros((sumNx-1),(sumNx-1)); % 散射源项矩阵
Q= zeros((sum(Nx)-1),1);
k_n(1:N_max+1)=ones; %有效增殖因子初值

D(1:sumNx)=zeros; %系数和截面分配到每个网格中
Sigma_t(1:sumNx)=zeros;
Sigma_s(1:sumNx)=zeros;
fission(1:sumNx)=zeros;
tep=1;
for j=1:1:n
    for i=tep:1:(tep+Nx(j)-1)
       fission(i)= fission_tep(1,j);
       D(i)= 1/(Sigma_t_tep(1,j)*3);
       Sigma_t(i)= Sigma_t_tep(1,j);
       Sigma_s(i)= Sigma_s_tep(1,j);
    end
    tep=i+1;
end

% xp(1:(sum(Nx)+1))=zeros; %网点位置
% for i=1:1:sum(Nx(1:1))
%     xp(i)=(i-1)*Delta(1);
% end
% for i=(sum(Nx(1:1))+1):1:sum(Nx(1:2))
%     xp(i)=xp(i-1)+Delta(2);
% end


%填迭代矩阵
%左侧边界条件
if left==1
    M(1,1) = Sigma_t(1)*Delta(1)/2 ...
        + Sigma_t(2)*Delta(2) + D(2)/Delta(2) ;
    M(1,2) =  -D(2)/Delta(2);
elseif left==2
    i=1;
    a(i) = -D(i)/Delta(i);
    b(i) = Sigma_t(i)*Delta(i)/2 + ...
        Sigma_t(i+1)*Delta(i+1)/2 + ...
        D(i+1)/Delta(i+1) + D(i)/Delta(i);
    c(i) = -D(i+1)/Delta(i+1);
    M(1,1) = b(i) + c(i)*2*D(i+1)/(Delta(i+1)+2*D(i+1));
    M(1,2) = a(i);
end

%填写M矩阵
a(1:sumNx)=zeros;
b(1:sumNx)=zeros;
c(1:sumNx)=zeros;
for i=2:1:(sumNx-2)
        a(i) = -D(i)/Delta(i);
        b(i) = Sigma_t(i)*Delta(i)/2 + ...
               Sigma_t(i+1)*Delta(i+1)/2 + ...
               D(i+1)/Delta(i+1) + D(i)/Delta(i);
        c(i) = -D(i+1)/Delta(i+1);
        M(i,i-1) = a(i);
        M(i,i) = b(i);
        M(i,i+1) = c(i);
end

%右侧边界条件
if right==1
    M((sumNx-1),(sumNx-1)) = Sigma_t(1)*Delta(1)/2 ...
        + Sigma_t(2)*Delta(2) + D(2)/Delta(2) ;
    M((sumNx-1),(sumNx-2)) =  -D(2)/Delta(2);
elseif right==2
    i=i+1;
    a(i) = -D(i)/Delta(i);
    b(i) = Sigma_t(i)*Delta(i)/2 + ...
        Sigma_t(i+1)*Delta(i+1)/2 + ...
        D(i+1)/Delta(i+1) + D(i)/Delta(i);
    c(i) = -D(i+1)/Delta(i+1);
    M((sumNx-1),(sumNx-1)) = b(i) + c(i)*2*D(i+1)/(Delta(i+1)+2*D(i+1));
    M((sumNx-1),(sumNx-2)) = a(i);
end

%填写F和S矩阵
for i=1:1:(sumNx-1)
    F(i,i) = (fission(i)*Delta(i)+fission(i+1)*Delta(i+1))/2;
end
for i=1:1:(sumNx-1)
    S(i,i) = (Sigma_s(i)*Delta(i)+Sigma_s(i+1)*Delta(i+1))/2;
end

%迭代
tic
for i=1:1:N_max
    f_n = F*phi_n(:,i);
    phi_n(:,i+1) = (M-S) \ (f_n/k_n(i));
    f_nn = F*phi_n(:,i+1);
    k_n(i+1) = norm(f_nn)/(norm(f_n)/k_n(i));
    if ((max(abs((f_nn-f_n)./f_nn)))<yip_f) && (abs((k_n(i+1)-k_n(i))/k_n(i+1))<yip_k)
        disp('符合收敛限的特征值结果：')
        disp(k_n(i+1))
        f_n=f_nn;
        %phi_n=phi_nn;
        disp('总迭代次数:')
        disp(i)
        break
    end
    disp(['第',num2str(i),'次迭代的特征值结果为'])
    disp(k_n(i+1))
    f_n=f_nn;
    %phi_n=phi_nn; 
end
disp('迭代用时')
toc
i_max=i;
if i==N_max
    disp('迭代次数达到最大次数未达到给定收敛标准')
end
disp('done')

DPLOT((sumNx-1):1)=zeros;
DPLOT(1,1)=Delta(1);
for i=2:1:(sumNx-1)
    DPLOT(i,1)=DPLOT(i-1,1)+Delta(i);
end
plot(DPLOT,phi_n(:,i_max+1),'LineWidth',2,'Color',"black");
grid on
xlabel("位置(cm)")
ylabel("中子通量(n·cm^{-2}·s^{-1})")
title("一维直角几何单群特征值问题中子通量分布结果")


gif(DPLOT,phi_n,i_max) %画GIF

kk=k_n(1:(i_max+1));
xlswrite('output-keff.xlsx',kk);
end
