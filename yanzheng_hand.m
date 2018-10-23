close all
clc
clear
syms pi L g0 D C tau h

%myrobot入口函数
%加入符号pi就不会出现很大的分数
%拉格朗日快速方法求动力学方程，D-H Standard，以一个五自由度机械臂编程
fold = 'D:\Program Files\MATLAB\R2017a\bin\UserWork\rvctools\myrobot\UVMS';%替换为本m文件所在目录
cd(fold);%定位到当前目录，为dir成功执行做准备

%% 改进D-H使用mykine 1，两者运动学方程不同

%% 标准D-H使用mykine 2
DOF = 4;
alpha = sym('alpha',[1,DOF]);
a = sym('a',[1,DOF]);
d = sym('d',[1,DOF]);
theta = sym('theta',[1,DOF]);

%% 两关节机械臂
a(1) = L;a(2) = L;
d(1) = 0;d(2) = 0;
alpha(1) = 0;alpha(2) = 0;

%% 正运动已验证
%解末端轨迹
% final_p = ones(4,4);
% final_p = 1;
% for i = 1:length(DOF)
%     T = myfkine(alpha(i),a(i),d(i),theta(i));
%     final_p = final_p * T ;
% end
% [m,n] = size(final_p);
% for i = 1:m
%     for j = 1:n
%         final_p(i,j) = simplify(final_p(i,j));
%     end
% end
% display(final_p);
%% 求工作空间
% myworkspace(alpha,a,d);
%% 运动逆解,myikine
% theta = myikine(final_p,L);
% display(theta);
%% 雅可比矩阵的计算
%% 拉格朗日推导动力学方程
%惯量矩阵函数
A = sym('A',[4,4,DOF]);
T = sym('T',[4,4,DOF]);

Ixx = sym('Ixx',[1,DOF]);
Iyy = sym('Iyy',[1,DOF]);
Izz = sym('Izz',[1,DOF]);
Ixy = sym('Ixy',[1,DOF]);
Ixz = sym('Ixz',[1,DOF]);
Iyz = sym('Iyz',[1,DOF]);
m = sym('m',[1,DOF]);
x = sym('x',[1,DOF]);
y = sym('y',[1,DOF]);
z = sym('z',[1,DOF]);
J = sym('J',[4,4,DOF]);

for i=1:DOF
    Ixy(i) = 0;
    Ixz(i) = 0;
    Iyz(i) = 0;
end
for i=1:DOF
    x(i)=-L/2;
    y(i) = 0;
    z(i) = 0;
end
Ixx(1) = 0;Ixy(1) = 0;Ixz(1) = 0;Iyz(1) = 0;
Ixx(2) = 0;Ixy(2) = 0;Ixz(2) = 0;Iyz(2) = 0;
Iyy(1) = (L^2)*m(1)/3;
Izz(1) = (L^2)*m(1)/3;
Iyy(2) = (L^2)*m(2)/3;
Izz(2) = (L^2)*m(2)/3;

%J,A,T
A(:,:,1)= [1,0,0,0;
    0,1,0,0;
    0,0,1,0;
    0,0,0,1];

for i=1:DOF
    J(:,:,i) = [(-Ixx(i)+Iyy(i)+Izz(i))/2, Ixy(i), Ixz(i), m(i)*x(i);
        Ixy(i), (Ixx(i)-Iyy(i)+Izz(i))/2, Iyz(i), m(i)*y(i);
        Ixz(i), Iyz(i), (Ixx(i)+Iyy(i)-Izz(i))/2, m(i)*z(i);
        m(i)*x(i), m(i)*y(i), m(i)*z(i), m(i)];
    A(:,:,i+1) = simplify(myfkine(alpha(i), a(i), d(i), theta(i), 2));
    
end
for i=1:DOF+1
    T(:,:,i) = [1,0,0,0;
        0,1,0,0;
        0,0,1,0;
        0,0,0,1];
    for j=1:i
        T(:,:,i) = simplify(T(:,:,i)*A(:,:,j));
    end
end


% Q,delta_j
Q = sym('Q',[4,4,DOF]);
delta = sym('delta',[4,4,DOF]);
U = sym('U',[4,4,DOF,DOF]);
for i=1:DOF
    Q(:,:,i) = [0,-1,0,0;
        1,0,0,0;
        0,0,0,0;
        0,0,0,0];
end
for j=1:DOF
    delta(:,:,j) =T(:,:,j)*Q(:,:,j)*inv(T(:,:,j));
end

for i=1:DOF
    for j=1:i
        U(:,:,i,j) = simplify( delta(:,:,j)*T(:,:,i+1) );
    end
end
% U(:,:,1,2) = simplify( delta(:,:,2)*T(:,:,2) );

%D
ddt = sym('ddt',[DOF,1]);
DD = sym('DD',[DOF,DOF]);
for i=1:DOF
    for j=1:DOF
        DD(i,j) = 0;
        for p=max([i,j]):DOF
            DD(i,j) = expand( simplify( DD(i,j)+trace(U(:,:,p,j)*J(:,:,p)*transpose(U(:,:,p,i)))) );
        end
    end
end
D = DD*ddt;

%U(i,j,k)
dt = sym('dt',[DOF,1]);
CC = sym('CC',[DOF,DOF,DOF]);
UU = sym('UU',[4,4,DOF,DOF,DOF]);
for i=1:DOF
    for k=1:i
        for j=1:k
            UU(:,:,i,j,k) = simplify( delta(:,:,j)*delta(:,:,k)*T(:,:,i+1) );
            UU(:,:,i,k,j) = UU(:,:,i,j,k);
        end
    end
end

% CC
for i=1:DOF
    for j=1:DOF
        for k=1:DOF
            CC(i,j,k) = 0;
            for p = max([i,j,k]):DOF
                CC(i,j,k) = simplify(CC(i,j,k)+trace(UU(:,:,p,j,k)*J(:,:,p)*transpose(U(:,:,p,i)))) ;
            end
        end
    end
end
C = sym('C',[DOF,1]);
for i = 1:DOF
    C(i) = 0;
    for j=1:DOF
        for k=1:DOF
            C(i) = simplify( C(i) + CC(i,j,k)*dt(j)*dt(k) );
        end
    end
end

%G
G = sym('G',[DOF,1]);
gx = 0;gy = -g0;gz = 0;% g0=9.81重力方向，y轴，D-H为标准方法。
g = [gx,gy,gz,0];
r = sym('r',[4,DOF]);
for i=1:DOF
    r(:,i) = transpose([x(i),y(i),z(i),1]);
end

for i=1:DOF
    G(i) = 0;
    for p=i:DOF
        G(i) = simplify(G(i) - m(p)*( g* simplify(U(:,:,p,i)*r(:,p)) ) );
    end
end

tau = expand( simplify(D + C ) ) + G