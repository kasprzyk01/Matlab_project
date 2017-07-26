clear all;
clc;

A = zeros(2,2);
A(1,1) = -12.6780;
A(1,2) = -0.0949;
A(2,1) = 1.5181e+03;
A(2,2) = 9.4622;

B = zeros(2,4);
B(1,1) = 1;
B(1,2) = 0;
B(1,3) = 0;
B(1,4) = 0;
B(2,1) = 0;
B(2,2) = -6.3038;
B(2,3) = 1;
B(2,4) = 1.8736; 
C = zeros(2,2);
C(1,1)=1;
C(1,2)=0;
C(2,1) = 0;
C(2,2)=1;

ss_stan = ss(A,B,C,0);

mian = [
    1,...
    - A(2,2) - A(1,1),...
    A(1,1)*A(2,2) - A(1,2)*A(2,1)
    ];

licz11 = [
    B(1,1)*C(1,1) + B(2,1)*C(1,2),...
    A(2,1)*B(1,1)*C(1,2) - A(2,2)*B(1,1)*C(1,1) + A(1,2)*B(2,1)*C(1,1) - A(1,1)*B(2,1)*C(1,2)
    ];

licz12 = [
    B(1,2)*C(1,1) + B(2,2)*C(1,2),...
    A(2,1)*B(1,2)*C(1,2) - A(2,2)*B(1,2)*C(1,1) + A(1,2)*B(2,2)*C(1,1) - A(1,1)*B(2,2)*C(1,2)
    ];

licz13 = [
    B(1,3)*C(1,1) + B(2,3)*C(1,2),...
    A(2,1)*B(1,3)*C(1,2) - A(2,2)*B(1,3)*C(1,1) + A(1,2)*B(2,3)*C(1,1) - A(1,1)*B(2,3)*C(1,2)
    ];

licz14 = [
    B(1,4)*C(1,1) + B(2,4)*C(1,2),...
    A(2,1)*B(1,4)*C(1,2) - A(2,2)*B(1,4)*C(1,1) + A(1,2)*B(2,4)*C(1,1) - A(1,1)*B(2,4)*C(1,2)
    ];

licz21 = [
    B(1,1)*C(2,1) + B(2,1)*C(2,2),...
    A(2,1)*B(1,1)*C(2,2) - A(2,2)*B(1,1)*C(2,1) + A(1,2)*B(2,1)*C(2,1) - A(1,1)*B(2,1)*C(2,2)
    ];

licz22 = [
    B(1,2)*C(2,1) + B(2,2)*C(2,2),...
    A(2,1)*B(1,2)*C(2,2) - A(2,2)*B(1,2)*C(2,1) + A(1,2)*B(2,2)*C(2,1) - A(1,1)*B(2,2)*C(2,2)
    ];

licz23 = [
    B(1,3)*C(2,1) + B(2,3)*C(2,2),...
    A(2,1)*B(1,3)*C(2,2) - A(2,2)*B(1,3)*C(2,1) + A(1,2)*B(2,3)*C(2,1) - A(1,1)*B(2,3)*C(2,2)
    ];

licz24 = [
    B(1,4)*C(2,1) + B(2,4)*C(2,2),...
    A(2,1)*B(1,4)*C(2,2) - A(2,2)*B(1,4)*C(2,1) + A(1,2)*B(2,4)*C(2,1) - A(1,1)*B(2,4)*C(2,2)
    ];


tf_system = tf( {licz11, licz12, licz13, licz14; licz21, licz22, licz23, licz24}, ...
    mian)


discrete_tf2 = c2d(tf_system,1,'zoh');
discrete_tf3 = c2d(tf_system,0.25,'zoh');
discrete_tf4 = c2d(tf_system,0.01,'zoh');
discrete_tf5 = c2d(tf_system,0.5,'zoh');
[Y1,T1] = step(tf_system);
[Y3,T3] = step(discrete_tf2);
[Y4,T4] = step(discrete_tf3);
[Y5,T5] = step(discrete_tf4);
[Y6,T6] = step(discrete_tf5);

figure(2);
plot(T1,Y1(:,2,1),'b','Linewidth',2);
hold on;
% stairs(T2,Y2(:,2,1),'r');
 %stairs(T4,Y4(:,2,1),'g','Linewidth',2);
 stairs(T5,Y5(:,2,1),'g','Linewidth',2);
 %stairs(T6,Y6(:,2,1),'g','Linewidth',2);
axis([0 2.4 0 50]);
grid on;
xlabel('Czas [s]');
ylabel('Odpowiedü');
print('Wykres3','-dpng','-r300');
discrete_ss = c2d(ss_stan,0.1,'zoh')


