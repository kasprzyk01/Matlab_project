function dxdt = linearyzacja(t , x , u)

global A;
global B;
dxdt = [ ...
   A(1,1)*x(1) + A(1,2)*x(2) + B(1,1)*u(1) + B(1,2)*u(2) + B(1,3)*u(3) + B(1,4)*u(4);
   A(2,1)*x(1) + A(2,2)*x(2) + B(2,1)*u(1) + B(2,2)*u(2) + B(2,3)*u(3) + B(2,4)*u(4);
];