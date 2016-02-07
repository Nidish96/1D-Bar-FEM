set xl 'x'
set yl 'u'
set grid
set zeroax lt -1

P = 1.0;
A(x) = 1.0;
E(x) = 1.0;
L = 1.0;

u(x) = ((P+L**2/2)*x-x**3/6)/(A(x)*E(x));

plot '<./Main' u 2:4 w lp t 'FEM Solution',u(x) w l t 'Actual Solution'
