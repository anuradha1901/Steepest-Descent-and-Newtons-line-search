x = linspace(-20,20);
y = linspace(-20,20);
[X,Y] = meshgrid(x,y);
Z = 5 - 5*X - 2*Y + 2*X.*X + 5*X.*Y + 6*Y.*Y;
contour3(X,Y,Z,1000)
colorbar 
title('Contour function of $f(x,y) = 5-5x-2y+2x^2+5xy+6y^2$','Interpreter','latex')
xlabel('X')
ylabel('Y')
zlabel('Value of the function')