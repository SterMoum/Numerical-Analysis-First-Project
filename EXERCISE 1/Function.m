x = linspace(-2,2,1000);
y = e.^(sin(x).^3)+x.^6-2*x.^4-x.^3-1;
plot(x,y);
title("Function Plot")
