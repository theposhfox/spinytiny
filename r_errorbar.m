
function r_errorbar(x,y,z,colorS)
%plot(x,y, 'color', colorS);
%z=error
hold on;
for i=1:length(x);
   line([x(i),x(i)], [y(i)-z(i), y(i)+z(i)], 'linewidth', 0.5, 'color', colorS);
end
