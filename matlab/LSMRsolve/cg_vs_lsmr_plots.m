
figure;
plot(T,normres-1.5551e4,'k*-','MarkerSize',18,'LineWidth',3);
hold on;
plot([4.5],norm(rescg)-1.5551e4,'rs','MarkerSize',18,'MarkerFaceColor','r');
xlabel('Iteration time (hr)','FontSize',16);
ylabel('norm(y - Ax) - 1.551e4','FontSize',16);


figure;
plot(T,normAtres-3318,'k*-','MarkerSize',18,'LineWidth',3);
hold on;
plot([4.5],norm(Atrescg)-3318,'rs','MarkerSize',18,'MarkerFaceColor','r');
xlabel('Iteration time (hr)','FontSize',16);
ylabel('norm( transpose(A)(y - Ax) )-3318','FontSize',16);


figure;
plot(T,normxdiff,'k*-','MarkerSize',18,'LineWidth',3);
xlabel('Iteration time (hr)','FontSize',16);
ylabel('norm( x(cg) - x(lsmr) )','FontSize',16);
