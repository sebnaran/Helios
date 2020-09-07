MType = 'Trig';
%MType = 'Vor';
%MType = 'Quad'
switch MType
    case 'Vor'
    
    case 'Quad'
    
    case 'Trig'

dx      = [0.101015,0.051886,0.0251418,0.0125255];
ElecErr = [0.011660836811078193,0.0037712742896791342,0.0008068230443378578,0.0002024975599295222];
MagnErr = [0.013888787706116087,0.008368253351531236,0.004707377162019469,0.0024400299532172714];
end

figure(1)
clf
loglog(dx,ElecErr,'p','LineWidth',3,'color','k')
hold on
loglog(dx,MagnErr,'x','LineWidth',3,'color','m')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
loglog(dx,ElecErr,'LineWidth',3,'color','k')
loglog(dx,MagnErr,'LineWidth',3,'color','m')
%leftaxis = 0;
%rightaxis = 10;
%upaxis = 45;
%downaxis = 0;

%axis([leftaxis rightaxis downaxis upaxis])
legend('E','B')
%legend('10 Nodes','50 Nodes','100 Nodes','300 Nodes','600 Nodes','900 Nodes')
%legend('Location','northwest')
