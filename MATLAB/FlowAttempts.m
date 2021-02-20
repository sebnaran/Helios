vordx    = [0.333333,0.12803687993289598,0.06772854614785964]/2;
quaddx   = [0.16666666666666666,0.08333333333333333,0.043478260869565216]/2;
trigdx   = [0.2,0.10101525445522107,0.05018856132284956]/2;
h        = quaddx;
Eq = 'StationaryStokes';
%Eq = 'StationaryMHDStokes';
%Eq = 'InEvoStokes';
%Eq = 'InEvoMHDStokes';
switch Eq
      case 'StationaryStokes'
voruerr  = [0.12794930927288048,0.0521724611922912,0.05266079931244949];
triguerr = [0.03759459564525446,0.009306522093831369,0.004895648724281222];
quaduerr = [0.04213950185413305,0.041612290938811004,0.05170521058942315];
vorperr  = [0.446902355391597,0.013092772500971833,0.00489457695521403];
trigperr = [0.03092108703746163,0.010514173511617892,0.004369464204592047];
quadperr = [0.004369464204592047,0.016578382200688704,0.008761431498975186];
set(gca,'FontSize',15)

%Plotting the presure
figure(1)
clf
loglog(vordx,vorperr,'s','LineWidth',3,'color','r');
hold on
loglog(quaddx,quadperr,'o','LineWidth',3,'color','k');
loglog(trigdx,trigperr,'d','Linewidth',3,'color','b');
loglog(vordx,vorperr,'LineWidth',2,'Color','r');
loglog(quaddx,quadperr,'--','LineWidth',2,'Color','k');
loglog(trigdx,trigperr,':','Linewidth',2,'Color','b');

set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold on
%Pick a basis point for the triangle
%xseed = 0.9*h(2)+h(1)*0.1;
%yseed = 0.00005;
%desiredSlope Of triangle
%slope = 2;
%Another x point
%xnext = h(1);


%b = log10(yseed/(xseed^slope));


%ynext = 10^(b)*xnext^slope;


%x = [xseed, xseed, xnext, xseed];
%y = [yseed, ynext, ynext, yseed];

%loglog(x,y,'LineWidth',2,'Color','k')

%text(xseed*0.7+0.3*xnext,ynext+0.003,'1','Fontsize',15)
%text(xseed-0.001,0.5*yseed+0.5*ynext,'2','Fontsize',15)





%title('Triangular Mesh')





%Now the ticks
%The x-axis
%leftaxis = 0.005;
%rightaxis = 0.13;
%upaxis = 0.012;
%downaxis = 10e-6;

%axis([leftaxis rightaxis downaxis upaxis])

legend('Vor','Quad', 'Trig')
legend('Location','northwest')



%Plotting The velocity Field
figure(2)
clf
loglog(vordx,voruerr,'s','LineWidth',3,'color','r');
hold on
loglog(quaddx,quaduerr,'o','LineWidth',3,'color','k');
loglog(trigdx,triguerr,'d','Linewidth',3,'color','b');
loglog(vordx,voruerr,'LineWidth',2,'Color','r');
loglog(quaddx,quaduerr,'--','LineWidth',2,'Color','k');
loglog(trigdx,triguerr,':','Linewidth',2,'Color','b');

set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold on
%Pick a basis point for the triangle
%xseed = 0.9*h(2)+h(1)*0.1;
%yseed = 0.00005;
%desiredSlope Of triangle
%slope = 2;
%Another x point
%xnext = h(1);


%b = log10(yseed/(xseed^slope));


%ynext = 10^(b)*xnext^slope;


%x = [xseed, xseed, xnext, xseed];
%y = [yseed, ynext, ynext, yseed];

%loglog(x,y,'LineWidth',2,'Color','k')

%text(xseed*0.7+0.3*xnext,ynext+0.003,'1','Fontsize',15)
%text(xseed-0.001,0.5*yseed+0.5*ynext,'2','Fontsize',15)





%title('Triangular Mesh')





%Now the ticks
%The x-axis
%leftaxis = 0.005;
%rightaxis = 0.13;
%upaxis = 0.012;
%downaxis = 10e-6;

%axis([leftaxis rightaxis downaxis upaxis])

legend('Vor','Quad', 'Trig')
legend('Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'StationaryMHDStokes'
 
voruerr  = [1,1,1];
triguerr = [1,1,1];
quaduerr = [0.4405982709,0.46497373,0.49033486657];
vorperr  = [1,1,1];
trigperr = [1,1,1];
quadperr = [26.3214467,60.357694,111.931484117];
set(gca,'FontSize',15)

%Plotting the presure
figure(1)
clf
loglog(vordx,vorperr,'s','LineWidth',3,'color','r');
hold on
loglog(quaddx,quadperr,'o','LineWidth',3,'color','k');
loglog(trigdx,trigperr,'d','Linewidth',3,'color','b');
loglog(vordx,vorperr,'LineWidth',2,'Color','r');
loglog(quaddx,quadperr,'--','LineWidth',2,'Color','k');
loglog(trigdx,trigperr,':','Linewidth',2,'Color','b');

set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold on
%Pick a basis point for the triangle
%xseed = 0.9*h(2)+h(1)*0.1;
%yseed = 0.00005;
%desiredSlope Of triangle
%slope = 2;
%Another x point
%xnext = h(1);


%b = log10(yseed/(xseed^slope));


%ynext = 10^(b)*xnext^slope;


%x = [xseed, xseed, xnext, xseed];
%y = [yseed, ynext, ynext, yseed];

%loglog(x,y,'LineWidth',2,'Color','k')

%text(xseed*0.7+0.3*xnext,ynext+0.003,'1','Fontsize',15)
%text(xseed-0.001,0.5*yseed+0.5*ynext,'2','Fontsize',15)





%title('Triangular Mesh')





%Now the ticks
%The x-axis
%leftaxis = 0.005;
%rightaxis = 0.13;
%upaxis = 0.012;
%downaxis = 10e-6;

%axis([leftaxis rightaxis downaxis upaxis])

legend('Vor','Quad', 'Trig')
legend('Location','northwest')



%Plotting The velocity Field
figure(2)
clf
loglog(vordx,voruerr,'s','LineWidth',3,'color','r');
hold on
loglog(quaddx,quaduerr,'o','LineWidth',3,'color','k');
loglog(trigdx,triguerr,'d','Linewidth',3,'color','b');
loglog(vordx,voruerr,'LineWidth',2,'Color','r');
loglog(quaddx,quaduerr,'--','LineWidth',2,'Color','k');
loglog(trigdx,triguerr,':','Linewidth',2,'Color','b');

set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold on
%Pick a basis point for the triangle
%xseed = 0.9*h(2)+h(1)*0.1;
%yseed = 0.00005;
%desiredSlope Of triangle
%slope = 2;
%Another x point
%xnext = h(1);


%b = log10(yseed/(xseed^slope));


%ynext = 10^(b)*xnext^slope;


%x = [xseed, xseed, xnext, xseed];
%y = [yseed, ynext, ynext, yseed];

%loglog(x,y,'LineWidth',2,'Color','k')

%text(xseed*0.7+0.3*xnext,ynext+0.003,'1','Fontsize',15)
%text(xseed-0.001,0.5*yseed+0.5*ynext,'2','Fontsize',15)





%title('Triangular Mesh')





%Now the ticks
%The x-axis
%leftaxis = 0.005;
%rightaxis = 0.13;
%upaxis = 0.012;
%downaxis = 10e-6;

%axis([leftaxis rightaxis downaxis upaxis])

legend('Vor','Quad', 'Trig')
legend('Location','northwest')

    case 'InEvoStokes'
vordx    = [0.12803687993289598,0.06772854614785964,0.0345033]/2;
quaddx   = [0.16666666666666666,0.08333333333333333,0.043478260869565216]/2;
trigdx   = [0.2,0.10101525445522107,0.05018856132284956]/2;    
 
voruerr  = [0.9727320998600277,0.037172623547815996,0.030170996957118477];
triguerr = [0.4888166877779046,0.008977497217483045,0.020830686482796547];
quaduerr = [0.022786335378031063,0.009227023868856113,0.10670777946340365];
vorperr  = [0.9774339262503107,0.9827250264558244,0.9842195287340206];
trigperr = [0.9593232815161192,0.979951907423957,0.9834830094666028];
quadperr = [0.9724366350438928,0.9816551236589589,0.9839169988252413];
set(gca,'FontSize',15)

%Plotting the presure
figure(1)
clf
loglog(vordx,vorperr,'s','LineWidth',3,'color','r');
hold on
loglog(quaddx,quadperr,'o','LineWidth',3,'color','k');
loglog(trigdx,trigperr,'d','Linewidth',3,'color','b');
loglog(vordx,vorperr,'LineWidth',2,'Color','r');
loglog(quaddx,quadperr,'--','LineWidth',2,'Color','k');
loglog(trigdx,trigperr,':','Linewidth',2,'Color','b');

set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold on
%Pick a basis point for the triangle
%xseed = 0.9*h(2)+h(1)*0.1;
%yseed = 0.00005;
%desiredSlope Of triangle
%slope = 2;
%Another x point
%xnext = h(1);


%b = log10(yseed/(xseed^slope));


%ynext = 10^(b)*xnext^slope;


%x = [xseed, xseed, xnext, xseed];
%y = [yseed, ynext, ynext, yseed];

%loglog(x,y,'LineWidth',2,'Color','k')

%text(xseed*0.7+0.3*xnext,ynext+0.003,'1','Fontsize',15)
%text(xseed-0.001,0.5*yseed+0.5*ynext,'2','Fontsize',15)





%title('Triangular Mesh')





%Now the ticks
%The x-axis
%leftaxis = 0.005;
%rightaxis = 0.13;
%upaxis = 0.012;
%downaxis = 10e-6;

%axis([leftaxis rightaxis downaxis upaxis])

legend('Vor','Quad', 'Trig')
legend('Location','northwest')



%Plotting The velocity Field
figure(2)
clf
loglog(vordx,voruerr,'s','LineWidth',3,'color','r');
hold on
loglog(quaddx,quaduerr,'o','LineWidth',3,'color','k');
loglog(trigdx,triguerr,'d','Linewidth',3,'color','b');
loglog(vordx,voruerr,'LineWidth',2,'Color','r');
loglog(quaddx,quaduerr,'--','LineWidth',2,'Color','k');
loglog(trigdx,triguerr,':','Linewidth',2,'Color','b');

set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold on
%Pick a basis point for the triangle
%xseed = 0.9*h(2)+h(1)*0.1;
%yseed = 0.00005;
%desiredSlope Of triangle
%slope = 2;
%Another x point
%xnext = h(1);


%b = log10(yseed/(xseed^slope));


%ynext = 10^(b)*xnext^slope;


%x = [xseed, xseed, xnext, xseed];
%y = [yseed, ynext, ynext, yseed];

%loglog(x,y,'LineWidth',2,'Color','k')

%text(xseed*0.7+0.3*xnext,ynext+0.003,'1','Fontsize',15)
%text(xseed-0.001,0.5*yseed+0.5*ynext,'2','Fontsize',15)





%title('Triangular Mesh')





%Now the ticks
%The x-axis
%leftaxis = 0.005;
%rightaxis = 0.13;
%upaxis = 0.012;
%downaxis = 10e-6;

%axis([leftaxis rightaxis downaxis upaxis])

legend('Vor','Quad', 'Trig')
legend('Location','northwest')

    case 'InEvoMHDStokes'
voruerr  = [0.12803146150978031,0.570162439213454,0.7738937317070286];
triguerr = [0.6832436932180734,24.30471169006181,7.709132227291314];
quaduerr = [0.459422207138292,1.4394024790349573,0.40365298335333777];
vorperr  = [2.5605880488235013,168.80942575660904,37.048307481973225];
trigperr = [3.7109040018213126,9.637150811152049,28.448960440767276];
quadperr = [5.038100230670561,47.79150424755166,4.767992011767707];
set(gca,'FontSize',15)

%Plotting the presure
figure(1)
clf
loglog(vordx,vorperr,'s','LineWidth',3,'color','r');
hold on
loglog(quaddx,quadperr,'o','LineWidth',3,'color','k');
loglog(trigdx,trigperr,'d','Linewidth',3,'color','b');
loglog(vordx,vorperr,'LineWidth',2,'Color','r');
loglog(quaddx,quadperr,'--','LineWidth',2,'Color','k');
loglog(trigdx,trigperr,':','Linewidth',2,'Color','b');

set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold on
%Pick a basis point for the triangle
%xseed = 0.9*h(2)+h(1)*0.1;
%yseed = 0.00005;
%desiredSlope Of triangle
%slope = 2;
%Another x point
%xnext = h(1);


%b = log10(yseed/(xseed^slope));


%ynext = 10^(b)*xnext^slope;


%x = [xseed, xseed, xnext, xseed];
%y = [yseed, ynext, ynext, yseed];

%loglog(x,y,'LineWidth',2,'Color','k')

%text(xseed*0.7+0.3*xnext,ynext+0.003,'1','Fontsize',15)
%text(xseed-0.001,0.5*yseed+0.5*ynext,'2','Fontsize',15)





%title('Triangular Mesh')





%Now the ticks
%The x-axis
%leftaxis = 0.005;
%rightaxis = 0.13;
%upaxis = 0.012;
%downaxis = 10e-6;

%axis([leftaxis rightaxis downaxis upaxis])

legend('Vor','Quad', 'Trig')
legend('Location','northwest')



%Plotting The velocity Field
figure(2)
clf
loglog(vordx,voruerr,'s','LineWidth',3,'color','r');
hold on
loglog(quaddx,quaduerr,'o','LineWidth',3,'color','k');
loglog(trigdx,triguerr,'d','Linewidth',3,'color','b');
loglog(vordx,voruerr,'LineWidth',2,'Color','r');
loglog(quaddx,quaduerr,'--','LineWidth',2,'Color','k');
loglog(trigdx,triguerr,':','Linewidth',2,'Color','b');

set(gca,'FontSize',15)
set(gca,'linewidth',2)
hold on
%Pick a basis point for the triangle
%xseed = 0.9*h(2)+h(1)*0.1;
%yseed = 0.00005;
%desiredSlope Of triangle
%slope = 2;
%Another x point
%xnext = h(1);


%b = log10(yseed/(xseed^slope));


%ynext = 10^(b)*xnext^slope;


%x = [xseed, xseed, xnext, xseed];
%y = [yseed, ynext, ynext, yseed];

%loglog(x,y,'LineWidth',2,'Color','k')

%text(xseed*0.7+0.3*xnext,ynext+0.003,'1','Fontsize',15)
%text(xseed-0.001,0.5*yseed+0.5*ynext,'2','Fontsize',15)





%title('Triangular Mesh')





%Now the ticks
%The x-axis
%leftaxis = 0.005;
%rightaxis = 0.13;
%upaxis = 0.012;
%downaxis = 10e-6;

%axis([leftaxis rightaxis downaxis upaxis])

legend('Vor','Quad', 'Trig')
legend('Location','northwest')
end