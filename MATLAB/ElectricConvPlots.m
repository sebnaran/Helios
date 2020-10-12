
    vordx    = [0.12803687993289598,0.06772854614785964,0.03450327796711771]/2;
    vorElecErr = [0.02410421435,0.005845916390,0.0007615];
    vorMagnErr = [0.0284955612,0.0146665194244,0.0073325];

    quaddx      = [0.16666666666666666,0.08333333333333333,0.043478260869565216]/2;
    quadElecErr = [0.0161072,0.00408560338,0.002365762389];
    quadMagnErr = [0.0046766442284,0.0031001,0.0023657623];


trigdx      = [0.101015,0.051886,0.0251418,0.0125255]/2;
trigElecErr = [0.011660836811078193,0.0037712742896791342,0.0008068230443378578,0.0002024975599295222];
trigMagnErr = [0.013888787706116087,0.008368253351531236,0.004707377162019469,0.0024400299532172714];


    
    


figure(1)
clf
loglog(vordx,vorElecErr,'o','LineWidth',3,'color','k')
hold on
loglog(quaddx,quadElecErr,'o','LineWidth',3,'color','b')
loglog(trigdx,trigElecErr,'o','LineWidth',3,'color','r')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
loglog(vordx,vorElecErr,'LineWidth',3,'color','k')
loglog(quaddx,quadElecErr,'LineWidth',3,'color','b')
loglog(trigdx,trigElecErr,'LineWidth',3,'color','r')


xseed = 9E-3;
yseed = 1E-3;
%desiredSlope Of triangle
    slope = 2;
%Another x point
    xnext = 3E-2;

b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];


loglog(x,y,'LineWidth',2,'Color','k')
legend('vor','quad','trig')
leftaxis = 0.009/2;
    rightaxis = 0.17/2;
    upaxis = 1E-1;
    downaxis = 1E-4;
axis([leftaxis rightaxis downaxis upaxis])


figure(2)
clf
loglog(vordx,vorMagnErr,'o','LineWidth',3,'color','k')
hold on
loglog(quaddx,quadMagnErr,'o','LineWidth',3,'color','b')
loglog(trigdx,trigMagnErr,'o','LineWidth',3,'color','r')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
loglog(vordx,vorMagnErr,'LineWidth',3,'color','k')
loglog(quaddx,quadMagnErr,'LineWidth',3,'color','b')
loglog(trigdx,trigMagnErr,'LineWidth',3,'color','r')


xseed = 9E-3;
yseed = 1E-2;
%desiredSlope Of triangle
    slope = 1;
%Another x point
    xnext = 3E-2;

b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];


loglog(x,y,'LineWidth',2,'Color','k')
legend('vor','quad','trig')
leftaxis = 0.009/2;
    rightaxis = 0.17/2;
    upaxis = 1E-1;
    downaxis = 1E-4;
axis([leftaxis rightaxis downaxis upaxis])
