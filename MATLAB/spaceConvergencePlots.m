vordx      = [0.209406,0.0939872,0.0452792,0.0225296];
vorElecErr = [0.00106894,0.000293946,8.95942e-5,2.2181e-05];
vorMagnErr = [0.0269721,0.0138868,0.00657099,0.0033001];
vorVeloErr = [0.000633585,0.000150453,3.45154e-5,8.5452e-06];
vorPresErr = [0.00714085,0.0020164,0.000429086,1.0623e-04];

quaddx      = [0.202031,0.101015,0.0505076,0.0250250];
quadElecErr = [0.00183731,0.000459349,0.000114851,2.8195e-05];
quadMagnErr = [0.0354471,0.0178215,0.00892299,0.00450001];
quadVeloErr = [0.00076285,0.000190922,4.77439e-05,1.1721e-05];
quadPresErr = [0.0080669,0.00203498,0.000509925,1.2518e-04];

hexadx      = [0.218956,0.103342,0.0542429,0.0271001];
hexaElecErr = [0.00257881,0.000660103,0.00019245,4.8037e-05];
hexaMagnErr = [0.0343722,0.0156144,0.00803157,0.0040012];
hexaVeloErr = [0.000933556,0.000237048,6.51147e-5,1.6253e-05];
hexaPresErr = [0.013759,0.00364099,0.00104277,2.6028e-04];

%Electric field plot
figure(1)
clf
loglog(vordx,vorElecErr,'o','LineWidth',3,'color','k')
hold on
loglog(quaddx,quadElecErr,'o','LineWidth',3,'color','r')
loglog(hexadx,hexaElecErr,'o','LineWidth',3,'color','b')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
loglog(vordx,vorElecErr,'LineWidth',3,'color','k')
loglog(quaddx,quadElecErr,'LineWidth',3,'color','r')
loglog(hexadx,hexaElecErr,'LineWidth',3,'color','b')


      
xseed = 0.8*vordx(4)+vordx(3)*0.2;
yseed = vorElecErr(4)+0.00006;
%desiredSlope Of triangle
slope = 2;
%Another x point
xnext = vordx(3);


b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'LineWidth',2,'Color','k')

set(gca,'XTick',10.^(-1))
leftaxis  = vordx(4)*0.95;
rightaxis = vordx(1)*1.1;
downaxis  = vorElecErr(4)*0.9;
upaxis    = vorElecErr(1)*2.6;
axis([leftaxis rightaxis downaxis upaxis])
legend('vor','quad','hexa')

%Magnetic field plot
figure(2)
clf
loglog(vordx,vorMagnErr,'o','LineWidth',3,'color','k')
hold on
loglog(quaddx,quadMagnErr,'o','LineWidth',3,'color','r')
loglog(hexadx,hexaMagnErr,'o','LineWidth',3,'color','b')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
loglog(vordx,vorMagnErr,'LineWidth',3,'color','k')
loglog(quaddx,quadMagnErr,'LineWidth',3,'color','r')
loglog(hexadx,hexaMagnErr,'LineWidth',3,'color','b')


      
xseed = 0.8*vordx(4)+vordx(3)*0.2;
yseed = vorMagnErr(4)*1.7;
%desiredSlope Of triangle
slope = 1;
%Another x point
xnext = vordx(3);


b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'LineWidth',2,'Color','k')

set(gca,'XTick',10.^(-1))
leftaxis  = vordx(4)*0.95;
rightaxis = vordx(1)*1.1;
downaxis  = vorMagnErr(4)*0.9;
upaxis    = vorMagnErr(1)*1.5;
axis([leftaxis rightaxis downaxis upaxis])
legend('vor','quad','hexa')

%Velocity plot
figure(3)
clf
loglog(vordx,vorVeloErr,'o','LineWidth',3,'color','k')
hold on
loglog(quaddx,quadVeloErr,'o','LineWidth',3,'color','r')
loglog(hexadx,hexaVeloErr,'o','LineWidth',3,'color','b')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
loglog(vordx,vorVeloErr,'LineWidth',3,'color','k')
loglog(quaddx,quadVeloErr,'LineWidth',3,'color','r')
loglog(hexadx,hexaVeloErr,'LineWidth',3,'color','b')


      
xseed = 0.8*vordx(4)+vordx(3)*0.2;
yseed = vorVeloErr(4)*2.6;
%desiredSlope Of triangle
slope = 2;
%Another x point
xnext = vordx(3);


b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'LineWidth',2,'Color','k')

set(gca,'XTick',10.^(-1))
leftaxis  = vordx(4)*0.95;
rightaxis = vordx(1)*1.1;
downaxis  = vorVeloErr(4)*0.9;
upaxis    = vorVeloErr(1)*1.9;
axis([leftaxis rightaxis downaxis upaxis])
legend('vor','quad','hexa')

%Pressure plot
figure(4)
clf
loglog(vordx,vorPresErr,'o','LineWidth',3,'color','k')
hold on
loglog(quaddx,quadPresErr,'o','LineWidth',3,'color','r')
loglog(hexadx,hexaPresErr,'o','LineWidth',3,'color','b')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
loglog(vordx,vorPresErr,'LineWidth',3,'color','k')
loglog(quaddx,quadPresErr,'LineWidth',3,'color','r')
loglog(hexadx,hexaPresErr,'LineWidth',3,'color','b')


      
xseed = 0.8*vordx(4)+vordx(3)*0.2;
yseed = vorPresErr(4)*3.3;
%desiredSlope Of triangle
slope = 2;
%Another x point
xnext = vordx(3);


b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'LineWidth',2,'Color','k')

set(gca,'XTick',10.^(-1))
leftaxis  = vordx(4)*0.95;
rightaxis = vordx(1)*1.1;
downaxis  = vorPresErr(4)*0.9;
upaxis    = vorPresErr(1)*2.2;
axis([leftaxis rightaxis downaxis upaxis])
legend('vor','quad','hexa')
