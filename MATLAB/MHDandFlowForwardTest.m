
  vordx    = [0.12803687993289598,0.06772854614785964,0.03450327796711771,0.017476749542968805,0.008787156237382746]/2;
  quaddx   = [0.16666666666666666,0.08333333333333333,0.043478260869565216,0.021739130434782608,0.010989010989010988]/2;
  trigdx    = [0.10101525445522107,0.05018856132284956,0.025141822757713456,0.012525468249897755,0.006261260829309998]/2;
  Eq = 'NewMom';
  Eq = 'OldMom';
  Eq = 'Divu';
  Eq = 'Faraday';
  Eq = 'Ampere-Ohm';
  switch Eq
      case 'NewMom'
          quadzLerr = [0.03522673023073798,0.006980113284582856,0.0020269543152644726,0.0005344146576953798,0.00014047945509407787];
          trigzLerr = [0.031221250277674455,0.007491345281056842,0.0018937665468182216,0.0004905136963706086,0.0001255133870815599];
          vorzLerr  = [0.017257647494405382,0.0039286279209871915,0.0011517394116343856,0.0003067400812248325,7.724589031578206e-05];
      
      
      figure(1)
clf
loglog(vordx,vorzLerr,'o','LineWidth',3,'color','k')
hold on
loglog(quaddx,quadzLerr,'o','LineWidth',3,'color','r')
loglog(trigdx,trigzLerr,'o','LineWidth',3,'color','b')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
loglog(vordx,vorzLerr,'LineWidth',3,'color','k')
loglog(quaddx,quadzLerr,'LineWidth',3,'color','r')
loglog(trigdx,trigzLerr,'LineWidth',3,'color','b')


      
xseed = 0.9*vordx(5)+vordx(4)*0.1;
yseed = vorzLerr(5)+0.0003;
%desiredSlope Of triangle
slope = 2;
%Another x point
xnext = vordx(4);


b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'LineWidth',2,'Color','k')
legend('vor','quad','trig')

      case 'OldMom'
           vorzLerr  = [0.18734420855049483,0.04713736352965822,0.012045765222006122,0.012045765222006122,0.0007575701884083207];
           quadzLerr  = [0.21834240463027907,0.05028271747416541,0.013262762051015418,0.0032899201518516987,0.0008391339525621642];
           trigzLerr = [0.11430722260468898,0.0314077236685773,0.008323251372461454,0.002161081460775735,0.0005164577617416136];
            figure(1)
clf
loglog(quaddx,quadzLerr,'o','LineWidth',3,'color','r')
hold on
%loglog(vordx,vorzLerr,'o','LineWidth',3,'color','k')

loglog(trigdx,trigzLerr,'o','LineWidth',3,'color','b')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
%loglog(vordx,vorzLerr,'LineWidth',3,'color','k')
loglog(quaddx,quadzLerr,'LineWidth',3,'color','r')
loglog(trigdx,trigzLerr,'LineWidth',3,'color','b')


      
xseed = 0.9*vordx(3)+vordx(2)*0.1;
yseed = vorzLerr(3)-0.004;
%desiredSlope Of triangle
slope = 2;
%Another x point
xnext = vordx(4);


b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'LineWidth',2,'Color','k')
%legend('vor','quad','trig')
legend('quad','trig')
      case 'Divu'
            vorzLerr  = [1.5682559248043064e-05,3.019296896700456e-06,3.28175828629172e-07,3.590693288629084e-08,4.365375062965227e-09];
            quadzLerr = [1.1396292117854075e-05,1.9556420624743963e-06,2.7543435160716325e-07,3.7842327913859874e-08,5.0362104989682035e-09];
            trigzLerr = [6.180669601655305e-05,8.666477019743677e-06,9.03723121580749e-07,1.053011747416689e-07,1.5653858436310615e-08];
      
      clf

loglog(quaddx,quadzLerr,'o','LineWidth',3,'color','r')
hold on
%loglog(vordx,vorzLerr,'o','LineWidth',3,'color','k')
loglog(trigdx,trigzLerr,'o','LineWidth',3,'color','b')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
%loglog(vordx,vorzLerr,'LineWidth',3,'color','k')
loglog(quaddx,quadzLerr,'LineWidth',3,'color','r')
loglog(trigdx,trigzLerr,'LineWidth',3,'color','b')

xseed = 0.9*vordx(5)+vordx(4)*0.1;
yseed = vorzLerr(5)+0.0000001;
%desiredSlope Of triangle
slope = 3;
%Another x point
xnext = vordx(4);


b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'LineWidth',2,'Color','k')
      

%legend('vor','quad','trig')
legend('quad','trig')
      case 'Faraday'
            vorzLerr    = [2.100879854860948e-08,2.701366740769061e-09,8.244676358813689e-11,2.742957387599599e-12,7.0118813932614196e-15];
            quadzLerr = [6.982114038831059e-06,1.1202065381797656e-07,2.8012020017401204e-09,3.6415928728210414e-11,7.879414053236425e-13];
            trigzLerr = [6.750578341432932e-08,1.461170835265526e-09,8.725553486004786e-12,2.1301229880981986e-15,1.0079606939290403e-14];
            clf
loglog(quaddx,quadzLerr,'o','LineWidth',3,'color','r')
hold on
%loglog(vordx,vorzLerr,'o','LineWidth',3,'color','k')

loglog(trigdx,trigzLerr,'o','LineWidth',3,'color','b')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
%loglog(vordx,vorzLerr,'LineWidth',3,'color','k')
loglog(quaddx,quadzLerr,'LineWidth',3,'color','r')
loglog(trigdx,trigzLerr,'LineWidth',3,'color','b')

xseed = 0.9*vordx(4)+vordx(3)*0.1;
yseed = vorzLerr(4)+1e-10;
%desiredSlope Of triangle
slope = 4;
%Another x point
xnext = vordx(3);


b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'LineWidth',2,'Color','k')
      

%legend('vor','quad','trig')
legend('quad','trig')
      case 'Ampere-Ohm'
          vorzLerr    = [0.013888609088342278,0.0021782575478403233,0.008300069883450101,0.0015028611794426045,0.0005373523478602525];
            quadzLerr = [0.0014317350324327912,0.0005725765264469921,0.00015395815102695474,4.0924059571475536e-05,1.0457951646922318e-05];
        trigzLerr     = [0.005265993916044133,0.0017321217902406076,0.0004654692346636819,0.00012526048693492533,3.1157973623239146e-05];
  
        clf

loglog(quaddx,quadzLerr,'o','LineWidth',3,'color','r')
hold on
%loglog(vordx,vorzLerr,'o','LineWidth',3,'color','k')
loglog(trigdx,trigzLerr,'o','LineWidth',3,'color','b')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
%loglog(vordx,vorzLerr,'LineWidth',3,'color','k')
loglog(quaddx,quadzLerr,'LineWidth',3,'color','r')
loglog(trigdx,trigzLerr,'LineWidth',3,'color','b')

xseed = 0.9*vordx(4)+vordx(3)*0.1;
yseed = vorzLerr(4)-1e-3;
%desiredSlope Of triangle
slope = 2;
%Another x point
xnext = vordx(3);


b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];

loglog(x,y,'LineWidth',2,'Color','k')
      

%legend('vor','quad','trig')
legend('quad','trig')
  end 

  %ForwardTest Flow
  %vorzLerr = [2E-27,1E-26,5E-26,4E-25,4E-24];

  %ForwardTest Flow 
  %quadzLerr = [4E-28,2E-27,1E-26,1E-25,1E-24];
 
  
  %ForwardTest Flow 
  %trigzLerr = [3E-27,5E-26,3E-25,4E-24,5E-23];