MType = 'Trig';
%MType = 'Vor';
%MType = 'Quad';
switch MType
    case 'Vor'
        dx    = [0.12803687993289598,0.06772854614785964,0.03450327796711771,0.017476749542968805,0.008787156237382746]/2;
        zLerr = [7.242154094888598e-05,5.207216907887613e-06,3.2367843338931834e-07,1.88541451340285e-08,1.3832059941876196e-09];
        zHerr = [0.00015928822837452117,1.1307537451266114e-05,8.494186296559292e-07,5.583176321266592e-08,3.685709071277188e-09];
        xseed = 0.9*dx(5)+dx(4)*0.1;
    yseed = (10^(-7.5));
%desiredSlope Of triangle
    slope = 2;
%Another x point
    xnext = dx(4);
    
    picfile = 'VoronoiMesh.png';
    leftaxis  = 0.008/2;
    rightaxis = 0.15/2;
    upaxis = 5E-4;
    downaxis = 5E-10;
    case 'Quad'
        dx    = [0.16666666666666666,0.08333333333333333,0.043478260869565216,0.021739130434782608,0.010989010989010988]/2;
        zLerr = [0.00017135388454758704,1.0409429030744377e-05,7.751338118033857e-07,4.894997251980726e-08,3.1816456136368743e-09];
        zHerr = [0.00021717128154863374,1.2987780303674867e-05,9.459868386407777e-07,6.106628491764354e-08,3.910463064471514e-09];
        
        xseed = 0.9*dx(5)+dx(4)*0.1;
    yseed = (10^(-7.5));
%desiredSlope Of triangle
    slope = 2;
%Another x point
    xnext = dx(4);
    
    picfile = 'PerturbedSquaresMesh.png';
    leftaxis = 0.009/2;
    rightaxis = 0.18/2;
    upaxis = 5E-4;
    downaxis = 5E-10;
    case 'Trig'  
    dx = [0.10101525445522107,0.05018856132284956,0.025141822757713456,0.012525468249897755,0.006261260829309998]/2;
    zLerr = [4.968788128856261e-05,2.7879765054361627e-06,3.110158703911736e-07,2.2914630637416167e-08,9.574605552842286e-10];
    zHerr = [5.389146388123933e-05,1.3070844661289271e-05,6.762815543126521e-07,4.2433729863233793e-08,2.4944268872673092e-09];
    
    xseed = 0.9*dx(5)+dx(4)*0.1;
    yseed = (10^(-7.5));
%desiredSlope Of triangle
    slope = 2;
%Another x point
    xnext = dx(4);
    
    picfile = 'TriangularMesh.png';
    leftaxis = 0.005/2;
    rightaxis = 0.13/2;
    upaxis = 1E-4;
    downaxis = 5E-10;
end
   

figure(1)
clf
loglog(dx,zLerr,'o','LineWidth',3,'color','k')
hold on
loglog(dx,zHerr,'o','LineWidth',3,'color','r')

set(gca,'FontSize',15)
set(gca,'linewidth',2)
loglog(dx,zLerr,'LineWidth',3,'color','k')
loglog(dx,zHerr,'LineWidth',3,'color','r')




b = log10(yseed/(xseed^slope));


ynext = 10^(b)*xnext^slope;


x = [xseed, xseed, xnext, xseed];
y = [yseed, ynext, ynext, yseed];
loglog(x,y,'LineWidth',2,'Color','k')
legend('L2','H1')
axis([leftaxis rightaxis downaxis upaxis])
test = imread(picfile); 
axes('position',[0.657 0.145 0.225 0.25]); 
imagesc(test)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

%legend('Location','northwest')