
h = figure;

run fnodeX.m
run fnodeY.m

mnodeX = reshape(nodeX,[51,51]);
mnodeY = reshape(nodeY,[51,51]);
for(i=0:20)
    fvel = join({'velx',num2str(i)});
    fvel = join({fvel{1},'.m'});
    fvel = fvel{1};
    fvel = erase(fvel,' ');
    mfile = dir(fvel);
    [trash, name] = fileparts(mfile(1).name);
    run(name);
    
    fvel = join({'vely',num2str(i)});
    fvel = join({fvel{1},'.m'});
    fvel = fvel{1};
    fvel = erase(fvel,' ');
    mfile = dir(fvel);
    [trash, name] = fileparts(mfile(1).name);
    run(name);

    %velx(1:3:end)
    %nodeX(1:3:end)
    quiver(nodeX(1:8:end),nodeY(1:8:end),velx(1:8:end),vely(1:8:end))
    axis([-0.2 1.2 -0.2 1.2])
    
    
    starty = 0:0.2:1;
    startx = 0.1*ones(size(starty));
    %startx = [0.5];
    %starty = [0.5];
    velx  = reshape(velx, [51,51]);
    vely  = reshape(vely, [51,51]);
    
    streamline(mnodeX',mnodeY',velx',vely',startx,starty)
    pause(0.5)
    
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
    if i == 0 
         imwrite(imind,cm,'coupledCavitiy.gif','gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,'coupledCavitiy.gif','gif','WriteMode','append'); 
    end 
    
end


