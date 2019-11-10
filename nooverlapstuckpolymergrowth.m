function [length] = nooverlapstuckpolymergrowth(na,nb,alpha,beta)
%Simulates polymer growth from two monomers with
%populations na and nb 
%bond angles alpha, beta
%This version does not allow any overlap but has no excluded areas

    ntot = na + nb;
    [~,ord] = sort(rand(ntot,1));% predetermines monomer ordering, if ord <= na, then monomer a chosen, else monomer b
    aord = ord<= na;
    angle = aord*(alpha-beta)+beta;  % determine bond angle based on ordering
    
    mxes=zeros(ntot,1); myes=mxes; % initialize arrays containing monomer locations (for plotting)
    nxes=zeros(ntot*2,1); nyes=nxes; % initialize arrays containing node + monomer locations (for line segment determination)
    
    i = 1; %find the first location
    [mx,my,n2x,n2y] = addmonomer(0,0,angle(i)); %add the first monomer
    mxes(i) = mx; myes(i) = my; %update monomer location array
    nxes(i*2-1:i*2) = [mx,n2x]; nyes(i*2-1:i*2) = [my,n2y]; %update node location array
    i = i + 1; %increase monomer iteration value
    sb = 1; %initialize how many monomers to remove when stuck (start with 1)
    imax = 1; %initialize max number of monomers placed
    while i <= ntot %iterate through, allows easy hopping back and forth for i while removing monomers
        [mx,my,n2x,n2y,ol] = addmonomernooverlap([0; nxes(1:(i-1)*2)],[0; nyes(1:(i-1)*2)],angle(i)); %add subsequent monomers
        if ol == 1 %check if overlap occured
            sb = sb + 1; %increase the number of monomers to remove (lowest is 2 which is the current unplaceable monomer and previous one)
            i = i - sb; %remove monomers
            if i <= 1 %check that we haven't started completely over, if yes then keep the first monomer
                i = 1;
                sb = 1;
            end
            continue
        end
        mxes(i) = mx; myes(i) = my; %update monomer locations
        nxes(i*2-1:i*2) = [mx,n2x]; nyes(i*2-1:i*2) = [my,n2y]; %update node and monomer locations
        i = i + 1;
        if i > imax %check if exceeded the previously set maximum chain length
            imax = i; %redefine the maximum chain length
            sb = 1; %reset step back parameter because if we were stuck in a trapped region, likely out now
        end
    end
    
    length = sqrt(nxes(ntot*2)^2+nyes(ntot*2)^2);
%     p=plot([0;nxes],[0;nyes],'-','Color','black'); %quality checks, plotting chain
%     daspect([1 1 1])
%     hold on
%     p=plot(mxes(find(ord>na)),myes(find(ord>na)),'o',...
%         'MarkerEdgeColor','red','MarkerFaceColor','red');
%     p=plot(mxes(find(ord<=na)),myes(find(ord<=na)),'o',...
%         'MarkerEdgeColor','blue','MarkerFaceColor','blue');
%     hold off
end

function [mx,my,n2x,n2y,ol] = addmonomernooverlap(nxes,nyes,theta)
%function to add a monomer at random given bond angle theta and check for
%overlap with lines formed by previous node locations nxes and nyes
    ol = 1; %initalize overlap value at 1, if not overlap it will be set to 0
    maxnit = 100; %only search for 100 locations, while we may miss some valid locations it's likely that we'd get stuck soon anyways if 100 can't find a good one
    i = 1; %intialize iteration step
    while ol == 1 && i <= maxnit %repeat searching for monomers until either found a location without overlap or reached max iterations
        ang = rand*2*pi ;%choose a random value for the angle

        mx = nxes(end)+cos(ang); my = nyes(end)+sin(ang); %find x,y locations for new monomer based on angle and previous location
  
        
        pn = round(rand)*2-1 ; %determine whether that monomer is bound on its left or right side
        if pn > 0 %depending on left or right side bound, determine location of end node based on bond angle
            n2x = mx-cos(theta-ang) ;
            n2y = my+sin(theta-ang) ;
        else
            n2x = mx+sin(theta+ang-pi/2) ;
            n2y = my-cos(theta+ang-pi/2) ;
        end
        ol = checkforoverlap(nxes,nyes,nxes(end),mx,nyes(end),my); %check if overlap occurs with first bond and all preceding bonds
        
        if ol == 1 %if overlap occurs, skip next check
            i = i + 1;
            continue
        end
        
        ol=checkforoverlap(nxes,nyes,mx,n2x,my,n2y); %check if overlap occurs with second bond and all preceding bonds (no need to worry about 1st bond of same monomer)
        i = i + 1;
    end
end