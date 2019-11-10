function [ol] = checkforoverlap(cxes,cyes,n1x,n2x,n1y,n2y)
%Function to check for overlapping between the line segment defined by
%points (n1x, n1y) and (n2x, n2y) and connected line segments defined by
%points (cxes(i), cyes(i)) and (cxes(i+1), cyes(i+1)
%Outputs 1 if overlap exists at all, 0 if it doesn't
    
    nx = [n1x n2x]; ny = [n1y n2y];
    [nLx,ind] = min(nx); nLy = ny(ind);
    [nRx,ind] = max(nx); nRy = ny(ind); %order points from left to right
    
    len = length(cxes);
    ol = 0; %start with overlap set to 0, change to 1 if overlap found
    for j = 1:(len-1) %check over all line segments (# of points - 1)
        [cLx,ind] = min(cxes(j:j+1)); cLy = cyes(j+ind-1);
        [cRx,ind] = max(cxes(j:j+1)); cRy = cyes(j+ind-1); %order points from left to right

        
        if (nLx >= cRx) || (nRx <= cLx) %Check if there's no x overlap, then lines can't overlap
            j = j + 1;
            continue
        end
        if (cRx-cLx) == 0) 
            cRx = cRx - 0.0000000001
        end
        if (nRx-nLx) == 0) %If one is vertical then make slightly not vertical, within error but removes NAN
            nRx = nRx - 0.0000000001
        end
        cm = (cLy - cRy) / (cLx - cRx); %Calculate slopes and intercepts of associated lines
        nm = (nLy - nRy) / (nLx - nRx);
        cb = cLy - cm*cLx;
        nb = nLy - nm*nLx;
        
        if (cm == nm) && (cb == nb) %If parallel and y intersect are same (x's overlap above), there's overlap
            ol = 1; break
        end
	    
        xint = (nb - cb) / (cm - nm); % location of lines (not segments) intersecting
        if (xint > max([cLx nLx])) && (xint < max([cRx nRx])) %check if location exists in segments
            ol = 1; break
        end
        
        j = j + 1;
    end
end