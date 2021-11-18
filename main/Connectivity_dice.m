function Connected_ROI = Connectivity_dice(corr, ROI,thres)
if nargin < 3; thres = .5;end

L = 1;
for m=1:length(ROI)
    for j = 1:length(ROI)
        if m>j
            rhoVal = corr(j,m);
            if rhoVal >= thres
                Connected_ROI(L,:) = [j,m,corr(j,m)];
                L = L + 1;
            end
        end
    end
end

if L == 1
    Connected_ROI = [0,0,0];
end