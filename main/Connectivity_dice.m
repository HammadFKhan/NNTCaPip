function Connected_ROI = Connectivity_dice(corr, ROI)


L = 1;
for m=1:length(ROI)
    for j = 1:length(ROI)
        if m>j
            rhoVal = corr(j,m);
            if rhoVal >= 0.15
                Connected_ROI(L,:) = [j,m,corr(j,m)];
                L = L + 1;
            end
        end
    end
end

if L == 1
    Connected_ROI = [0,0,0];
end