function Connected_ROI = Connectivity(PVAL, ROI)


L = 1;
for m=1:length(ROI)
    for j = 1:length(ROI)
        if m>j
            Pval = PVAL(j,m);
            if Pval < 0.01
                Connected_ROI(L,:) = [j,m,PVAL(j,m)];
                L = L + 1;
            end
        end
    end
end

if L == 1
    Connected_ROI = [0,0,1];
end




        

