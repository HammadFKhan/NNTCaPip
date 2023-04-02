function Connected_ROI = Connectivity_dice(corr,thres)
if nargin < 3; thres = .5;end

[r,c] = find(tril(corr,-1)>thres);
if ~isempty(r)
for i = 1:length(r)
    v(i,1) = corr(r(i),c(i));
end

Connected_ROI = [r c v];
else
    Connected_ROI = [];
end

% L = 1;
% for m=1:length(ROI)
%     for j = 1:length(ROI)
%         if m>j
%             rhoVal = corr(j,m);
%             if rhoVal >= thres
%                 Connected_ROI(L,:) = [j,m,corr(j,m)];
%                 L = L + 1;
%             end
%         end
%     end
% end
% 
% if L == 1
%     Connected_ROI = [0,0,0];
% end