function M2 = memMapmotionCorrection(name)
gcp;
tic; memMap = read_file(name); toc; % read the file (optional, you can also pass the path in the function instead of Y)
memMap = single(memMap);   
% convert to single precision 
T = size(memMap,ndims(memMap));
memMap = memMap - min(memMap(:));
%% Downsample data
disp('Downsampling data...');
memMap = downsample_data(memMap,'spacetime',2,1); %tsub,ssub

%% now try non-rigid motion correction (also in parallel)
% disp('Starting first motion correction pass...')
% options_nonrigid = NoRMCorreSetParms('d1',size(memMap,1),'d2',size(memMap,2),'grid_size',[64,64],'bin_width',200,'max_shift',[30 30],'max_dev',[8 8],'us_fac',50,'init_batch',200);
% tic; [M1,shifts2,template2,options_nonrigid] = normcorre_batch(memMap,options_nonrigid); toc
disp('Starting second motion correction pass...')
options_nonrigid = NoRMCorreSetParms('d1',size(memMap,1),'d2',size(memMap,2),'grid_size',[64,64],'bin_width',200,'max_shift',[15 15],'max_dev',[2 2],'us_fac',50,...
    'init_batch',200,'output_type','tiff');
tic; [M2,shifts2,template2,options_nonrigid] = normcorre_batch(memMap,options_nonrigid); toc
%% compute metrics

% nnY = quantile(memMap(:),0.005);
% mmY = quantile(memMap(:),0.995);
% 
% [cY,mY,vY] = motion_metrics(memMap,10);
% [cM2,mM2,vM2] = motion_metrics(M2,10);
% T = length(cY);
% %% plot metrics
% figure;
%     ax1 = subplot(2,2,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
%     ax3 = subplot(2,2,2); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
%     linkaxes([ax1,ax3],'xy')
%% plot a movie with the results

% figure;
% for t = 1:1:T
%     subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
%     subplot(122);imagesc(M2(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
%     set(gca,'XTick',[],'YTick',[]);
%     drawnow;
%     pause(0.02);
% end