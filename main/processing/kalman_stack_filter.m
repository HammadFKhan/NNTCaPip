function Img_p= kalman_stack_filter(Img,G, V)
%     Kalman filter with each frame regarded as one measurement
%     Parameters
%     ----------
%     G : scalar
%         filter gain
%     V : scalar
%         estimated variance
%
%     Returns
%     -------
%     None.
    if nargin <2 , G = 0.8; end
    if nargin <3 , V = 0.05; end
    dimz = size(Img, 3);
    dimx = size(Img, 2);
    dimy = size(Img, 1);
    Img_s = reshape(Img, [dimx*dimy, dimz]);
    %initialization
    Ik = Img_s(:, 1); % use the first image as prediction seed
    Ek = V; % use the estimated variance
    Img_p = zeros(size(Img),'uint16');
    Img_p(:,:,1) = Img(:,:,1);
    % iteration
    for k = 1:dimz-1
        % correction
        Mk = Img_s(:, k+1); % current measurement
        K = Ek/(Ek+V);   %kalman gain
        Ik = G*Ik+(1-G)*Mk+K*(Mk-Ik);  % updated correction
        Ek = Ek*(1-K); % updated estimation of variance
        %prediction
        Img_p(:,:,k+1) = reshape(Ik, [dimx, dimy]);
    end
end







