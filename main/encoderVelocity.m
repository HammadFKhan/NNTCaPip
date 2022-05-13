function Velocity = encoderVelocity(VR_data)
binWidth = 4;
for trial = 1:length(VR_data.Position)
    position = abs(VR_data.Position{trial});
    edgesPos = 0:1:abs(ceil(VR_data.Position{1}(end)));
    VRposition =  discretize(position,edgesPos);
    time = VR_data.Time{trial};
    edgesTime = 0:binWidth:time(end); %Bin width on the track
    VRtimebin = discretize(time,edgesTime);
    
    % Checks binning routine
    checkBin = find(diff(VRtimebin)>1);
    count = 0;
    while checkBin>0
        for check = 1:size(checkBin,1)
            checkValue = checkBin(check)+1; %Adds 1 to the location of the skipped cell
            VRtimebin(checkValue,1) = VRtimebin(checkValue,1)-1;
        end
        checkBin = find(diff(VRtimebin)>1);
        count = count+1;
    end
    disp(['Binning routine adjusted: ' num2str(count)]);
    clear count
    % Checks for NaN values
    nanValues = find(isnan(VRtimebin)); 
    VRtimebin(isnan(VRtimebin)) = VRtimebin(nanValues(1)-1,1); % Assigns last value to NaN
    
    VRtime = VR_data.Time{trial};
    dP = diff(position);
    dP = vertcat(0,dP);
   
    % Custom moving velocity window
    for i = 1:VRtimebin(end)
        findBinT = find(VRtimebin==i);
        sumPosition = sum(dP(findBinT,1)); %#ok<FNDSB>
        Vel(i,1) = i;
        Vel(i,2) = sumPosition;
    end
    Velocity = Vel;
end