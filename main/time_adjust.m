function time = time_adjust(num_images,frame_rate)
time = num_images/frame_rate;
time = linspace(0,time,num_images);
end
