lfplowfreq = ;
lfphighfreq = ;
spikelowfreq = ;
spikehighfreq = ;
stdmult = ;
detthresh = ;
snrthresh = ;
window = ;
windowpre= ;
windowpost= ;
windowL=max(windowpre,windowpost); %%%% the larger side of window.
windowT=windowpre+windowpost;  %%%%% stand for total window.
numchan = numel(Wideband_data);
endsample = min(cellfun(@numel, Wideband_data));

% Calculate filter coefficients
[blfp, alfp] = butter(1, [lfplowfreq lfphighfreq]/(samprate*0.5));
[b, a] = butter(1, [spikelowfreq spikehighfreq]/(samprate*0.5));
 % Determine timing for data
sizewindow=ceil(window/windowT*samprate);
spiketimeaxis = 1000*(0:1:(windowT*sizewindow-1))/samprate;
datafilter = cell(1,numchan);
datafilterlfp = cell(1,numchan);
    for channel = 1:numel(Wideband_data)
        datafilter{channel} = filter(b,a,double(Wideband_data{channel}(:, 1:endsample)), [],2);
        datafilterlfp{channel} = decimate(filter(blfp,alfp,double(Wideband_data{channel}(:, 1:endsample)), [],2), 10);
    end
   
for channel = 1:numchan
	datafilterch=datafilter{channel};
    datafilterlfpch=datafilterlfp{channel};
	ppnoiseraw(channel) =stdmult*std(datafilterch);
    nthresh = mean(datafilterch) - detthresh*std(datafilterch);
    holdermin = find(datafilterch < nthresh);
    
    chantemp = [];
    chantempholder = [];
    placeholder = 0;
    
       for i = 1:numel(holdermin)
            % Is it too early?
            if holdermin(i) <= windowL*sizewindow+1
                % And is it too close to the end?
            elseif holdermin(i) <= (numel(datafilterch) - windowL*sizewindow)
                % Figure out minimum value and point at which it occurs in
                % snippet window
                [mintemp, mintempposition] = min(datafilterch(holdermin(i)-windowpre*sizewindow:holdermin(i)+windowpost*sizewindow-1));
                % Adjust minimum position to be the index for the whole
                % array, not just the small snippet fed into it.
                mintempposition = holdermin(i) - windowpre*sizewindow-1 + mintempposition;
                % If this minimum happens to be the holder value in all
                % minimums, then continue
                if (mintemp == datafilterch(holdermin(i)));
                    placeholder = placeholder +1;
                    chantemp(placeholder, :) =  datafilterch(mintempposition-windowpre*sizewindow:mintempposition+windowpost*sizewindow-1);
                    chantempholder(size(chantempholder,2)+1:size(chantempholder,2)+windowT*sizewindow) = (mintempposition-windowpre*sizewindow):(mintempposition+windowpost*sizewindow-1);
                    sorts.timestamp{channel}(placeholder) = (mintempposition+1)/samprate;
                end
            end
       end
        
        % Verify that candidates were found, if so, save in snip
        if size(chantemp,1)>0
            sorts.snip{channel}(1:size(chantemp,1), :) = chantemp;
            chanminholder{channel}(1:length(chantempholder)) = chantempholder;
        else
            % Setup conditions for set diff below for situations with no
            % spikes
            sorts.snip{channel}(1, 1:2*sizewindow+1) = zeros(1,2*sizewindow+1);
            chanminholder{channel}(1) = 0;
        end
        clear chantemp
        %subtracts out signal and calculates the noise for each channel
        indextemp = 1:numel(datafilterch);
        % Grab all non-spike data
        nosignal4 = datafilterch(setdiff(indextemp,chanminholder{channel}(chanminholder{channel} > 0)));
        % Calculate noise on the channel
        ppnoise(channel) =  stdmult*std(nosignal4);
end