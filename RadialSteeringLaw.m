function [output, condOutput] = RadialSteeringLaw(subID)
    % Seed rng stream
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
    
    % Ensure a compatible version of PTB is installed
    AssertOpenGL;
    
    % Cache m-file for later, time-critical use use
    GetSecs();
    
    % Color constants
    BG_COLOR = [128 128 128];
    HIGH_COLOR = [255 255 255];
    LOW_COLOR = [64 64 64];
    TX_COLOR = [255 255 255];
    VERBOSE = false;
    
    % Prompt for subject if not provided
    if(~exist('subID','var'))
        subID = inputdlg('Enter subject ID','Subject Information');
        if(isempty(subID))
            disp('Aborting...');
            return;
        end
        subID = str2double(subID{1});
    end

    % # of points per condition
    numPositions = 25;
    pixelUnit = 48;
    
    % Define target sizes and distances w/r/t pixelUnit
    targetSizes = [1/2, 1, 3/2, 2];
    distances = [4,8,12,16];

    % Generate the correct chord sequence for this # of positions
    sequence = [];
    for i=1:floor(numPositions / 2)
        sequence = [sequence, i, ceil(numPositions / 2) + i];
    end
    if(mod(numPositions,2) ~= 0) % odd # of positions
        sequence = [sequence, ceil(numPositions / 2), 1];
    else
        sequence = [sequence, ceil(numPositions / 2) + 1, 1];
    end
    
    % Create randomized design matrix
    designMatrix = fullfact([length(targetSizes),length(distances)]);
    conditions = [targetSizes(designMatrix(:,1));distances(designMatrix(:,2))]';
    conditions = conditions(randperm(length(conditions)),:);

    % Initialize output matrices
    output = zeros((length(sequence)-1) * length(conditions),13);
    condOutput = zeros(length(conditions),7);
    
    try
        % Open experimental screen on last screen reported by Screens
        screen = max(Screen('Screens'));
        [windowPtr, rect] = Screen('OpenWindow',screen, BG_COLOR);
        [mx,my] = RectCenter(rect);
        
        % Enforce desired mouse cursor display
        ShowCursor('Arrow');

        % k = condition index
        for k=1:length(conditions)
            
            % Extract condition parameters
            targetWidth = pixelUnit * conditions(k,1);
            distance = pixelUnit * conditions(k,2);
            
            radius = targetWidth / 2;

            % Dynamically calculate locations of each start/target
            angles = 0:360/numPositions:360-(360/numPositions);
            angle_offset = -90;
            angles = angles + angle_offset;
            config_radius = distance / 2;

            locations = [cosd(angles); sind(angles)]' * config_radius;
            dirVector = zeros(length(sequence) - 1,2);
            normVector = zeros(length(sequence) - 1,2);
            trueDistance = zeros(length(sequence) - 1,1);

            % Calculate direction unit vectors and true distances
            for i=1:length(sequence) - 1
                a = locations(sequence(i),:);
                b = locations(sequence(i+1),:);
                
                trueDistance(i) = dist(a,b); % These should all be the same for odd
                dirVector(i,:) = (b - a) / trueDistance(i);
                normVector(i,:) = dirVector(i,:) * [0 1;-1 0];
            end
            
            % Translate locations to absolute screen coordinates
            locations(:,1) = locations(:,1) + mx;
            locations(:,2) = locations(:,2) + my;

            
            polys = cell(length(sequence) - 1,3);
            % {:,1} Envelope
            % {:,2} Start
            % {:,3} End
            
            % Calculate all polys
            for i=1:length(sequence) - 1
                targetNode = locations(sequence(i),:);
                nextNode = locations(sequence(i+1),:);
                vecNorm = normVector(i,:);
                vecDir = dirVector(i,:);
                
                polys{i,1} = [
                    targetNode - vecNorm * radius;
                    targetNode + vecNorm * radius;
                    nextNode + vecNorm * radius;
                    nextNode - vecNorm * radius];
                polys{i,2} = [
                    targetNode - vecNorm * radius;
                    targetNode - vecNorm * radius - vecDir * pixelUnit;
                    targetNode + vecNorm * radius - vecDir * pixelUnit;
                    targetNode + vecNorm * radius];
                polys{i,3} = [
                    nextNode - vecNorm * radius;
                    nextNode - vecNorm * radius + vecDir * pixelUnit;
                    nextNode + vecNorm * radius + vecDir * pixelUnit;
                    nextNode + vecNorm * radius];
            end

            % j = sequence index
            j = 1;
            while(j < length(sequence))
                [startX,startY] = GetMouse(windowPtr);

                % First trial, show instructions
                if(j == 1) 
                    DrawFormattedText(windowPtr,'Click the highlighted target to begin','center',rect(4) - 32,TX_COLOR);
                end
                
                complete = false;
                while(~complete)
                    % Render steering envelope
                    Screen('FramePoly', windowPtr, LOW_COLOR, polys{j,1}, 1);
                    % Render start area
                    Screen('FillPoly', windowPtr, HIGH_COLOR, polys{j,2}, 1);
                    % Render end area
                    Screen('FillPoly', windowPtr, LOW_COLOR, polys{j,3}, 1);
                    Screen('Flip',windowPtr);

                    % Wait for click in start area
                    waitForStart = true;
                    while(waitForStart)
                        [clicks,x,y] = GetClicks(windowPtr,0);
                        if(inpoly([x,y],polys{j,2}))
                            waitForStart = false;
                        end
                    end

                    Screen('FramePoly', windowPtr, LOW_COLOR, polys{j,1}, 1);
                    Screen('FillPoly', windowPtr, LOW_COLOR, polys{j,2}, 1);
                    Screen('FillPoly', windowPtr, HIGH_COLOR, polys{j,3}, 1);
                    Screen('Flip',windowPtr);

                    % Poll mouse continuously, start timing when mouse enters envelope,
                    % end when mouse enters end area OR reset if mouse is released 
                    % or moves out of bounds

                    enteredEnv = false;
                    disp('Starting trial');
                    while(true)
                        [x,y,buttons] = GetMouse(windowPtr);
                        
                        % Check if mb down, otherwise reset
                        if(~any(buttons))
                            Screen('FillPoly', windowPtr, LOW_COLOR, polys{j,1}, 1);
                            Screen('FillPoly', windowPtr, LOW_COLOR, polys{j,2}, 1);
                            Screen('FillPoly', windowPtr, LOW_COLOR, polys{j,3}, 1);
                            Screen('Flip',windowPtr);
                            disp('  MB released! Reset');
                            break;
                        end
                        
                        % If in envelope, check if we are in end zone
                        if(enteredEnv && inpoly([x,y],polys{j,3}))
                            disp('Entered end zone. Complete');
                            endTime = GetSecs();
                            complete = true;
                            break;
                        end
                        
                        % In envelope
                        if(inpoly([x,y],polys{j,1}))
                            % If time hasn't started, start
                            if(~enteredEnv)
                                Screen('FramePoly', windowPtr, LOW_COLOR, polys{j,2}, 1);
                                Screen('FramePoly', windowPtr, HIGH_COLOR, polys{j,1}, 1);
                                Screen('FillPoly', windowPtr, HIGH_COLOR, polys{j,3}, 1);
                                Screen('Flip',windowPtr);
                                disp('Entered envelope! Start time');
                                enteredEnv = true;
                                startTime = GetSecs();
                            end
                        elseif(enteredEnv) % Had entered but left...
                            Screen('FillPoly', windowPtr, LOW_COLOR, polys{j,1}, 1);
                            Screen('FillPoly', windowPtr, LOW_COLOR, polys{j,2}, 1);
                            Screen('FillPoly', windowPtr, LOW_COLOR, polys{j,3}, 1);
                            Screen('Flip',windowPtr);
                            disp('Diverged from envelope. Reset');
                            WaitSecs(0.5);
                            break;
                        end
                        WaitSecs(.01);
                    end
                end
                
                disp(['Time elapsed: ' num2str(endTime - startTime)]);
                
                % j = 1 begins sequence, non-informational.
%                 if(j > 1)
%                     % Calculate some values...
%                     nominalDistance = trueDistance(j - 1);
%                     nominalWidth = targetWidth;
%                     nominalID = log2((trueDistance(j - 1)/nominalWidth)+1);
%                     movementTime = clickTime - StimulusOnsetTime;
%                     actualDistance = dist([startX,startY],[endX,endY]);
%                     moveVector = [endX - startX,endY - startY];
%                     projectedDistance = dot(moveVector,dirVector(j - 1,:));
%                     
%                     % Display them
%                     if(VERBOSE)
%                         disp(['Trial ' num2str(j - 1)]);
%                         disp(['Dist     = ' num2str(nominalDistance)]);
%                         disp(['T. Width = ' num2str(nominalWidth)]);
%                         disp(['ID       = ' num2str(nominalID)]);
%                         disp(['MT       = ' num2str(movementTime) 's']);
%                         disp(['A. Dist  = ' num2str(actualDistance) 'px']);
%                         disp(['P. Dist  = ' num2str(projectedDistance) 'px']);
%                         disp('----------------------------------');
%                     end
%                     
%                     % Record them in output matrix
%                     ind = ((k - 1) * (length(sequence) - 1)) + j - 1;
%                     output(ind,1) = subID;
%                     output(ind,2) = k;
%                     output(ind,3) = ind;
%                     output(ind,4) = nominalDistance;
%                     output(ind,5) = nominalWidth;
%                     output(ind,6) = nominalID;
%                     output(ind,7) = movementTime;
%                     output(ind,8) = actualDistance;
%                     output(ind,9) = projectedDistance;
%                     output(ind,10) = startX;
%                     output(ind,11) = startY;
%                     output(ind,12) = endX;
%                     output(ind,13) = endY;
%                 elseif(remedialTrial)
%                     remedialTrial = false;
%                 end % End of non-first target
%                 
                j = j + 1;
            end % End of sequence
            
%             % Calculate start and end indices for this condition only
%             startIndex = (k - 1) * (length(sequence) - 1) + 1;
%             endIndex = k * (length(sequence) - 1);
%             condOutput(k,1) = subID;
%             condOutput(k,2) = k;
%             condOutput(k,3) = mean(output(startIndex:endIndex,7));
%             condOutput(k,4) = mean(output(startIndex:endIndex,9));
%             condOutput(k,5) = std(output(startIndex:endIndex,9)) * 4.133;
%             condOutput(k,6) = log2((condOutput(k,4) / condOutput(k,5)) + 1);
%             condOutput(k,7) = condOutput(k,6) / condOutput(k,3);
% 
%             if(VERBOSE)
%                 disp('Condition Summary');
%                 disp(['D (A)     = ' num2str(distance)]);
%                 disp(['W         = ' num2str(targetWidth)]);
%                 disp(['ID        = ' num2str(log2((distance/targetWidth) + 1))]);
%                 disp(['IDe       = ' num2str(condOutput(k,6))]);
%                 disp(['MT        = ' num2str(condOutput(k,3))]);
%                 disp(['TPiso p/e = ' num2str(condOutput(k,7))]);
%                 disp('==================================');
%             end
            % Post-condition screen
            DrawFormattedText(windowPtr,'Press any key to continue','center','center',TX_COLOR);
            Screen('Flip',windowPtr);
            KbWait([],3);
        end % End of condition
    catch E
        % On error
        Screen('CloseAll');
        disp(E.message);
        disp('You are in debug mode.  Type ''return'' to exit');
        disp('*** DATA HAS NOT BEEN SAVED: TYPE ''return'' TO SAVE ***');
        keyboard;
    end

    Screen('CloseAll');

    % Calculate some subject values
    throughput = mean(condOutput(:,7));
    [B,BINT,R,RINT,STATS] = regress(condOutput(:,3), [ones(size(condOutput(:,3),1)) condOutput(:,6)]);
    
    % Display
    disp('Subject Summary');
    disp(['TPiso      = ' num2str(throughput) ' bits/sec']);
    disp('Regression');
    disp(['        MT = ' num2str(B(1)) ' + ' num2str(B(2)) ' * IDe']);
    disp(['       R^2 = ' num2str(STATS(1))]);
    disp(['   F(1,' num2str(size(condOutput,3) - 1) ') = ' num2str(STATS(2))]);
    disp(['         p = ' num2str(STATS(3))]);
    disp(['TPreg      = ' num2str(1 / B(2)) ' bits/sec']);
    
    plot(condOutput(:,6),condOutput(:,3),'b.');
    hold on;
    plot([0 max(condOutput(:,6))],[B(1), B(1) + B(2) * max(condOutput(:,6))],'r:');
    
    % Record
    dirname = './data';
    file_sep = '_';
    file_base = num2str(subID,'%03d');
    file_suffix_output = 'trials';
    file_suffix_condOutput = 'conditions';
    file_ext = '.csv';
    
    file_output_norm = [dirname '/' file_base file_sep file_suffix_output];
    file_condOutput_norm = [dirname '/' file_base file_sep file_suffix_condOutput];
    
    file_output = file_output_norm;
    file_condOutput = file_condOutput_norm;
    attempt = 0;
    while(exist([file_output file_ext],'file') || exist([file_condOutput file_ext],'file'))
        file_output = [file_output_norm file_sep num2str(attempt,'%03d')];
        file_condOutput = [file_condOutput_norm file_sep num2str(attempt,'%03d')];
        attempt = attempt + 1;
    end
    
    csvwrite([file_output file_ext],output);
    csvwrite([file_condOutput file_ext],condOutput);
    
    % Return output and condOutput
end

function dist = dist(u,v) % 2d euclidean distance between u and v
    dist = sqrt(sum((u - v) .^ 2));
end

function dotprod = dot(u,v) % Dot product of u and v (2d only)
    dotprod = u(1) * v(1) + u(2) * v(2);
end

% function projVector = project(u,v) % Project u onto v
%     projVector = v * dot(u,v);
% end