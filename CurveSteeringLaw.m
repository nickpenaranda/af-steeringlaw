% TODO: Render arc

function [output, condOutput] = CurveSteeringLaw(subID)
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
    numPositions = 12;
    pixelUnit = 32;
    
    % Define experiment params w/r/t pixelUnit
    pathWidths = [1, 3/2, 2];
    pathRadii = [4,8,12];
    
    minRadius = min(pathRadii);
    maxDist = 2 * pi * minRadius * (3/4);
    pathDists = [maxDist, maxDist / 2, maxDist / 4];
    
    angles = 0:360/numPositions:360-(360/numPositions);
    angle_offset = -90;
    angles = angles + angle_offset;

    % x = r * cos(theta)
    % y = r * sin(theta)
    % arcLength = theta * r
    % theta = arcLength / r
    
    % Generate the correct chord sequence for this # of positions
%     sequence = [];
%     for i=1:floor(numPositions / 2)
%         sequence = [sequence, i, ceil(numPositions / 2) + i];
%     end
%     if(mod(numPositions,2) ~= 0) % odd # of positions
%         sequence = [sequence, ceil(numPositions / 2), 1];
%     else
%         sequence = [sequence, ceil(numPositions / 2) + 1, 1];
%     end
    
    % Create randomized design matrix
    designMatrix = fullfact([length(pathWidths),length(pathRadii),length(pathDists)]);
    conditions = [
        pathWidths(designMatrix(:,1));
        pathRadii(designMatrix(:,2));
        pathDists(designMatrix(:,3))
        ]';
        
    conditions = conditions(randperm(length(conditions)),:);

    % Initialize output matrices
    %output = zeros((length(sequence)-1) * length(conditions),13); % FIXME
    %condOutput = zeros(length(conditions),7); %FIXME
    
    try
        % Open experimental screen on last screen reported by Screens
        Screen('Preference', 'SkipSyncTests', 1);
        screen = max(Screen('Screens'));
        [windowPtr, rect] = Screen('OpenWindow',screen, BG_COLOR);
        [mx,my] = RectCenter(rect);
        
        % Enforce desired mouse cursor display
        ShowCursor('Arrow');

        % k = condition index
        for k=1:length(conditions)
            
            % Extract condition parameters
            pathWidth = pixelUnit * conditions(k,1);
            pathRadius = pixelUnit * conditions(k,2);
            distance = pixelUnit * conditions(k,3);
            trialAngle = (distance / pathRadius) * 180 / pi;
            posDelta = 360 / numPositions;
            
            % Produce sequence for this condition
            sequence = [];
            for i=0:numPositions-1
                sequence = [sequence,i+1,mod(i+ceil(trialAngle / posDelta),numPositions)+1];
            end
            
            outerStarts = [cosd(angles); sind(angles)]' * (pathRadius + (pathWidth / 2));
            innerStarts = [cosd(angles); sind(angles)]' * (pathRadius - (pathWidth / 2));
            outerCWEnds = [cosd(angles + trialAngle); sind(angles + trialAngle)]' ...
                * (pathRadius + (pathWidth / 2));
            innerCWEnds = [cosd(angles + trialAngle); sind(angles + trialAngle)]' ...
                * (pathRadius - (pathWidth / 2));
            outerCCWEnds = [cosd(angles - trialAngle); sind(angles - trialAngle)]' ...
                * (pathRadius + (pathWidth / 2));
            innerCCWEnds = [cosd(angles - trialAngle); sind(angles - trialAngle)]' ...
                * (pathRadius - (pathWidth / 2));
            
            % Translate node locations to absolute screen coordinates
            outerStarts(:,1)  = outerStarts(:,1) + mx;
            innerStarts(:,1)  = innerStarts(:,1) + mx;
            outerStarts(:,2)  = outerStarts(:,2) + my;
            innerStarts(:,2)  = innerStarts(:,2) + my;
            outerCWEnds(:,1)  = outerCWEnds(:,1) + mx;
            innerCWEnds(:,1)  = innerCWEnds(:,1) + mx;
            outerCWEnds(:,2)  = outerCWEnds(:,2) + my;
            innerCWEnds(:,2)  = innerCWEnds(:,2) + my;
            outerCCWEnds(:,1)  = outerCCWEnds(:,1) + mx;
            innerCCWEnds(:,1)  = innerCCWEnds(:,1) + mx;
            outerCCWEnds(:,2)  = outerCCWEnds(:,2) + my;
            innerCCWEnds(:,2)  = innerCCWEnds(:,2) + my;
            
            % j = sequence index
            j = 1;
            while(j < length(sequence))
                % startVec is always CW unit vector orthogonal to radius
                startVec = [cosd(angles(sequence(j)) + 90), sind(angles(sequence(j)) + 90)];
                startAngle = angles(sequence(j));
                
                outerStart = outerStarts(sequence(j),:);
                innerStart = innerStarts(sequence(j),:);
                
                % Generate start/end zone information
                if(mod(j,2)) % Odd, movement is CW
                    polyStart = [
                        outerStart;
                        innerStart;
                        innerStart - (startVec * pixelUnit);
                        outerStart - (startVec * pixelUnit)];
                    
                    outerEnd = outerCWEnds(sequence(j),:);
                    innerEnd = innerCWEnds(sequence(j),:);
                    
                    endVec = [cosd(angles(sequence(j)) + 90 + trialAngle), ... 
                              sind(angles(sequence(j)) + 90 + trialAngle)];
                    endAngle = angles(sequence(j)) + trialAngle;
                    polyEnd = [
                        innerEnd;
                        outerEnd;
                        outerEnd + (endVec * pixelUnit);
                        innerEnd + (endVec * pixelUnit)];
                    fixedStartAngle = startAngle;
                else % Even, movement is CCW
                    polyStart = [
                        innerStart;
                        outerStart;
                        outerStart + (startVec * pixelUnit);
                        innerStart + (startVec * pixelUnit)];

                    outerEnd = outerCCWEnds(sequence(j),:);
                    innerEnd = innerCCWEnds(sequence(j),:);
                    
                    endVec = [cosd(angles(sequence(j)) + 90 - trialAngle), ... 
                              sind(angles(sequence(j)) + 90 - trialAngle)];
                    endAngle = angles(sequence(j)) - trialAngle;
                    polyEnd = [
                        outerEnd;
                        innerEnd;
                        innerEnd - (endVec * pixelUnit);
                        outerEnd - (endVec * pixelUnit)];
                    fixedStartAngle = endAngle;
                end
                    

                % First trial, show instructions
                if(j == 1) 
                    DrawFormattedText(windowPtr,'Click the highlighted target to begin','center',rect(4) - 32,TX_COLOR);
                end
                
                
                complete = false;
                while(~complete)
                    % Render steering envelope
                    Screen('DrawArc', windowPtr, LOW_COLOR, ...
                        [mx - pathRadius - (pathWidth / 2), ...
                         my - pathRadius - (pathWidth / 2), ...
                         mx + pathRadius + (pathWidth / 2), ...
                         my + pathRadius + (pathWidth / 2)], ...
                        fixedStartAngle + 90, trialAngle);
                    Screen('DrawArc', windowPtr, LOW_COLOR, ...
                        [mx - pathRadius + (pathWidth / 2), ...
                         my - pathRadius + (pathWidth / 2), ...
                         mx + pathRadius - (pathWidth / 2), ...
                         my + pathRadius - (pathWidth / 2)], ...
                        fixedStartAngle + 90, trialAngle);
                    
                    % Render start area
                    Screen('FillPoly', windowPtr, HIGH_COLOR, polyStart, 1);
                    % Render end area
                    Screen('FillPoly', windowPtr, LOW_COLOR, polyEnd, 1);
                    
                    Screen('Flip',windowPtr);

                    % Wait for click in start area
                    waitForStart = true;
                    while(waitForStart)
                        [clicks,x,y] = GetClicks(windowPtr,0);
                        if(inpoly([x,y],polyStart))
                            waitForStart = false;
                        end
                    end

                    Screen('DrawArc', windowPtr, LOW_COLOR, ...
                        [mx - pathRadius - (pathWidth / 2), ...
                         my - pathRadius - (pathWidth / 2), ...
                         mx + pathRadius + (pathWidth / 2), ...
                         my + pathRadius + (pathWidth / 2)], ...
                        fixedStartAngle + 90, trialAngle);
                    Screen('DrawArc', windowPtr, LOW_COLOR, ...
                        [mx - pathRadius + (pathWidth / 2), ...
                         my - pathRadius + (pathWidth / 2), ...
                         mx + pathRadius - (pathWidth / 2), ...
                         my + pathRadius - (pathWidth / 2)], ...
                        fixedStartAngle + 90, trialAngle);

                    Screen('FillPoly', windowPtr, LOW_COLOR, polyStart, 1);
                    Screen('FillPoly', windowPtr, HIGH_COLOR, polyEnd, 1);
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
                            %Screen('FillPoly', windowPtr, LOW_COLOR, polys{j,1}, 1);
                            Screen('FillPoly', windowPtr, LOW_COLOR, polyStart, 1);
                            Screen('FillPoly', windowPtr, LOW_COLOR, polyEnd, 1);
                            Screen('Flip',windowPtr);
                            disp('  MB released! Reset');
                            break;
                        end
                        
                        % If in envelope, check if we are in end zone
                        if(enteredEnv && inpoly([x,y],polyEnd))
                            disp('Entered end zone. Complete');
                            endTime = GetSecs();
                            complete = true;
                            break;
                        end
                        
                        % In envelope
                        if(inarc(x,y,mx,my, ...
                                pathRadius - (pathWidth / 2), ...
                                pathRadius + (pathWidth / 2), ...
                                fixedStartAngle, trialAngle))
                            % If time hasn't started, start
                            if(~enteredEnv)
                                Screen('FramePoly', windowPtr, LOW_COLOR, polyStart, 1);
                                Screen('DrawArc', windowPtr, HIGH_COLOR, ...
                                    [mx - pathRadius - (pathWidth / 2), ...
                                     my - pathRadius - (pathWidth / 2), ...
                                     mx + pathRadius + (pathWidth / 2), ...
                                     my + pathRadius + (pathWidth / 2)], ...
                                    fixedStartAngle + 90, trialAngle);
                                Screen('DrawArc', windowPtr, HIGH_COLOR, ...
                                    [mx - pathRadius + (pathWidth / 2), ...
                                     my - pathRadius + (pathWidth / 2), ...
                                     mx + pathRadius - (pathWidth / 2), ...
                                     my + pathRadius - (pathWidth / 2)], ...
                                    fixedStartAngle + 90, trialAngle);
                                Screen('FillPoly', windowPtr, HIGH_COLOR, polyEnd, 1);
                                Screen('Flip',windowPtr);
                                disp('Entered envelope! Start time');
                                enteredEnv = true;
                                startTime = GetSecs();
                            end
                        elseif(enteredEnv) % Had entered but left...
                            %Screen('FillPoly', windowPtr, LOW_COLOR, polys{j,1}, 1);
                            Screen('FillPoly', windowPtr, LOW_COLOR, polyStart, 1);
                            Screen('FillPoly', windowPtr, LOW_COLOR, polyEnd, 1);
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
    dirname = './exp4data';
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