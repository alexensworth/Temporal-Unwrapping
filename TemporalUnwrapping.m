%{
    input:  Four dimensional, spatially unwrapped data, with the echoes in
            the 4th dimension. May or may not have temporal wraps.

            A mask, where the first 3 dimensions must be the same size as  
            the spatially unwrapped data, 4th dimension (echoes) not 
            necessary.

    output: Four dimensional, spatially and temporally unwrapped data
%}


function [ TempUnwrappedPhase] = TemporalUnwrapping(SpatUnwrappedPhase, Mask)

%Determine the number of echoes:
numechoes = size(SpatUnwrappedPhase,4);

%Allocate space for future matrices to decrease computation time:
TheMedians=zeros(1,numechoes);
diff = zeros(1,numechoes-1);

% Use the median of an entire echo as the metric

for i=1:numechoes
    TheMedians(i)=median(nonzeros(SpatUnwrappedPhase(:,:,:,i)));
end


% Determine the difference between median values for each echo

for j = 1:(numechoes-1) 
    diff(j)=TheMedians(j+1)-TheMedians(j);
end



% Systematically go through the echoes to determine where a wrap has 
% occured (triggered by a difference between echoes greater than pi)

jj=1;
while jj<numechoes
    
    %if the difference is large enough, this indicates a temporal wrap
    
    if diff(jj) > pi 
        
        % Remove 2pi from each point, multiply by the mask, write the data 
        % in this slice on the existing dataset
        SpatUnwrappedPhase(:,:,:,jj+1)=(SpatUnwrappedPhase(:,:,:,jj+1)-(2*pi)).*Mask(:,:,:); 
        
        % Re-evaluate the Medians again to ensure the wrap has been solved
        
        for i=1:numechoes
            TheMedians(i) = median(nonzeros(SpatUnwrappedPhase(:,:,:,i))); 
        end
        for j = 1:(numechoes-1)
            diff(j)=TheMedians(j+1)-TheMedians(j);
        end
        
    elseif diff(jj) < -pi 
        
        % The same as above, but for a negative phase difference - add 2pi
        SpatUnwrappedPhase(:,:,:,jj+1)=(SpatUnwrappedPhase(:,:,:,jj+1)+(2*pi)).*Mask(:,:,:);
        
        % Need to re-evaluate the Medians again and repeat.
        
        for i=1:numechoes
            TheMedians(i) = median(nonzeros(SpatUnwrappedPhase(:,:,:,i)));
        end
        for j = 1:(numechoes-1)
            diff(j)=TheMedians(j+1)-TheMedians(j);
        end
    end
    
    % If the new wrap is not solved, meaning a large phase difference,
    % repeat the process for the same echo.
    
    if diff(jj) > pi || diff(jj) < -pi
        jj=jj-1;
    end
    
    % Go to the next echo, will exit the while loop once this exceeds the
    % number of echoes present.
    jj=jj+1; 
    
end

% The once only spacially unwrapped phase is now temporally unwrapped too.

TempUnwrappedPhase=SpatUnwrappedPhase;

end
