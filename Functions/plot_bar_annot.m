% Adopted by Jonas Thiele
% Created in 24.06.2021
% Purpose: Plot bar plot with annotations
function plot_bar_annot(data,labels)

%%%
% inputs: data = data for bar plot
%         labels = labels for each bar
%%%



    %figure; % create new figure
    hbar = bar(data);    % create bar plot
    % get the data for all the bars that were plotted
    x = get(hbar,'XData');
    y = get(hbar,'YData');
    ygap = 0.025;  % specify vertical gap between the bar and label
    ylimits = get(gca,'YLim');
    set(gca,'YLim',[ylimits(1),ylimits(2)+0.2*max(y)]); % increase y limit for label
    
    for i = 1:length(x) % loop over each bar 
            xpos = x(i); % set x position for the text label
            if data(i)>=0
                ypos = y(i) + ygap; % set y position, including gap
            else
                ypos = y(i) - ygap;
            end
            htext = text(xpos,ypos,labels{i});          % add text label
            set(htext,'VerticalAlignment','bottom',...  % adjust properties
                      'HorizontalAlignment','center', 'fontsize',18)
    end
    
end