function changedWBDFcn(src,ev,flag)

    hFig = ev.AffectedObject;            %# figure handle
    currFcn = ev.NewValue;               %# current callback function
    delete(src);                         %# delete event listener
    if nargin < 3, flag = false; end     %# determine flag

    %# hijack WindowButtonDownFcn function
    set(hFig, 'WindowButtonDownFcn',{@wbdFcn,currFcn,flag})


    %# callback function
    function wbdFcn(o,e,currFcn,flag)
        %# skip anything but single-clicks
        if strcmpi(get(hFig,'SelectionType'),'open')
            return
        end

        %# evaluate previous callback function
        hgfeval(currFcn)  %# getline('FirstButtonDown'),getline('NextButtonDown')

        %# repeat process after first click
        if flag
            addlistener(handle(hFig), 'WindowButtonDownFcn', ...
                'PostSet', {@changedWBDFcn,true});
        end
    end
end