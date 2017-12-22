function framegeneral(hObject, imagearray)

% try
    im = get(gca, 'Children');

    frame = round(get(sliderhandle, 'Value'));

    set(sliderhandle, 'Value', frame);

    set(im, 'CData', imagearray(:,:,frame));
% catch
%     disp('set the scroll callback as an anonymouse function with inputs SLIDERHANDLE and IMAGEARRAY')
% end


