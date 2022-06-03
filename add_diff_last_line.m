function [y_target] = add_diff_last_line(y2,y1,x2,x1,x_traget)
    %finds the last line of an array in the direction dim (dim=1 columns dim=2 rows)
    % y2 is value in the last line y1 is the value one row befor it
    % x2 is position in the last line x1 is the position one row befor it
    % x_target is the postion we are looking for 
    % all varibels are vectors

    m=(y2-y1)./(x2-x1); %slop of linear fit
    b=(y1.*x2-y2.*x1)./(x2-x1); %free term in linear line
    y_target=x_traget.*m+b;
end

