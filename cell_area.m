function [Area] = cell_area(Rc,Rm)
Area=2*pi*Rc.^2.*(1+(1-Rm.^2*Rc^-2).^0.5);
end