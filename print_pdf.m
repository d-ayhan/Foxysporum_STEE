function print_pdf(h,filename)
% Saves the figure as seen on screen in PDF
h.Units = 'Inches';
pos = h.Position;
h.PaperPositionMode = 'Auto';
h.PaperUnits = 'Inches';
h.PaperSize = [pos(3), pos(4)];
h.Renderer = 'Painters';
print(h,filename,'-dpdf','-r0')
end