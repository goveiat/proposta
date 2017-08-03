function [] = exportPDF(fig, alg, ref, nmPlot)
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    nmFig = ['figuras/' alg '_' nmPlot '_' num2str(ref) '.pdf'];
    print(fig,nmFig,'-dpdf','-r0')
end