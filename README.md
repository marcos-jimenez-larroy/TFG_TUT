# TFG_TUT

He añadido dos carpetas. RooResults y figures.

En mi ui (/home/marcosj/HGam/ntuplereader) están los scripts para correr sobre las muestras (AnalyzeMxAOD.C sobre el MC de señal, Analyze_data.C sobre datos y AnalyzeBkg.C sobre MC background). Me falta modificarlo para hacer el mapping que vimos todavía.

En la carpeta RootResults están los histogramas de datos, MC señal y MC background después de seleccionar los eventos y el MC de background inclusivo repesado por los eventos de yj.

En la carpeta de figures hay varios archivos que solo tienen funciones que uso luego (ATLAsUtils.C y demás). Los 'importantes' son MakeFit.C, que hace el Fit del MC de la señal con una DSCB; macro_figure_ptyy.C, donde está el binado de pt que vimos y utilizaremos por ahora; figure_compare_bkg.C, el que se ha utilizado para comparar MC de background solo de yy con los datos y para repesarlo para tener en cuenta los eventos yj y el de SpuriousSignal.C donde se calcula qué función aproxima mejor el background. Eso es todo por ahora.
