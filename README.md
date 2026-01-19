# Articulo log-normal sesgado ictpc cdmx
Este repositorio contiene:

- Código fuente de R, Stan y conjunto de datos para reproducir el análisis estadístico. Para instalar Stan en R, siga las instrucciones en *https://mc-stan.org/cmdstanr/*. El conjunto de datos se dispone en formato `.csv` y también como objeto de R (`Xyregion-CDMX.RData`) para facilidad en el análisis.


- Código LaTeX para generar un `.pdf` con el documento.

Para generar un documento MS Word (`.docx`), puede emplearse `pandoc`. Sin embargo, no reproduce todos los elementos. Por ejemplo, el comando empleado es:


`pandoc Articulo-final.tex --filter pandoc-crossref --bibliography=referencias-short.bib --citeproc --csl apa.csl -o Articulo-final0.docx`

Un excelente tutorial para esta tarea: *https://youtu.be/28yfsJ1sj8I*


Saul Arturo Ortiz Muñoz (*ortiz.saul@colpos.mx*)
