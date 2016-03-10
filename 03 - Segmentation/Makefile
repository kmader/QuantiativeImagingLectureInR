# For Linux or any machines with gcc compiler
pandoc = /usr/local/bin/pandoc
Rscript = /usr/local/bin/Rscript

mainFile = Segmentation 
presFile = $mainFile.Rpres
mdFile = $mainFile.md

# make the updated script
all: printable pdfversion

clean:
	/bin/rm printable.html 
    
printable:  Segmentation.Rpres
	$(Rscript) --vanilla -e "library(knitrBootstrap);printed<-T;knit_bootstrap('Segmentation.Rpres',output='printable.html')"

rocco:  Segmentation.Rpres
	$(Rscript) --vanilla -e "library(knitr);printed<-T;rocco('Segmentation.Rpres','rocco.html')"

pdfversion:  Segmentation.md
	$(pandoc) -o printable.pdf Segmentation.md
