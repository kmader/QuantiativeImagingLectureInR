# For Linux or any machines with gcc compiler
pandoc = /usr/local/bin/pandoc
Rscript = /usr/bin/Rscript

mainFile = ShapeAnalysis
presFile = $(wildcard *.Rpres)
mdFile = $(mainFile).md

# make the updated script
all: printable

clean:
	/bin/rm printable.html 
    
printable:  $(presFile)
	$(Rscript) --vanilla -e "library(knitrBootstrap);printed<-T;knit_bootstrap('ShapeAnalysis.Rpres',output='printable.html',chooser=c('boot','code'))"

rocco:  $presFile
	$(Rscript) --vanilla -e "library(knitr);printed<-T;rocco('AdvSegmentation.Rpres','rocco.html')"

pdfversion:  $mdFile
	$(pandoc) -o printable.pdf Segmentation.md
