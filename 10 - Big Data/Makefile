# For Linux or any machines with gcc compiler
pandoc = /usr/local/bin/pandoc
Rscript = /usr/bin/Rscript

mainFile = BigData
presFile = $(wildcard *.Rpres)
mdFile = $(mainFile).md

# make the updated script
all: printable

clean:
	/bin/rm 11-handout.html 
    
printable:  $(presFile)
	$(Rscript) --vanilla -e "library(knitrBootstrap);printed<-T;knit_bootstrap('BigData.Rpres',output='11-handout.html',chooser=c('boot','code'))"

rocco:  $presFile
	$(Rscript) --vanilla -e "library(knitr);printed<-T;rocco('BigData.Rpres','rocco.html')"

pdfversion:  $mdFile
	$(pandoc) -o printable.pdf Segmentation.md
