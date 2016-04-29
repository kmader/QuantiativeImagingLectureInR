# For Linux or any machines with gcc compiler
pandoc = /usr/local/bin/pandoc
Rscript = /usr/bin/Rscript

mainFile = Dynamic
presFile = $(wildcard *.Rpres)
mdFile = $(mainFile).md

# make the updated script
all: printable

clean:
	/bin/rm 10-handout.html
    
printable:  $(presFile)
	$(Rscript) --vanilla -e "library(knitrBootstrap);printed<-T;knit_bootstrap('Dynamic.Rpres',output='10-handout.html',chooser=c('boot','code'))"

rocco:  $presFile
	$(Rscript) --vanilla -e "library(knitr);printed<-T;rocco('Dynamic.Rpres','rocco.html')"

pdfversion:  $mdFile
	$(pandoc) -o printable.pdf Segmentation.md
