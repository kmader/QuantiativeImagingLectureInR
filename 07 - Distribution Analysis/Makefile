# For Linux or any machines with gcc compiler
pandoc = /usr/local/bin/pandoc
Rscript = /usr/bin/Rscript

mainFile = DistributionAnalysis
presFile = $(wildcard *.Rpres)
mdFile = $(mainFile).md

# make the updated script
all: printable

clean:
	/bin/rm printable.html 
    
printable:  $(presFile)
	$(Rscript) --vanilla -e "library(knitrBootstrap);printed<-T;knit_bootstrap('DistributionAnalysis.Rpres',output='08-Handout.html',chooser=c('boot','code'),markdown_options=c('mathjax', 'use_xhtml'))"

rocco:  $presFile
	$(Rscript) --vanilla -e "library(knitr);printed<-T;rocco('DistributionAnalysis.Rpres','rocco.html')"

pdfversion:  $mdFile
	$(pandoc) -o printable.pdf DistributionAnalysis.md
