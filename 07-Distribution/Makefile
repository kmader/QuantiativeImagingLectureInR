# For Linux or any machines with gcc compiler
pandoc = /usr/local/bin/pandoc
Rscript = /usr/local/bin/Rscript

mainFile = ComplexObjects
presFile = $(wildcard *.Rpres)
mdFile = $(mainFile).md

# make the updated script
all: printable

clean:
	/bin/rm 07-handout.html 
    
printable:  $(presFile)
	$(Rscript) --vanilla -e "library(knitrBootstrap);knit_bootstrap('DistributionAnalysis.Rpres',output='07a-handout.html')"

rocco:  $presFile
	$(Rscript) --vanilla -e "library(knitr);printed<-T;rocco('AdvSegmentation.Rpres','rocco.html')"

pdfversion:  $mdFile
	$(pandoc) -o printable.pdf Segmentation.md
