# For Linux or any machines with gcc compiler
pandoc = /usr/local/bin/pandoc
Rscript = /usr/local/bin/Rscript

mainFile = Introduction
presFile = $(wildcard *.Rpres)
tempRmd = temp.Rmd
mdFile = $(mainFile).md
mdHeader = ../common/yamlHeader.Rmd

# make the updated script
all: printable

clean:
	/bin/rm $(tempRmd)
	/bin/rm printable.html

yaml: $(presFile)
	/bin/cat $(mdHeader) > $(tempRmd)
	/bin/cat $(presFile) >> $(tempRmd)
   
printable:  $(tempRmd)
	$(Rscript) --vanilla -e "library(rmarkdown);printed<-T;render('temp.Rmd','knitrBootstrap::bootstrap_document',output_file='01-Handout.html')"

rocco:  $presFile
	$(Rscript) --vanilla -e "library(knitr);printed<-T;rocco('Introduction.Rpres','rocco.html')"

pdfversion:  $mdFile
	$(pandoc) -o printable.pdf Segmentation.md
