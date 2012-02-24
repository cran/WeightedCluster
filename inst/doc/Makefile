

RSCRIPT="$(R_HOME)/bin/Rscript"



all: pdfclean rotateFiles WeightedCluster.pdf WeightedCluster-fr.pdf WeightedClusterPreview.pdf clean

pdfclean:
	rm -f WeightedCluster-fr.pdf WeightedCluster.pdf WeightedClusterPreview.pdf


rotateFiles: pdfclean
	cp -f ./WeightedCluster/WeightedCluster-fake.Rnw WeightedCluster.Rnw
	cp -f ./WeightedCluster/WeightedCluster-fr-fake.Rnw WeightedCluster-fr.Rnw
	cp -f ./WeightedClusterPreview/WeightedClusterPreview-fake.Rnw WeightedClusterPreview.Rnw
	cp -f WeightedCluster/randB.RData .
	cp -f WeightedCluster/wardtree.png .
	cp -f WeightedCluster/seqtree.png .
	cp -f WeightedClusterPreview/previewhctree.png .
	

clean: WeightedCluster.pdf WeightedCluster-fr.pdf WeightedClusterPreview.pdf
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary
	rm -rf graphics
	rm -f randB.RData wardtree.png seqtree.png
	rm -f previewhctree.png
	
	
	
WeightedClusterPreview.pdf: rotateFiles WeightedClusterPreview/WeightedClusterPreview.Rnw
	cp -f WeightedClusterPreview/WeightedClusterPreview.Rnw .
	$(RSCRIPT) --vanilla -e "library(knitr); knit2pdf('WeightedClusterPreview.Rnw')"
	cp -f WeightedClusterPreview/WeightedClusterPreview-fake.Rnw WeightedClusterPreview.Rnw

WeightedCluster.pdf: rotateFiles WeightedCluster/WeightedCluster.Rnw
	cp -f WeightedCluster/WeightedCluster.Rnw .
	$(RSCRIPT) --vanilla -e "library(knitr); knit2pdf('WeightedCluster.Rnw')"
	cp -f WeightedCluster/WeightedCluster-fake.Rnw WeightedCluster.Rnw

WeightedCluster-fr.pdf: rotateFiles WeightedCluster/WeightedCluster-fr.Rnw
	cp -f WeightedCluster/WeightedCluster-fr.Rnw .
	$(RSCRIPT) --vanilla -e "library(knitr); knit2pdf('WeightedCluster-fr.Rnw')"
	cp -f WeightedCluster/WeightedCluster-fr-fake.Rnw WeightedCluster-fr.Rnw

	
compactPDF:  WeightedCluster-fr.pdf  WeightedCluster.pdf WeightedClusterPreview.pdf
	$(RSCRIPT) --vanilla -e "library(tools); compactPDF('./', gs_quality='ebook')"