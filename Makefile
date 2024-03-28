PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)
BIOCVER := RELEASE_3_18

all: rd check clean

alldocs: rd readme mkdocs

crd:
	Rscript -e 'Rcpp::compileAttributes()'

rd: crd
	Rscript -e 'roxygen2::roxygenise(".")'

readme:
	Rscript -e 'rmarkdown::render("README.Rmd")'

readme2:
	Rscript -e 'rmarkdown::render("README.Rmd", "html_document")'

build:
	#cd ..;\
	#R CMD build $(PKGSRC)
	Rscript -e 'devtools::build()'
	
build2:
	cd ..;\
	R CMD build --no-build-vignettes $(PKGSRC)

install:
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: #build
	#cd ..;\
	#Rscript -e 'rcmdcheck::rcmdcheck("$(PKGNAME)_$(PKGVERS).tar.gz")'
	Rscript -e 'devtools::check()'

debug: rd build2 install

check2: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

bioccheck:
	cd ..;\
	Rscript -e 'BiocCheck::BiocCheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

rmrelease:
	git branch -D $(BIOCVER)

release:
	git checkout $(BIOCVER);\
	git fetch --all


clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

update:
	git fetch --all;\
	git checkout devel;\
	git merge upstream/devel;\
	git merge origin/devel


push:
	git push upstream devel;\
	git push origin devel


pages:
	Rscript -e 'rmarkdown::render("gh-pages/index.Rmd")'

publish:
	cd gh-pages;\
	git add .; git commit -m 'update'; git push





