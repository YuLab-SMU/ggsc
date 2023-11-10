gitbook:
	Rscript -e 'library(bookdown); render_book("index.Rmd", "gitbook")'

bs4book:
	Rscript -e 'library(bookdown); render_book("index.Rmd", "bs4_book")'

pdfbook:
	Rscript -e 'library(bookdown); render_book("index.Rmd", "pdf_book")'

epub:
	Rscript -e 'library(bookdown); render_book("index.Rmd", "epub_book")'

clean:
	Rscript -e 'bookdown::clean_book()';\
	rm -rf _bookdown_files

serve:
	Rscript -e 'library(bookdown); serve_book()'

publish:
	cd gh-pages;\
	git add .;\
	git commit -m 'update';\
	git push
	