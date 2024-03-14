# deploy.Makefile

# Target to deploy HTML files to gh-pages branch
deploy:
	mkdir -p _book_deploy
	cp -r _book/* _book_deploy/
	cd _book_deploy && \
	git init && \
	git checkout -b gh-pages || git checkout gh-pages && \
	git add --all && \
	git commit -m "Update HTML files" && \
	git push -f origin gh-pages
