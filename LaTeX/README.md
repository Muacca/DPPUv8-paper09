# Build the English manuscript

The inputs include the vector figure. From this directory run:

~~~sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
~~~

The third pass settles page and table-of-contents references if the contents page reflows on the second pass. Standard TeX Live packages suffice. The distributed PDF is at the repository root.
