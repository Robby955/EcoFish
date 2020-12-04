This folder allows anyone to render and recreate my thesis. To create the thesis, you need the image files and all the latex files. Alternatively, you can render most of the files by running the R codes. The plots and saving the plot is commented out in the code, so you would need to remove the comment section and save the plot to the correct chapter.

This folder alone recreates the thesis. You create the file using the mainthesisUVIC.tex document. You need the image chapters saved in the same folder along with the frontmatter. 

If when running mainthesisUVIC.text creates the thesis with a missing table of contents, try rendering again. If it still does not render, you need to change the line in mainthesisUVIC.text from \bibliography{UvicBib} to \bibliography{UvicBib.bib}. 

If you are using a Mac and TeXShop,you will need to navigate to macros and select BibTex.

To recreate the thesis, you will need to have Latex installed and the correct packages updated.
