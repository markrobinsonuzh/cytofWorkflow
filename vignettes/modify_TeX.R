
z <- readLines("cytofWorkflow.tex")
z <- paste(z, collapse="\n")
z <- gsub("\\\\textbf\\{Figure[ \n]?[0-9]{1,2}.}", "", z)

write.table(z, "main.tex", row.names = FALSE, col.names = FALSE, quote = FALSE)

## In LaTeX before building, comment out the first line, add the other two ..
#%\DefineVerbatimEnvironment{Highlighting}{Verbatim}{commandchars=\\\{\}}
#\DefineVerbatimEnvironment{Highlighting}{Verbatim}{commandchars=\\\{\},fontsize=\small}
#\setlength{\parskip}{.35em}

## comment out pagestyle
#%\pagestyle{main}

