fns = list.files(root,
           "*-mirbase-ready.counts", recursive = T)[1:14]

coldata = data.frame(row.names = dirname(fns),
                     condition = substr(dirname(fns), 1, 2))

ids <- IsomirDataSeqFromFiles(file.path(root, fns),
                              coldata = coldata)
