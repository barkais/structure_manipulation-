## Substitute many structures at once.

# start with a series of already substituted structures.
# Define the changed unit (different conformer essential unit) as a part
# of an xyz file (e.g. atoms 17-38 in group 1)
# replace all of the atoms in the substituted molecules.


dir <- 'Renumbering/'
coor.atoms <- '2 9 3'
keep.atoms <- 16:39
remove.atoms.from.sub <- 16

sub.with.common <- function(dir, coor.atoms, keep.atoms, remove.atoms.from.sub = NULL) {
  coor.trans.dir(dir, coor.atoms)
  setwd(dir)
  basic <- data.table::fread('basic_tc.xyz')
  unlink('basic.xyz')
  file.rename('basic_tc.xyz', 'basic.xyz')
  molecules <- list.files(pattern = 'tc.xyz', full.names = F, recursive = F)
  basic <- tibble::rowid_to_column(basic)
  attach.atoms <- basic[keep.atoms, ]
  atom.num.basic <- attach.atoms$rowid
  for (molecule in molecules) {
    xyz <- data.table::fread(molecule)
    if (!is.null(remove.atoms.from.sub)) xyz <- xyz[-remove.atoms.from.sub,]
    part.1 <- xyz[!(as.numeric(unlist(row.names(xyz))) %in% atom.num.basic),]
    part.2 <- setdiff(xyz, part.1)
    xyz <- rbind(part.1, attach.atoms[, 2:5], part.2)
    num.atoms <- nrow(xyz)
    m <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
    m[1, 1] <- num.atoms
    m[is.na(m)] <- ""
    names(m) <- names(xyz)
    xyz <- rbind(m, xyz)
    xyz[xyz == "0"] <- "0.0"
    new_xyz <- knitr::kable(xyz,
                            format = "simple", row.names = F, col.names = NULL
    )
    new_xyz <- new_xyz[-1]
    new_xyz <- new_xyz[-length(new_xyz)]
    clean.name <- tools::file_path_sans_ext(molecule)
    write(new_xyz, paste(stringr::str_replace(clean.name, '_tc', '_num'), ".xyz", sep = ""))
    unlink(stringr::str_replace(molecule,'_tc',''))
    unlink(molecule)
  }
  setwd('..')
}


