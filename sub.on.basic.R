## Substitute many structures at once.

# start with a series of already substituted structures.
# Define the changed unit (different conformer essential unit) as a part
# of an xyz file (e.g. atoms 17-38 in group 1)
# replace all of the atoms in the substituted molecules.


dir <- 'structure00/'
coor.atoms <- '7 4 9'
keep.atoms <- '7 8 9 10 11 12 13 14 15 16 17'

sub.on.basic <- function(dir, coor.atoms, keep.atoms) {
  basic <- (list.files(path = dir, pattern = 'basic'))
  coor.trans.dir(dir, coor.atoms)
  setwd(dir)
  unlink(basic)
  file.rename(list.files(pattern = 'basic'), 'basic.xyz')
  molecules <- list.files(pattern = 'tc.xyz', full.names = F, recursive = F)
  atoms <- strsplit(keep.atoms, " ")
  unlisted.atoms <- unlist(atoms)
  numeric.atoms <- as.numeric(unlisted.atoms)
  repeat.unit <- data.table::fread(list.files(pattern = 'basic'))
  replace.atoms <- seq(1:nrow(repeat.unit))[-numeric.atoms]
  for (molecule in molecules) {
    xyz <- data.table::fread(molecule)
    xyz[replace.atoms, ] <- repeat.unit[replace.atoms, ]
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
    write(new_xyz, paste(stringr::str_replace(clean.name, '_tc', '_reoriented'), ".xyz", sep = ""))
    unlink(stringr::str_replace(molecule,'_tc',''))
    unlink(molecule)
  }
  name_changer('.','_reoriented','')
}

for (dir in list.files()) {
  sub.on.basic(dir, coor.atoms, keep.atoms)
}

sub.on.basic('structure00/', coor.atoms, keep.atoms)
sub.on.basic('two/', coor.atoms, keep.atoms)
sub.on.basic('three', coor.atoms, keep.atoms)
sub.on.basic('four', coor.atoms, keep.atoms)
sub.on.basic('five', coor.atoms, keep.atoms)
sub.on.basic('six', coor.atoms, keep.atoms)
sub.on.basic('seven', coor.atoms, keep.atoms)
sub.on.basic('eight', coor.atoms, keep.atoms)
sub.on.basic('nine', coor.atoms, keep.atoms)
sub.on.basic('ten', coor.atoms, keep.atoms)
sub.on.basic('eleven/', coor.atoms, keep.atoms)
sub.on.basic('twelve/', coor.atoms, keep.atoms)
sub.on.basic('thirteen/', coor.atoms, keep.atoms)
sub.on.basic('fourteen', coor.atoms, keep.atoms)
