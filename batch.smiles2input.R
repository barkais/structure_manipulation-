# generate multiple smiles files from a csv list

list.of.packages <- c("data.table", "stringr","tibble","dplyr","knitr",
                      "tools","pracma","qpcR","filesstrings","diversitree")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages)

smiles.from.csv <- function(csv.file) {
  dir.create('smiles')
  smiles.list <- data.table::fread(csv.file, header = F)
  for (i in 1:nrow(smiles.list)) {
   write(as.character(smiles.list[i,2]),paste(as.character(smiles.list[i,1]), '.smi',sep = ''))
  }
  filesstrings::file.move(list.files(pattern = ".smi"), "smiles", overwrite = T)
}

## After getting the smiles directory - transfer it to the cluster and follow the following steps
# 1. cd smiles
# 2. cp //gpfs0/gaus/users/barkais/scripts/smiles_to_structure.sh .
# 3. transfer back to local computer
# 4. In RStudio console do - setwd('structures_input_for_R/')
# 5. Run the script below 
#      |
#      |
#      |
#      |
#      v


batch.SMILES2input <- function(path) {
  options(scipen = 999)
  dir.create("../new_xyz_files", showWarnings = F)
  dir.create("../failed", showWarnings = F)
  molecules <- list.dirs(full.names = F, recursive = F)
  molecules <- molecules[!grepl("basic", molecules)]
  flag <- TRUE
  basic_atoms <- data.table::fread("basic/basic_atoms.csv")
  basic_bonds <- data.table::fread("basic/basic_bonds.csv", fill = T)
  basic_connectivity <- matrix(ncol = ncol(basic_bonds), nrow = nrow(basic_bonds))
  for (i in 1:dim(basic_bonds)[1]) {
    for (j in 1:dim(basic_bonds)[2]) {
      if (is.na(basic_bonds[i, ..j]) == FALSE) {
        basic_connectivity[i, j] <- basic_atoms$V2[basic_atoms$V1 == as.numeric(basic_bonds[i, ..j])]
      } else {
        basic_connectivity[i, j] <- NA
      }
    }
  }
  basic_ordered <- matrix(ncol = 1, nrow = nrow(basic_connectivity))
  for (i in 1:nrow(basic_connectivity)) {
    basic_ordered[i, ] <- stringr::str_remove_all(paste(sort(do.call(paste0, list(basic_connectivity[i, 1:5]))), collapse = ""), "NA")
  }
  answer <- readline('Would you like to choose the shared unit? T/F ')
  if (answer == 'T') {
    print(cbind(basic_ordered, basic_bonds))
    common.unit <<- as.numeric(readline('Please indicate the location of your chosen unit? '))
  }
  
  for (molecule in molecules) {
    tryCatch(
      expr = {
        setwd(molecule)
        # atomic indecies information of the basic structure from its pdb file
        basic_atoms <- data.table::fread("../basic/basic_atoms.csv")
        sub_atoms <- data.table::fread(list.files(pattern = "atoms.csv"))
        
        # connectivity information of the substitued structure from its pdb file
        basic_bonds <- data.table::fread("../basic/basic_bonds.csv", fill = T)
        sub_bonds <- data.table::fread(list.files(pattern = "bonds.csv"), fill = T)
        
        basic_connectivity <- matrix(ncol = ncol(basic_bonds), nrow = nrow(basic_bonds))
        for (i in 1:dim(basic_bonds)[1]) {
          for (j in 1:dim(basic_bonds)[2]) {
            if (is.na(basic_bonds[i, ..j]) == FALSE) {
              basic_connectivity[i, j] <- basic_atoms$V2[basic_atoms$V1 == as.numeric(basic_bonds[i, ..j])]
            } else {
              basic_connectivity[i, j] <- NA
            }
          }
        }
        
        
        sub_connectivity <- matrix(ncol = ncol(sub_bonds), nrow = nrow(sub_bonds))
        for (i in 1:dim(sub_bonds)[1]) {
          for (j in 1:dim(sub_bonds)[2]) {
            if (is.na(sub_bonds[i, ..j]) == FALSE) {
              sub_connectivity[i, j] <- sub_atoms$V2[sub_atoms$V1 == as.numeric(sub_bonds[i, ..j])]
            } else {
              sub_connectivity[i, j] <- NA
            }
          }
        }
        
        sub_ordered <- matrix(ncol = 1, nrow = nrow(sub_connectivity))
        for (i in 1:nrow(sub_connectivity)) {
          sub_ordered[i, ] <- stringr::str_remove_all(paste(sort(do.call(paste0, list(sub_connectivity[i, 1:5]))), collapse = ""), "NA")
        }
        
        basic_ordered <- matrix(ncol = 1, nrow = nrow(basic_connectivity))
        for (i in 1:nrow(basic_connectivity)) {
          basic_ordered[i, ] <- stringr::str_remove_all(paste(sort(do.call(paste0, list(basic_connectivity[i, 1:5]))), collapse = ""), "NA")
        }
        
        similarity.list <- list()
        for (i in 1:nrow(sub_ordered)) {
          similarity.list[[i]] <- which(basic_ordered[, 1] == sub_ordered[i, ])[1]
        }
        simi.indicator <- tibble::rowid_to_column(data.frame(matrix(similarity.list)))
        simi.indicator <- dplyr::filter(simi.indicator, matrix.similarity.list. != "NULL")
        simi.indicator <- dplyr::filter(simi.indicator, matrix.similarity.list. != "NA")
        simi.indicator <- simi.indicator[!rev(duplicated(rev(simi.indicator$matrix.similarity.list.))) , ]
        names(simi.indicator) <- c("sub", "basic")
        simi.indicator$basic <- unlist(simi.indicator$basic)
        basic_connectivity <- data.frame(basic_connectivity, stringsAsFactors = F)
        names(basic_connectivity) <- names(basic_bonds)
        sub_connectivity <- data.frame(sub_connectivity, stringsAsFactors = F)
        names(sub_connectivity) <- names(sub_bonds)
        list.basic <- list()
        list.sub <- list()
        for (i in 1:nrow(simi.indicator)) {
          list.basic[[i]] <- rbind(basic_connectivity[simi.indicator$basic[i], ], basic_bonds[simi.indicator$basic[i], ])
          list.sub[[i]] <- rbind(sub_connectivity[simi.indicator$sub[i], ], sub_bonds[simi.indicator$sub[i], ])
        }
        
        for (i in 1:length(list.basic)) {
          list.basic[[i]] <- data.frame(data.table::transpose(data.table::transpose(list.basic[[i]])[complete.cases(data.table::transpose(list.basic[[i]])), ]))
          names(list.basic[[i]]) <- as.character(list.basic[[i]][1, ])
          list.basic[[i]] <- list.basic[[i]][-1, ]
          dups <- which(duplicated(data.table::transpose(list.basic[[i]])))
          if (length(dups) > 0) {
            dups.first <- min(which(duplicated(data.table::transpose(list.basic[[i]])))) - 1
            names(list.basic[[i]])[dups] <- paste(names(list.basic[[i]])[dups], "_dup", sep = "")
            names(list.basic[[i]])[dups.first] <- paste(names(list.basic[[i]])[dups.first], "_dup", sep = "")
            list.basic[[i]] <- list.basic[[i]][, -dups]
          }
          names(list.basic[[i]]) <- make.unique(names(list.basic[[i]]))
          
          list.sub[[i]] <- data.frame(data.table::transpose(data.table::transpose(list.sub[[i]])[complete.cases(data.table::transpose(list.sub[[i]])), ]))
          names(list.sub[[i]]) <- as.character(list.sub[[i]][1, ])
          list.sub[[i]] <- list.sub[[i]][-1, ]
          dups <- which(duplicated(data.table::transpose(list.sub[[i]])))
          if (length(dups) > 0) {
            dups.first <- min(which(duplicated(data.table::transpose(list.sub[[i]])))) - 1
            names(list.sub[[i]])[dups] <- paste(names(list.sub[[i]])[dups], "_dup", sep = "")
            names(list.sub[[i]])[dups.first] <- paste(names(list.sub[[i]])[dups.first], "_dup", sep = "")
            list.sub[[i]] <- list.sub[[i]][, -dups]
          }
          names(list.sub[[i]]) <- make.unique(names(list.sub[[i]]))
          
          if (length(list.basic[[i]]) != length(list.sub[[i]])) {
            list.basic[[i]] <- NA
            list.sub[[i]] <- NA
          } else {
            data.table::setcolorder(list.sub[[i]], names(list.basic[[i]]))
          }
          Filter(Negate(anyNA), list.basic)
          Filter(Negate(anyNA), list.sub)
        }
        swap.dfs <- list()
        for (i in 1:length(list.basic)) {
          swap.dfs[[i]] <- rbind(list.sub[[i]], list.basic[[i]])
        }
        if (answer == 'T') {
          swap.list <- character()
          for (i in 1:length(swap.dfs)) {
            basic.unit <- data.frame(basic_bonds[common.unit])
            basic.unit <- basic.unit[ , colSums(is.na(basic.unit)) == 0]
            if (all(basic.unit[1, ] %in% swap.dfs[[i]][2, ])) {
              sub.com.unit <- swap.dfs[[i]]
              for (j in 1:ncol(sub.com.unit)) {
                swap.list <- paste(swap.list, sub.com.unit[1, j], sub.com.unit[2, j])
              }
            }
          }
        } else {
          hetero.ind <- vector()
          for (i in 1:length(swap.dfs)) {
            hetero.ind[[i]] <- all(grepl('C|H', colnames(swap.dfs[[i]])))
          }
          if (any(hetero.ind)) swap.dfs <- swap.dfs[which(!hetero.ind)]
          max.ncol <- which.max(lapply(swap.dfs, function(x) max(ncol(x))))
          swap.list <- character()
          for (i in 1:ncol(swap.dfs[[max.ncol]])) {
            swap.list <- paste(swap.list, swap.dfs[[max.ncol]][1, i], swap.dfs[[max.ncol]][2, i])
          }
        }
        
        swap.list <- trimws(swap.list)
        
        
        # swap shared unit numbering in the original substituted xyz
        co.sys <- data.table::fread(list.files(pattern = ".xyz"), header = F,
                                    colClasses = c("character", "numeric", "numeric", "numeric"))
        temp.co.sys <- data.frame(matrix(ncol = 4, nrow = nrow(co.sys)))
        swap.vec <- strsplit(swap.list, " ")
        unlisted.svec <- unlist(swap.vec)
        numeric.svec <- as.numeric(unlisted.svec)
        paired <- split(
          numeric.svec,
          ceiling(seq_along(numeric.svec) / 2)
        )
        replacers <- vector()
        for (i in 1:length(paired)) {
          temp.co.sys[paired[[i]][2], ] <- co.sys[paired[[i]][1], ]
          replacers <- append(replacers, paired[[i]][1])
        }
        swapped.atoms <- sort(replacers)
        missing <- unique(which(is.na(temp.co.sys), arr.ind = T)[, 1])
        temp.co.sys[missing, ] <- co.sys[-swapped.atoms, ]
        num.atoms <- nrow(temp.co.sys)
        m <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
        m[1, 1] <- num.atoms
        m[is.na(m)] <- ""
        names(m) <- names(temp.co.sys)
        temp.co.sys <- rbind(m, temp.co.sys)
        swapped <- knitr::kable(temp.co.sys,
                                format = "simple", row.names = F, col.names = NULL, align = "r"
        )
        swapped <- swapped[-1]
        swapped <- swapped[-length(swapped)]
        write(swapped, "temp.xyz")
        
        
        # transform coordinate systems to shared unit atoms
        max.col <- which.max(lapply(list.basic, function(x) sum(lengths(x))))
        if (answer == 'T') {
          tc.coordinates <- paste(as.character(sub.com.unit[2,1:3]), collapse = ' ')
        } else {
          tc.coordinates <- paste(as.character(swap.dfs[[max.ncol]][2,1:3]), collapse = ' ')  
        }
        
        coor.trans("temp.xyz", tc.coordinates)
        coor.trans("../basic/basic.xyz", tc.coordinates)
        basic.co.sys <- data.table::fread("../basic/basic_tc.xyz", header = F, colClasses = c("character", "numeric", "numeric", "numeric"))
        sub.co.sys <- data.table::fread("temp_tc.xyz", header = F, colClasses = c("character", "numeric", "numeric", "numeric"))
        
        mag <- function(vector) {
          sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
        }
        swap.align <- matrix(ncol = nrow(sub.co.sys), nrow = nrow(basic.co.sys))
        for (i in 1:nrow(swap.align)) {
          for (j in 1:ncol(swap.align)) {
            swap.align[i, j] <- mag(basic.co.sys[i, 2:4] - sub.co.sys[j, 2:4])
          }
        }
        swap.matrix <- data.frame(which(swap.align < 0.5, arr.ind = T))
        for (i in 1:nrow(swap.matrix)) {
          if (swap.matrix[i, 1] == swap.matrix[i, 2]) {
          }
        }
        
        final.swap.list <- character()
        for (i in 1:nrow(swap.matrix)) {
          final.swap.list <- paste(final.swap.list, swap.matrix[i, 2], swap.matrix[i, 1])
        }
        final.swap.list <- trimws(final.swap.list)
        final.co.sys <- data.frame(matrix(ncol = 4, nrow = nrow(sub.co.sys)))
        swap.vec <- strsplit(final.swap.list, " ")
        unlisted.svec <- unlist(swap.vec)
        numeric.svec <- as.numeric(unlisted.svec)
        paired <- split(
          numeric.svec,
          ceiling(seq_along(numeric.svec) / 2)
        )
        replacers <- vector()
        for (i in 1:length(paired)) {
          final.co.sys[paired[[i]][2], ] <- sub.co.sys[paired[[i]][1], ]
          replacers <- append(replacers, paired[[i]][1])
        }
        swapped.atoms <- sort(replacers)
        missing <- unique(which(is.na(final.co.sys), arr.ind = T)[, 1])
        final.co.sys[missing, ] <- sub.co.sys[-swapped.atoms, ]
        num.atoms <- nrow(final.co.sys)
        m <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
        m[1, 1] <- num.atoms
        m[is.na(m)] <- ""
        names(m) <- names(final.co.sys)
        final.co.sys <- rbind(m, final.co.sys)
        final.co.sys[final.co.sys == "0"] <- "0.0"
        swapped <- knitr::kable(final.co.sys,
                                format = "simple", row.names = F, col.names = NULL, align = "r"
        )
        swapped <- swapped[-1]
        swapped <- swapped[-length(swapped)]
        write(swapped, paste(stringr::str_remove(tools::file_path_sans_ext(molecule), "/"), "_fixed", ".xyz", sep = ""))
        
        
        ####### find missing matches #######
        fin.cs <- data.table::fread(paste(stringr::str_remove(tools::file_path_sans_ext(molecule), "/"), "_fixed", ".xyz", sep = ""),
                                    header = F, colClasses = c("character", "numeric", "numeric", "numeric")
        )
        replacies <- vector()
        for (i in 1:length(paired)) {
          replacies <- append(replacies, paired[[i]][2])
        }
        
        triang <- as.numeric(unlist(strsplit(tc.coordinates, " ")))
        
        basic.dis.mat <- data.frame(matrix(nrow = nrow(basic.co.sys), ncol = 4))
        basic.dis.mat[, 1] <- basic.co.sys[, 1]
        
        for (i in 1:nrow(basic.dis.mat)) {
          basic.dis.mat[i, 2] <- round(mag(basic.co.sys[i, 2:4]), digits = 3)
          basic.dis.mat[i, 3] <- round(mag(basic.co.sys[i, 2:4] - basic.co.sys[triang[2], 2:4]), digits = 3)
          basic.dis.mat[i, 4] <- round(mag(basic.co.sys[i, 2:4] - basic.co.sys[triang[3], 2:4]), digits = 3)
        }
        
        farthest <- as.matrix(basic.co.sys[which.max(basic.dis.mat$X2), 2:4])
        far.mag <- mag(farthest)
        
        basic.ang.mat <- data.frame(matrix(nrow = nrow(basic.dis.mat), ncol = 2))
        basic.ang.mat[, 1] <- basic.dis.mat[, 1]
        
        for (i in 1:nrow(basic.ang.mat)) {
          vec <- as.matrix(basic.co.sys[i, 2:4])
          vec.mag <- mag(vec)
          suppressWarnings(basic.ang.mat[i, 2] <- abs((pracma::dot(vec, farthest) / (vec.mag * far.mag)) * (180 / pi)))
        }
        
        sub.dis.mat <- data.frame(matrix(nrow = nrow(fin.cs), ncol = 4))
        sub.dis.mat[, 1] <- fin.cs[, 1]
        
        for (i in 1:nrow(sub.dis.mat)) {
          sub.dis.mat[i, 2] <- round(mag(fin.cs[i, 2:4]), digits = 3)
          sub.dis.mat[i, 3] <- round(mag(fin.cs[i, 2:4] - fin.cs[triang[2], 2:4]), digits = 3)
          sub.dis.mat[i, 4] <- round(mag(fin.cs[i, 2:4] - fin.cs[triang[3], 2:4]), digits = 3)
        }
        
        sub.ang.mat <- data.frame(matrix(nrow = nrow(sub.dis.mat), ncol = 2))
        sub.ang.mat[, 1] <- sub.dis.mat[, 1]
        
        for (i in 1:nrow(sub.ang.mat)) {
          vec <- as.matrix(fin.cs[i, 2:4])
          vec.mag <- mag(vec)
          suppressWarnings(sub.ang.mat[i, 2] <- abs((pracma::dot(vec, farthest) / (vec.mag * far.mag)) * (180 / pi)))
        }
        
        basic.dis.mat <- suppressWarnings(tibble::rowid_to_column(basic.dis.mat, "atom"))
        basic.dis.mat <- basic.dis.mat[-replacies, ]
        sub.dis.mat <- suppressWarnings(tibble::rowid_to_column(sub.dis.mat, "atom"))
        sub.dis.mat <- sub.dis.mat[-replacies, ]
        
        basic.ang.mat <- suppressWarnings(tibble::rowid_to_column(basic.ang.mat, "atom"))
        basic.ang.mat <- basic.ang.mat[-replacies, ]
        sub.ang.mat <- suppressWarnings(tibble::rowid_to_column(sub.ang.mat, "atom"))
        sub.ang.mat <- sub.ang.mat[-replacies, ]
        
        if (dim(basic.ang.mat)[1] != 0) {
          missing.pairs.ang <- matrix(nrow = nrow(basic.ang.mat), ncol = nrow(sub.ang.mat))
          for (i in 1:nrow(missing.pairs.ang)) {
            for (j in 1:ncol(missing.pairs.ang)) {
              missing.pairs.ang[i, j] <- abs(abs(basic.ang.mat[i, 3]) - abs(sub.ang.mat[j, 3]))
            }
          }
          row.names(missing.pairs.ang) <- row.names(basic.dis.mat)
          colnames(missing.pairs.ang) <- row.names(sub.dis.mat)
        }
        if (dim(basic.dis.mat)[1] != 0) {
          missing.pairs <- matrix(nrow = nrow(basic.dis.mat), ncol = nrow(sub.dis.mat))
          for (i in 1:nrow(missing.pairs)) {
            for (j in 1:ncol(missing.pairs)) {
              if (missing.pairs.ang[i, j] > 20) {
                missing.pairs[i, j] <- 100
              } else {
                missing.pairs[i, j] <- mag(basic.dis.mat[i, 3:5] - sub.dis.mat[j, 3:5])
              }
              if (missing.pairs[i, j] > 0.65) {
                missing.pairs[i, j] <- 100
              }
            }
          }
          for (i in 1:nrow(missing.pairs)) {
            for (j in 1:ncol(missing.pairs)) {
              if (missing.pairs[i, j] != min(missing.pairs[i,])) {
                missing.pairs[i, j] <- 100
              }
            }
          }
          row.names(missing.pairs) <- row.names(basic.dis.mat)
          colnames(missing.pairs) <- row.names(sub.dis.mat)
        }
        
        if (exists("missing.pairs") && all(missing.pairs == 100)) rm(missing.pairs)
        
        if (exists("missing.pairs")) {
          col.names <- colnames(missing.pairs)
          ids <- row.names(missing.pairs)
          out <- length(which(rowMeans(missing.pairs) == 100))
          missing.pairs <- missing.pairs[-which(rowMeans(missing.pairs) == 100), ]
          if (out != 0 && is.null(dim(missing.pairs))) {
            missing.pairs <- t(matrix(missing.pairs))
            colnames(missing.pairs) <- col.names
            rownames(missing.pairs) <- ids[1:nrow(missing.pairs)]
          }
          rogue.columns <- list()
          for (i in 1:ncol(missing.pairs)) {
            rogue.columns[i] <- length(which(missing.pairs[,i] < 100))
          }
          if (any(rogue.columns > 1)) {
            cols.values <- apply(missing.pairs, 2, unique)
            cols.values <- lapply(cols.values, length)
            bad.fruit <- which(missing.pairs[, cols.values > 2] != min(missing.pairs[, cols.values > 2]) && 
                                 missing.pairs[, cols.values > 2] != 100)
            missing.pairs <- missing.pairs[-bad.fruit,]
            if (is.null(dim(missing.pairs))) {
              missing.pairs <- t(matrix(missing.pairs))
              colnames(missing.pairs) <- col.names
              rownames(missing.pairs) <- ids[1:nrow(missing.pairs)]
            }
          }
          min.values <- apply(missing.pairs, 1, which.min)
          
          if (length(min.values) > 0) {
            for (i in 1:length(min.values)) {
              min.values[[i]] <- colnames(missing.pairs)[as.numeric(min.values[[i]])]
            }
            leave.alone <- names(min.values)[which(duplicated(min.values))]
            min.values[which(duplicated(min.values))] <- leave.alone
            new.swap.mat <- tibble::rownames_to_column(data.frame(min.values))
            new.swap.list <- character()
            for (i in 1:nrow(new.swap.mat)) {
              new.swap.list <- paste(new.swap.list, new.swap.mat[i, 2], new.swap.mat[i, 1])
            }
            new.swap.list <- trimws(new.swap.list)
            
            joined.co.sys <- qpcR:::cbind.na(basic.co.sys, fin.cs)
            dist.mat <- matrix(ncol = 1, nrow = nrow(joined.co.sys))
            for (i in 1:nrow(dist.mat)) {
              dist.mat[i, 1] <- mag(joined.co.sys[i, 2:4] - joined.co.sys[i, 6:8])
            }
            
            indexed <- which(dist.mat < 0.5, arr.ind = T)
            
            fixed.co.sys <- data.frame(matrix(ncol = 4, nrow = nrow(fin.cs)))
            swap.vec <- strsplit(new.swap.list, " ")
            unlisted.svec <- unlist(swap.vec)
            numeric.svec <- as.numeric(unlisted.svec)
            paired <- split(
              numeric.svec,
              ceiling(seq_along(numeric.svec) / 2)
            )
            
            indexed.T <- unique(sort(as.numeric(c(indexed[, 1], replacies))))
            
            if (length(paired) > 0) {
              for (i in 1:length(paired)) {
                if (paired[[i]][1] %in% indexed.T) next
                fixed.co.sys[paired[[i]][2], ] <- fin.cs[paired[[i]][1], ]
              }
            }
            fixed.co.sys[indexed.T, ] <- fin.cs[indexed.T, ]
            names(fixed.co.sys) <- names(fin.cs)
            repl <- suppressMessages(dplyr::setdiff(fin.cs, dplyr::semi_join(fixed.co.sys, fin.cs)))
            miss <- unique(which(is.na(fixed.co.sys), arr.ind = T)[, 1])
            fixed.co.sys[miss, ] <- repl
            num.atoms <- nrow(fixed.co.sys)
            m <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
            m[1, 1] <- num.atoms
            m[is.na(m)] <- ""
            names(m) <- names(fixed.co.sys)
            fixed.co.sys <- rbind(m, fixed.co.sys)
            fixed.co.sys[fixed.co.sys == "0"] <- "0.0"
            swapped <- knitr::kable(fixed.co.sys,
                                    format = "simple", row.names = F, col.names = NULL, align = "r"
            )
            swapped <- swapped[-1]
            swapped <- swapped[-length(swapped)]
            write(swapped, paste(stringr::str_remove(tools::file_path_sans_ext(molecule), "/"), "_fixed", ".xyz", sep = ""))
          }
          suppressWarnings(rm(missing.pairs, missing.pairs.ang))
        }
        unlink(list.files(pattern = "temp"))
        filesstrings::file.move(list.files(pattern = "fixed"), "../../new_xyz_files/", overwrite = T)
        setwd("..")
      }, error = function(e){
        flag <<- FALSE
        print(basename(getwd()))
        unlink(list.files(pattern = "temp"))
        filesstrings::file.move(list.files(pattern = "fixed"), "../../failed", overwrite = T)
        suppressWarnings(rm(missing.pairs, missing.pairs.ang))
        setwd('..')
      }
    ) 
    if (!flag) next()
  }
  setwd("basic")
  filesstrings::file.move("basic_tc.xyz", "../../new_xyz_files/", overwrite = T)
  setwd("../..")
  home <- getwd()
  setwd('new_xyz_files')
  name_changer("_fixed", "")
  name_changer("_tc", "")
  setwd(home)
}

diversitree::set.defaults(batch.SMILES2input, getwd())

# 6. Run the command (in RStudio console) - batch.SMILES2input()

# We now have a new folder with xyz files, numbered and ready for xTB computation 
# It is best to rename this folder before moving on 
# You might notice that there's also a folder names "Failed" - which in the better cases is empty and 
# is only there to catch unsucceful renumberings

# To continue - move the new xyz folder to the cluster 
# cd (get into) this folder on the clsuter and:
# cp //gpfs0/gaus/users/barkais/scripts/ .

