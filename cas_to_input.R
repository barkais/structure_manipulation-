# imported data from boronic acids table
install.packages('RSelenium')

BAs <- data.table::fread('boronic_acids_table.csv',header = T)
mols <- BAs[grep("^C", BAs$`*Molecular Formula*`), ]
name.mols <- BAs[!grep("^C", BAs$`*Molecular Formula*`), ]
new.bas <- cbind(name.mols[1:nrow(mols),], mols)
new.bas <- new.bas[,c(1,4,6)]
names(new.bas) <- c('Name', 'Formula', 'CAS')
new.bas$CAS[new.bas$CAS == ''] <- NA
new.bas <- dplyr::mutate(new.bas, smiles = rep(NA, nrow(new.bas)))
# Docker must be up and running - see https://cran.r-project.org/web/packages/RSelenium/vignettes/basics.html
# you should restart the docker webDriver before running the following 
remDr <- RSelenium::remoteDriver(
  port = 5900L,
  browser= "chrome",
)
remDr$open()
for (i in 1:nrow(new.bas)) {
  if (!is.na(new.bas$CAS[i])) {
    address <- paste("https://cactus.nci.nih.gov/chemical/structure/", new.bas$CAS[i], '/smiles', sep = '')
    remDr$navigate(address)
    smiles <- remDr$findElement('css selector', 'body')
    insert <- smiles$getElementText()
    if (insert == "Page not found (404)") {
      address <- paste("https://cactus.nci.nih.gov/chemical/structure/", new.bas$Name[i], '/smiles', sep = '')
      remDr$navigate(address)
      smiles <- remDr$findElement('css selector', 'body')
      insert <- smiles$getElementText()
    }
  } else {
    address <- paste("https://cactus.nci.nih.gov/chemical/structure/", new.bas$Name[i], '/smiles', sep = '')
    remDr$navigate(address)
    smiles <- remDr$findElement('css selector', 'body')
    insert <- smiles$getElementText()
  }
  new.bas$smiles[i] <- unlist(insert)
}

# Restart docker webserver

# needs.fixing <- new.bas[new.bas$smiles == "Page not found (404)",]
# needs.fixing <- needs.fixing[!is.na(needs.fixing$CAS),]
# remDr <- RSelenium::remoteDriver(
#   port = 5900L,
#   browser= "chrome",
# )
# suppressMessages(remDr$open())
# for (i in 1:nrow(needs.fixing)) {
#   address <- paste('https://pubchem.ncbi.nlm.nih.gov/#query=', needs.fixing$CAS[26], sep = '')
#   remDr$navigate(address)
#   iupac <- remDr$findElement('css selector', 
#                              '#featured-results > div > div.box-shadow > div > div:nth-child(1) > div:nth-child(2) > div:nth-child(4) > div > span > span.breakword > span')
#   iupac.replace <- iupac$getElementText()
#   address <- paste("https://cactus.nci.nih.gov/chemical/structure/", iupac.replace, '/smiles', sep = '')
#   remDr$navigate(address)
#   smiles <- remDr$findElement('css selector', 'body')
#   insert <- smiles$getElementText()
#   needs.fixing$smiles[i] <- insert
# }

new.bas <- new.bas[!grepl('\\+', new.bas$smiles), ] 
new.bas <- new.bas[new.bas$smiles != "Page not found (404)", ]
new.bas[1,1] <- 'basic'

new.bas.acids <- new.bas
for (i in 1:nrow(new.bas.acids)) {
  new.bas.acids$Name[i] <- stringr::str_replace_all(new.bas.acids$Name[i],"[[:space:][:punct:][:cntrl:]]", "")
  new.bas.acids$Name[i] <- stringr::str_replace_all(new.bas.acids$Name[i],"4", "p_")
  new.bas.acids$Name[i] <- stringr::str_replace_all(new.bas.acids$Name[i],"2", "o_")
  new.bas.acids$Name[i] <- stringr::str_replace_all(new.bas.acids$Name[i],"6", "o_")
  new.bas.acids$Name[i] <- stringr::str_replace_all(new.bas.acids$Name[i],"3", "m_")
  new.bas.acids$Name[i] <- stringr::str_replace_all(new.bas.acids$Name[i],"5", "m_")
  new.bas.acids$Name[i] <- stringr::str_replace_all(new.bas.acids$Name[i],"[[:digit:]]", "")
}

for (i in which(duplicated(new.bas.acids[,1]))) {
  new.bas.acids$Name[i] <- paste(new.bas.acids$Name[i], '_a', sep = '')
}
for (i in 1:nrow(new.bas.acids)) {
  write(as.character(new.bas.acids[i,4]),paste(as.character(new.bas.acids[i,1]), '.smi',sep = ''))
}

